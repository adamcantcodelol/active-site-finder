"""MBRC Active Site Finder."""
from __future__ import annotations

import base64
import html
import os

import pandas as pd
import requests
import streamlit as st
import streamlit.components.v1 as components

import active_site as asite
import foldseek_client as fs

st.set_page_config(page_title="MBRC Active Site Finder", page_icon="🧬", layout="wide")

st.markdown("""
<style>
.block-container{max-width:1400px;padding-top:1.4rem;padding-bottom:4rem}
.mbrc-header{display:flex;align-items:center;gap:12px;margin-bottom:.25rem}
.mbrc-m{font-size:2.05rem;font-weight:800;line-height:1}
.mbrc-logo{width:58px;height:58px;object-fit:contain;flex:0 0 58px}
.mbrc-title{font-size:1.8rem;font-weight:650;line-height:1.05;letter-spacing:-.03em}
.mbrc-subtitle{color:#64748b;font-size:1rem;margin:.35rem 0 1.5rem}
.metric-card{border:1px solid #dbe3ee;border-radius:12px;padding:16px 18px;background:#f8fafc;min-height:105px}
.metric-label{color:#64748b;font-size:.82rem;margin-bottom:5px}.metric-value{font-size:1.55rem;font-weight:700}.metric-note{color:#64748b;font-size:.78rem;margin-top:3px}
.site-callout{border:1px solid #cbd5e1;border-radius:12px;padding:18px;background:#f8fafc;margin:8px 0 16px}.site-title{font-size:1.15rem;font-weight:700;margin-bottom:4px}.site-residues{font-size:1.05rem;font-weight:650;word-break:break-word}.small-muted{color:#64748b;font-size:.84rem}
</style>
""", unsafe_allow_html=True)

base = os.path.dirname(__file__)
logo_html = ""
for path, mime in [
    (os.path.join(base, "assets", "mbrc_logo.svg"), "image/svg+xml"),
    (os.path.join(base, "assets", "mbrc_logo.png"), "image/png"),
]:
    if os.path.exists(path):
        with open(path, "rb") as f:
            logo_html = f'<img class="mbrc-logo" src="data:{mime};base64,{base64.b64encode(f.read()).decode()}">'
        break

st.markdown(
    '<div class="mbrc-header"><div class="mbrc-m">M</div>' + logo_html +
    '<div class="mbrc-title">Active Site Finder</div></div>'
    '<div class="mbrc-subtitle">Structure-first active-site prediction — ranked by the closest RMSD match.</div>',
    unsafe_allow_html=True,
)
st.markdown(
    "Enter a **PDB ID**. The app searches experimental PDB structures with Foldseek, "
    "ranks them by **RMSD (lower is better)**, and transfers known ligand-binding "
    "residues from the closest structural matches to your protein."
)

pdb_id = st.text_input("PDB ID", placeholder="e.g. 4HHB", max_chars=4).strip().upper()
run_clicked = st.button("Find active site", type="primary")


def render_viewer(pdb_text: str, residues: list[int], height: int = 720) -> None:
    active = ",".join(str(x) for x in sorted(set(residues)))
    pdb = pdb_text.replace("\\", "\\\\").replace("`", "\\`").replace("</script>", "<\\/script>")
    doc = """<!doctype html><html><head><meta charset="utf-8"><style>
html,body{margin:0;width:100%;height:100%;overflow:hidden;background:#fff}#viewer{width:100%;height:calc(100% - 46px)}#legend{height:46px;display:flex;align-items:center;padding:0 12px;color:#475569;font:13px Arial;border-top:1px solid #e2e8f0;box-sizing:border-box}
</style><script src="https://3Dmol.org/build/3Dmol-min.js"></script></head><body>
<div id="viewer"></div><div id="legend">Gray = protein &nbsp;&nbsp;|&nbsp;&nbsp; Orange = predicted active-site residues</div>
<script>
const v=$3Dmol.createViewer(document.getElementById('viewer'),{backgroundColor:'white'});
const m=v.addModel(`__PDB__`,'pdb');m.setStyle({},{cartoon:{color:'lightgray'}});const active='__ACTIVE__';
if(active){m.setStyle({resi:active},{cartoon:{color:'orange'},stick:{radius:.18,colorscheme:'orangeCarbon'}});m.addStyle({resi:active},{sphere:{radius:.35,colorscheme:'orangeCarbon'}});}
v.zoomTo();v.render();setTimeout(()=>v.resize(),100);window.addEventListener('resize',()=>v.resize());
</script></body></html>"""
    components.html(doc.replace("__PDB__", pdb).replace("__ACTIVE__", active), height=height, scrolling=False)


if run_clicked:
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        st.error("Enter a valid 4-character PDB ID, such as 4HHB.")
        st.stop()
    try:
        with st.spinner(f"Fetching {pdb_id} and finding structural matches..."):
            ticket, query_pdb = fs.submit_search_by_pdb_id(pdb_id, mode="tmalign", databases=["pdb100"])
        query_chain = asite.guess_first_chain_id(query_pdb)

        status = st.empty()
        fs.poll_until_complete(ticket, max_wait_seconds=300, on_status=lambda s,e: status.info(f"Foldseek search: {s.lower()} — {e}s elapsed"))
        status.empty()

        with st.spinner("Calculating RMSD for the structural matches..."):
            hits = fs.fetch_results(ticket, databases=["pdb100"])
            fs.populate_missing_rmsd(hits, query_pdb, query_chain_id=query_chain, max_hits=50)

        if not hits:
            st.warning("Foldseek returned no structural matches. Try another PDB entry.")
            st.stop()

        rmsd_hits = [h for h in hits if h.rmsd is not None]
        if rmsd_hits:
            best = rmsd_hits[0]
            st.subheader("Closest structural match")
            c1,c2,c3,c4=st.columns(4)
            with c1: st.markdown(f'<div class="metric-card"><div class="metric-label">PDB / chain</div><div class="metric-value">{html.escape(best.target_id)}</div><div class="metric-note">lowest RMSD match</div></div>',unsafe_allow_html=True)
            with c2: st.markdown(f'<div class="metric-card"><div class="metric-label">RMSD</div><div class="metric-value">{best.rmsd:.2f} Å</div><div class="metric-note">lower is closer</div></div>',unsafe_allow_html=True)
            with c3:
                tm = f"{best.tm_score:.3f}" if best.tm_score is not None else "n/a"
                st.markdown(f'<div class="metric-card"><div class="metric-label">TM-score</div><div class="metric-value">{tm}</div><div class="metric-note">structural similarity</div></div>',unsafe_allow_html=True)
            with c4:
                ident = f"{best.seq_identity*100:.1f}%" if best.seq_identity is not None else "n/a"
                st.markdown(f'<div class="metric-card"><div class="metric-label">Sequence identity</div><div class="metric-value">{ident}</div><div class="metric-note">aligned residues</div></div>',unsafe_allow_html=True)
            if best.description: st.caption(best.description)
        else:
            st.error("The search returned matches, but no valid RMSD could be calculated. The app will not pretend that a TM-score-ranked hit is an RMSD winner.")

        st.subheader("Structural matches — ranked by RMSD")
        shown_hits = rmsd_hits[:25] if rmsd_hits else hits[:25]
        st.dataframe(pd.DataFrame([{
            "Rank":i+1,"PDB / Chain":h.target_id,"Description":h.description or "n/a",
            "RMSD (Å)":f"{h.rmsd:.2f}" if h.rmsd is not None else "n/a",
            "TM-Score":f"{h.tm_score:.3f}" if h.tm_score is not None else "n/a",
            "Seq. Identity":f"{h.seq_identity*100:.1f}%" if h.seq_identity is not None else "n/a",
        } for i,h in enumerate(shown_hits)]),use_container_width=True,hide_index=True)

        st.subheader("Predicted active site")
        with st.spinner("Mapping experimentally observed ligand-binding residues..."):
            predicted = asite.predict_active_site(hits, query_pdb_text=query_pdb, query_chain_id=query_chain, top_n_hits=15)
        if not predicted:
            st.info("No active-site residues could be transferred from the experimental structural matches. The structural search itself completed successfully.")
        else:
            confident=[r for r in predicted if r.support_count>=2]
            site=confident if confident else predicted[:15]
            text=", ".join(f"{r.query_resname or '?'}{r.display_resnum}" for r in site)
            st.markdown(f'<div class="site-callout"><div class="site-title">Most strongly supported residues</div><div class="site-residues">{html.escape(text)}</div><div class="small-muted">Highest-supported residue has evidence from {site[0].support_count} structural template(s).</div></div>',unsafe_allow_html=True)
            st.dataframe(pd.DataFrame([{"Residue":f"{r.query_resname or '?'}{r.display_resnum}","Query position":r.display_resnum,"Homologs agreeing":r.support_count,"Supporting structures":", ".join(r.supporting_hits[:6])} for r in site]),use_container_width=True,hide_index=True)

        st.subheader("3D structure")
        render_viewer(query_pdb,[r.query_resnum for r in predicted])

    except fs.FoldseekError as exc:
        st.error(f"Structural search failed: {exc}")
    except requests.exceptions.RequestException as exc:
        st.error(f"Network request failed: {exc}")
    except Exception as exc:
        st.error(f"Unexpected error: {type(exc).__name__}: {exc}")
