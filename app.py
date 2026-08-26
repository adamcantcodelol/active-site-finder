"""
MBRC Active Site Finder

Primary goal:
    Find the closest structural matches to a query PDB structure by RMSD,
    then transfer experimentally observed ligand-binding residues from the
    closest matching PDB structures onto the query.
"""
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


st.set_page_config(
    page_title="MBRC Active Site Finder",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
<style>
.block-container{max-width:1400px;padding-top:1.5rem;padding-bottom:4rem}
.mbrc-header{display:flex;align-items:center;gap:10px;margin-bottom:.25rem}
.mbrc-logo{width:64px;height:64px;object-fit:contain;flex:0 0 64px}
.mbrc-brand-m{font-size:2.05rem;font-weight:800;line-height:1.05;letter-spacing:-.04em}
.mbrc-title{font-size:1.8rem;font-weight:650;line-height:1.05;letter-spacing:-.03em;margin-left:0}
.mbrc-subtitle{color:#64748b;font-size:1rem;margin:.35rem 0 1.5rem 74px}
.metric-card{border:1px solid #dbe3ee;border-radius:12px;padding:16px 18px;background:#f8fafc;min-height:105px}
.metric-label{color:#64748b;font-size:.82rem;margin-bottom:5px}.metric-value{font-size:1.55rem;font-weight:700}.metric-note{color:#64748b;font-size:.78rem;margin-top:3px}
.site-callout{border:1px solid #cbd5e1;border-radius:12px;padding:18px;background:#f8fafc;margin:8px 0 16px}.site-title{font-size:1.15rem;font-weight:700;margin-bottom:4px}.site-residues{font-size:1.05rem;font-weight:650;word-break:break-word}.small-muted{color:#64748b;font-size:.84rem}
.viewer-shell{border:1px solid #e2e8f0;border-radius:12px;overflow:hidden;background:#fff}
</style>
""",
    unsafe_allow_html=True,
)

base_dir = os.path.dirname(__file__)
logo_candidates = [
    (os.path.join(base_dir, "assets", "mbrc_logo.svg"), "image/svg+xml"),
    (os.path.join(base_dir, "assets", "mbrc_logo.png"), "image/png"),
    (os.path.join(base_dir, "mbrc_logo.png"), "image/png"),
]
logo_html = ""
for logo_path, mime in logo_candidates:
    if os.path.exists(logo_path):
        with open(logo_path, "rb") as f:
            logo_b64 = base64.b64encode(f.read()).decode("ascii")
        logo_html = f'<img class="mbrc-logo" src="data:{mime};base64,{logo_b64}" alt="BRC shield logo">'
        break
if not logo_html:
    logo_html = '<div class="mbrc-logo" aria-hidden="true"></div>'

header_html = (
    '<div class="mbrc-header">'
    + logo_html
    + '<div class="mbrc-brand-m">M</div>'
    + '<div class="mbrc-title">Active Site Finder</div>'
    + '</div>'
    + '<div class="mbrc-subtitle">Structure-first active-site prediction — ranked by the closest RMSD match.</div>'
)
st.markdown(header_html, unsafe_allow_html=True)

st.markdown(
    "Enter a **PDB ID**. The app searches experimental PDB structures with Foldseek, "
    "ranks them by **RMSD (lower is better)**, and transfers known ligand-binding "
    "residues from the closest structural matches to your protein."
)

pdb_id = st.text_input("PDB ID", placeholder="e.g. 4HHB", max_chars=4).strip().upper()
run_clicked = st.button("Find active site", type="primary")


def render_viewer(query_pdb: str, active_site_resnums: list[int], height: int = 720) -> None:
    """Render 3Dmol without putting JavaScript braces inside a Python f-string."""
    def escape_js_template(text: str) -> str:
        return (
            text.replace("\\", "\\\\")
            .replace("`", "\\`")
            .replace("</script>", "<\\/script>")
        )

    resi_selector = ",".join(str(n) for n in sorted(set(active_site_resnums)))
    pdb_js = escape_js_template(query_pdb)

    # Build the HTML with placeholders. This avoids Python f-string brace
    # parsing entirely, which was the source of the Streamlit SyntaxError.
    html_doc = """
<!doctype html>
<html>
<head>
<meta charset="utf-8">
<style>
html,body{margin:0;padding:0;width:100%;height:100%;overflow:hidden;background:#fff;font-family:Arial,sans-serif}
#viewer3d{width:100%;height:calc(100% - 46px);position:relative}
#legend{height:46px;display:flex;align-items:center;padding:0 12px;color:#475569;font-size:13px;box-sizing:border-box;border-top:1px solid #e2e8f0;background:#fff}
</style>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
</head>
<body>
<div id="viewer3d"></div>
<div id="legend">Gray = protein &nbsp;&nbsp;|&nbsp;&nbsp; Orange = predicted active-site residues</div>
<script>
const el=document.getElementById("viewer3d");
const viewer=$3Dmol.createViewer(el,{backgroundColor:"white"});
const model=viewer.addModel(`__PDB__`,"pdb");
model.setStyle({}, {cartoon:{color:"lightgray"}});
const active="__ACTIVE__";
if(active){
  model.setStyle({resi:active}, {cartoon:{color:"orange"}, stick:{radius:0.18,colorscheme:"orangeCarbon"}});
  model.addStyle({resi:active}, {sphere:{radius:0.35,colorscheme:"orangeCarbon"}});
}
viewer.zoomTo();
viewer.render();
setTimeout(()=>viewer.resize(),100);
window.addEventListener("resize",()=>viewer.resize());
</script>
</body>
</html>
"""
    html_doc = html_doc.replace("__PDB__", pdb_js).replace("__ACTIVE__", resi_selector)
    components.html(html_doc, height=height, scrolling=False)


if run_clicked:
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        st.error("Enter a valid 4-character PDB ID, such as 4HHB.")
        st.stop()

    try:
        with st.spinner(f"Fetching {pdb_id} and finding the closest structural matches..."):
            ticket_id, query_pdb_text = fs.submit_search_by_pdb_id(
                pdb_id, mode="tmalign", databases=["pdb100"]
            )

        status_box = st.empty()
        fs.poll_until_complete(
            ticket_id,
            max_wait_seconds=300,
            on_status=lambda status, elapsed: status_box.info(
                f"Foldseek search: {status.lower()} — {elapsed}s elapsed"
            ),
        )
        status_box.empty()

        with st.spinner("Reading structural matches and calculating the RMSD ranking..."):
            hits = fs.fetch_results(ticket_id, databases=["pdb100"])
            fs.populate_missing_rmsd(hits, query_pdb_text, max_hits=50)

        if not hits:
            st.warning("Foldseek returned no structural matches. Try another PDB entry.")
            st.stop()

        rmsd_hits = [h for h in hits if h.rmsd is not None]
        if rmsd_hits:
            best = rmsd_hits[0]
            st.subheader("Closest structural match")
            c1, c2, c3, c4 = st.columns(4)
            with c1:
                st.markdown(
                    f'<div class="metric-card"><div class="metric-label">PDB / chain</div>'
                    f'<div class="metric-value">{html.escape(best.target_id)}</div>'
                    '<div class="metric-note">lowest RMSD match</div></div>',
                    unsafe_allow_html=True,
                )
            with c2:
                st.markdown(
                    f'<div class="metric-card"><div class="metric-label">RMSD</div>'
                    f'<div class="metric-value">{best.rmsd:.2f} Å</div>'
                    '<div class="metric-note">lower is closer</div></div>',
                    unsafe_allow_html=True,
                )
            with c3:
                tm = f"{best.tm_score:.3f}" if best.tm_score is not None else "n/a"
                st.markdown(
                    f'<div class="metric-card"><div class="metric-label">TM-score</div>'
                    f'<div class="metric-value">{tm}</div>'
                    '<div class="metric-note">structural similarity</div></div>',
                    unsafe_allow_html=True,
                )
            with c4:
                ident = f"{best.seq_identity * 100:.1f}%" if best.seq_identity is not None else "n/a"
                st.markdown(
                    f'<div class="metric-card"><div class="metric-label">Sequence identity</div>'
                    f'<div class="metric-value">{ident}</div>'
                    '<div class="metric-note">aligned residues</div></div>',
                    unsafe_allow_html=True,
                )
            if best.description:
                st.caption(best.description)
        else:
            st.warning(
                "Foldseek returned matches, but RMSD could not be calculated from the available "
                "alignments. No RMSD-based winner is shown."
            )

        st.subheader("Structural matches — ranked by RMSD")
        table_hits = rmsd_hits[:25] if rmsd_hits else hits[:25]
        table_df = pd.DataFrame(
            [
                {
                    "Rank": i + 1,
                    "PDB / Chain": h.target_id,
                    "Description": h.description or "n/a",
                    "RMSD (Å)": f"{h.rmsd:.2f}" if h.rmsd is not None else "n/a",
                    "TM-Score": f"{h.tm_score:.3f}" if h.tm_score is not None else "n/a",
                    "Seq. Identity": f"{h.seq_identity * 100:.1f}%" if h.seq_identity is not None else "n/a",
                }
                for i, h in enumerate(table_hits)
            ]
        )
        st.dataframe(table_df, use_container_width=True, hide_index=True)

        st.subheader("Predicted active site")
        query_chain = asite.guess_first_chain_id(query_pdb_text)
        with st.spinner(
            "Transferring known ligand-binding sites from the closest-RMSD experimental structures..."
        ):
            predicted = asite.predict_active_site(
                hits,
                query_pdb_text=query_pdb_text,
                query_chain_id=query_chain,
                top_n_hits=15,
            )

        if not predicted:
            st.info(
                "No active-site residues could be transferred from the structural templates. "
                "The search completed successfully, but the closest-RMSD structures may lack "
                "deposited ligand annotations or usable residue alignments."
            )
        else:
            confident = [r for r in predicted if r.support_count >= 2]
            shown = confident if confident else predicted[:15]
            residue_text = ", ".join(
                f"{r.query_resname or '?'}{r.display_resnum}" for r in shown
            )
            strongest = shown[0]
            st.markdown(
                f'<div class="site-callout"><div class="site-title">Most strongly supported residues</div>'
                f'<div class="site-residues">{html.escape(residue_text)}</div>'
                f'<div class="small-muted">Highest-supported residue has evidence from '
                f'{strongest.support_count} structural template(s).</div></div>',
                unsafe_allow_html=True,
            )
            site_df = pd.DataFrame(
                [
                    {
                        "Residue": f"{r.query_resname or '?'}{r.display_resnum}",
                        "Query position": r.display_resnum,
                        "Homologs agreeing": r.support_count,
                        "Supporting structures": ", ".join(r.supporting_hits[:6])
                        + (" ..." if len(r.supporting_hits) > 6 else ""),
                    }
                    for r in shown
                ]
            )
            st.dataframe(site_df, use_container_width=True, hide_index=True)
            st.caption(
                "Multiple independent experimental structures agreeing on the same query residue "
                "provide stronger support."
            )

        st.subheader("3D structure")
        render_viewer(query_pdb_text, [r.query_resnum for r in predicted], height=720)

    except fs.FoldseekError as exc:
        st.error(f"Structural search failed: {exc}")
    except requests.exceptions.RequestException as exc:
        st.error(f"Network request failed: {exc}")
    except Exception as exc:
        st.error(f"Unexpected error: {type(exc).__name__}: {exc}")
