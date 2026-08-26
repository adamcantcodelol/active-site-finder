"""MBRC Active Site Finder."""
from __future__ import annotations

import base64
import html
import os
import urllib.parse

import pandas as pd
import requests
import streamlit as st
import streamlit.components.v1 as components

import active_site as asite
import foldseek_client as fs

st.set_page_config(page_title="MBRC Active Site Finder", page_icon="🧬", layout="wide")

st.markdown("""
<style>
.block-container{max-width:1400px;padding-top:2.2rem;padding-bottom:4rem}
.mbrc-header{display:flex;align-items:center;gap:16px;margin-bottom:.3rem;min-height:86px;overflow:visible}
.mbrc-logo{width:230px;height:82px;object-fit:contain;object-position:center;display:block;flex:0 0 230px}
.mbrc-title{font-size:1.8rem;font-weight:650;line-height:1.05;letter-spacing:-.03em}
.mbrc-subtitle{color:#64748b;font-size:1rem;margin:.35rem 0 1.5rem}
.metric-card{border:1px solid #dbe3ee;border-radius:12px;padding:16px 18px;background:#f8fafc;min-height:105px}
.metric-label{color:#64748b;font-size:.82rem;margin-bottom:5px}.metric-value{font-size:1.55rem;font-weight:700}.metric-note{color:#64748b;font-size:.78rem;margin-top:3px}
.site-callout{border:1px solid #cbd5e1;border-radius:12px;padding:18px;background:#f8fafc;margin:8px 0 16px}.site-title{font-size:1.15rem;font-weight:700;margin-bottom:4px}.site-residues{font-size:1.05rem;font-weight:650;word-break:break-word}.small-muted{color:#64748b;font-size:.84rem}
.alignment-card{border:1px solid #dbe3ee;border-radius:12px;padding:18px;background:#fff;margin:8px 0 16px;overflow-x:auto}
.alignment-title{font-weight:700;font-size:1.05rem;margin-bottom:4px}.alignment-meta{color:#64748b;font-size:.84rem;margin-bottom:14px}
.aln-row{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:14px;line-height:1.8;white-space:pre}.aln-label{display:inline-block;width:86px;font-family:Arial,sans-serif;font-weight:700;color:#475569}.aln-match{color:#166534;font-weight:700}
.match-card{border:1px solid #dbe3ee;border-radius:12px;padding:16px 18px;background:#fff;margin:8px 0}.match-header{font-weight:700;font-size:1rem;margin-bottom:3px}.match-meta{color:#64748b;font-size:.82rem;margin-bottom:10px}.match-pattern{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:.92rem;line-height:1.75;word-break:break-word}.match-template{color:#475569}.match-query{font-weight:650}
.viewer-shell{border:1px solid #dbe3ee;border-radius:12px;overflow:hidden;background:white}
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
            logo_html = f'<img class="mbrc-logo" alt="MBRC logo" src="data:{mime};base64,{base64.b64encode(f.read()).decode()}">'
        break

st.markdown(
    '<div class="mbrc-header">' + logo_html +
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


def render_viewer(pdb_id: str, height: int = 760) -> None:
    safe_id = urllib.parse.quote(pdb_id.lower(), safe="")
    url = f"https://molstar.org/viewer/?pdb={safe_id}"
    st.markdown('<div class="viewer-shell">', unsafe_allow_html=True)
    components.iframe(url, height=height, scrolling=False)
    st.markdown("</div>", unsafe_allow_html=True)


def _best_unique_homologs(hits: list[fs.Hit], limit: int = 12) -> list[fs.Hit]:
    """Select exactly one best Foldseek chain for each unique PDB protein."""
    best_by_pdb: dict[str, fs.Hit] = {}
    for hit in hits:
        if not hit.q_aln or not hit.t_aln:
            continue
        old = best_by_pdb.get(hit.pdb_id)
        if old is None or (
            hit.rmsd is not None and (old.rmsd is None or hit.rmsd < old.rmsd)
        ):
            best_by_pdb[hit.pdb_id] = hit
    result = list(best_by_pdb.values())
    result.sort(key=lambda h: (
        h.rmsd is None,
        h.rmsd if h.rmsd is not None else float("inf"),
        -(h.tm_score if h.tm_score is not None else -1.0),
    ))
    return result[:limit]


def _site_correspondences(hit: fs.Hit, query_pdb: str, query_chain: str | None) -> list[dict]:
    """Map only the homolog's experimentally observed active-site residues."""
    target_pdb = asite.fetch_target_structure(hit)
    if not target_pdb or not hit.q_aln or not hit.t_aln:
        return []

    target_residues = asite.parse_ca_residues(target_pdb, hit.chain_id)
    query_residues = asite.parse_ca_residues(query_pdb, query_chain)
    site_residues = asite._template_site_residues(target_pdb, hit.chain_id)
    if not target_residues or not query_residues or not site_residues:
        return []

    site_keys = {(r["resnum"], r["insertion_code"]) for r in site_residues}
    ti, _ = fs._alignment_start_score(hit.t_aln, target_residues, hit.t_start)
    qi, _ = fs._alignment_start_score(hit.q_aln, query_residues, hit.q_start)
    out: list[dict] = []

    for qc, tc in zip(hit.q_aln, hit.t_aln):
        if qc != "-" and tc != "-" and 0 <= qi < len(query_residues) and 0 <= ti < len(target_residues):
            target = target_residues[ti]
            key = (target["resnum"], target["insertion_code"])
            if key in site_keys:
                query = query_residues[qi]
                out.append({
                    "template_chain": hit.chain_id or "?",
                    "template_resnum": target["resnum"],
                    "template_resname": target["resname"],
                    "query_chain": query_chain or "?",
                    "query_resnum": query["resnum"],
                    "query_resname": query["resname"],
                    "exact": target["resname"].upper() == query["resname"].upper(),
                })
        if qc != "-": qi += 1
        if tc != "-": ti += 1
    return out


def render_match_results(hits: list[fs.Hit], query_pdb: str, query_chain: str | None) -> None:
    homologs = _best_unique_homologs(hits, limit=12)
    if not homologs:
        st.info("No structural homolog alignments were returned.")
        return

    st.markdown(
        "**One best homolog per protein.** Each PDB protein is represented once, "
        "using its lowest-RMSD chain. The pattern below shows the experimentally "
        "defined binding-site residue in the homolog and the residue it corresponds "
        "to in your query, like SPRITE."
    )

    displayed = 0
    for hit in homologs:
        pairs = _site_correspondences(hit, query_pdb, query_chain)
        if not pairs:
            continue
        displayed += 1
        rmsd = f"{hit.rmsd:.2f} Å" if hit.rmsd is not None else "n/a"
        lines = []
        for p in pairs:
            lines.append(
                f'<div class="match-pattern">'
                f'<span class="match-template">{p["template_chain"]}{p["template_resnum"]} '
                f'{html.escape(p["template_resname"])}</span>'
                f' matches '
                f'<span class="match-query">{p["query_chain"]}{p["query_resnum"]} '
                f'{html.escape(p["query_resname"])}</span>'
                f'{" ✓" if p["exact"] else ""}'
                f'</div>'
            )
        st.markdown(
            f'<div class="match-card">'
            f'<div class="match-header">{html.escape(hit.pdb_id)} — '
            f'{html.escape(hit.description or "structural homolog")}</div>'
            f'<div class="match-meta">Chain {html.escape(hit.chain_id or "?")} · '
            f'RMSD {rmsd} · {len(pairs)} active-site residue match(es)</div>'
            + "".join(lines)
            + '</div>',
            unsafe_allow_html=True,
        )

    if displayed == 0:
        st.info("Homologs were found, but their deposited structures did not contain parseable SITE records or ligand-contact residues to display.")


if run_clicked:
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        st.error("Enter a valid 4-character PDB ID, such as 4HHB.")
        st.stop()

    try:
        with st.spinner(f"Fetching {pdb_id} and finding structural matches..."):
            ticket, query_pdb = fs.submit_search_by_pdb_id(
                pdb_id, mode="tmalign", databases=["pdb100"]
            )
        query_chain = asite.guess_first_chain_id(query_pdb)

        status = st.empty()
        fs.poll_until_complete(
            ticket,
            max_wait_seconds=300,
            on_status=lambda s, e: status.info(
                f"Foldseek search: {s.lower()} — {e}s elapsed"
            ),
        )
        status.empty()

        with st.spinner("Calculating RMSD for the structural matches..."):
            hits = fs.fetch_results(ticket, databases=["pdb100"])
            fs.populate_missing_rmsd(
                hits, query_pdb, query_chain_id=query_chain, max_hits=50
            )

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
                    f'<div class="metric-note">lowest RMSD match</div></div>',
                    unsafe_allow_html=True,
                )
            with c2:
                st.markdown(
                    f'<div class="metric-card"><div class="metric-label">RMSD</div>'
                    f'<div class="metric-value">{best.rmsd:.2f} Å</div>'
                    f'<div class="metric-note">lower is closer</div></div>',
                    unsafe_allow_html=True,
                )
            with c3:
                tm = f"{best.tm_score:.3f}" if best.tm_score is not None else "n/a"
                st.markdown(
                    f'<div class="metric-card"><div class="metric-label">TM-score</div>'
                    f'<div class="metric-value">{tm}</div>'
                    f'<div class="metric-note">structural similarity</div></div>',
                    unsafe_allow_html=True,
                )
            with c4:
                ident = f"{best.seq_identity * 100:.1f}%" if best.seq_identity is not None else "n/a"
                st.markdown(
                    f'<div class="metric-card"><div class="metric-label">Sequence identity</div>'
                    f'<div class="metric-value">{ident}</div>'
                    f'<div class="metric-note">aligned residues</div></div>',
                    unsafe_allow_html=True,
                )
            if best.description:
                st.caption(best.description)
        else:
            st.error(
                "The search returned matches, but no valid RMSD could be calculated. "
                "The app will not pretend that a TM-score-ranked hit is an RMSD winner."
            )
            best = None

        st.subheader("Structural matches — ranked by RMSD")
        shown_hits = rmsd_hits[:25] if rmsd_hits else hits[:25]
        st.dataframe(
            pd.DataFrame([
                {
                    "Rank": i + 1,
                    "PDB / Chain": h.target_id,
                    "Description": h.description or "n/a",
                    "RMSD (Å)": f"{h.rmsd:.2f}" if h.rmsd is not None else "n/a",
                    "TM-Score": f"{h.tm_score:.3f}" if h.tm_score is not None else "n/a",
                    "Seq. Identity": f"{h.seq_identity * 100:.1f}%" if h.seq_identity is not None else "n/a",
                }
                for i, h in enumerate(shown_hits)
            ]),
            use_container_width=True,
            hide_index=True,
        )

        st.subheader("Active-site match results")
        render_match_results(hits, query_pdb, query_chain)

        st.subheader("Predicted active site")
        with st.spinner("Mapping experimentally observed ligand-binding residues..."):
            predicted = asite.predict_active_site(
                hits,
                query_pdb_text=query_pdb,
                query_chain_id=query_chain,
                top_n_hits=15,
            )
        if not predicted:
            st.info(
                "No active-site residues could be transferred from the experimental "
                "structural matches. The structural search itself completed successfully."
            )
        else:
            confident = [r for r in predicted if r.support_count >= 2]
            site = confident if confident else predicted[:15]
            text = ", ".join(f"{r.query_resname or '?'}{r.display_resnum}" for r in site)
            st.markdown(
                f'<div class="site-callout"><div class="site-title">Most strongly supported residues</div>'
                f'<div class="site-residues">{html.escape(text)}</div>'
                f'<div class="small-muted">Highest-supported residue has evidence from '
                f'{site[0].support_count} structural template(s).</div></div>',
                unsafe_allow_html=True,
            )
            st.dataframe(
                pd.DataFrame([
                    {
                        "Residue": f"{r.query_resname or '?'}{r.display_resnum}",
                        "Query position": r.display_resnum,
                        "Homologs agreeing": r.support_count,
                        "Supporting structures": ", ".join(r.supporting_hits[:6]),
                    }
                    for r in site
                ]),
                use_container_width=True,
                hide_index=True,
            )

        st.subheader("3D structure")
        render_viewer(pdb_id)

    except fs.FoldseekError as exc:
        st.error(f"Structural search failed: {exc}")
    except requests.exceptions.RequestException as exc:
        st.error(f"Network request failed: {exc}")
    except Exception as exc:
        st.error(f"Unexpected error: {type(exc).__name__}: {exc}")
