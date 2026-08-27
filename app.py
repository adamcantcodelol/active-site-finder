"""MBRC Active Site Finder - compact SPRITE-style active-site matching."""
from __future__ import annotations

import html
import urllib.parse

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

import active_site as asite
import foldseek_client as fs

st.set_page_config(page_title="MBRC Active Site Finder", page_icon="🧬", layout="wide")

st.markdown("""
<style>
.block-container{max-width:1400px;padding-top:2.7rem;padding-bottom:4rem}
.mbrc-header{display:flex;align-items:center;gap:14px;min-height:94px;overflow:visible}
.mbrc-logo{width:226px;height:90px;display:block;overflow:visible;flex:none}
.mbrc-title{font-size:1.8rem;font-weight:650;line-height:1.05;letter-spacing:-.03em}
.mbrc-subtitle{color:#64748b;font-size:1rem;margin:.25rem 0 1.5rem}
.metric-card{border:1px solid #dbe3ee;border-radius:12px;padding:16px 18px;background:#f8fafc;min-height:105px}
.metric-label{color:#64748b;font-size:.82rem;margin-bottom:5px}.metric-value{font-size:1.55rem;font-weight:700}.metric-note{color:#64748b;font-size:.78rem;margin-top:3px}
.sprite-card{border:1px solid #dbe3ee;border-radius:12px;background:#fff;padding:18px 20px;margin:8px 0 20px;overflow-x:auto}
.sprite-head{font-size:1.1rem;font-weight:700;margin-bottom:3px}.sprite-meta{color:#64748b;font-size:.84rem;margin-bottom:14px}
.sprite-table{border-collapse:collapse;width:100%;min-width:560px;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
.sprite-table th{text-align:left;color:#64748b;font-family:Arial,Helvetica,sans-serif;font-size:.78rem;font-weight:600;border-bottom:1px solid #e2e8f0;padding:7px 9px}
.sprite-table td{padding:9px;border-bottom:1px solid #eef2f7;font-size:.93rem}.sprite-table tr:last-child td{border-bottom:0}
.sprite-exact{font-weight:700}.sprite-arrow{color:#64748b;padding:0 10px}.sprite-note{font-family:Arial,Helvetica,sans-serif;color:#64748b;font-size:.8rem;margin-top:10px}
.site-callout{border:1px solid #cbd5e1;border-radius:12px;padding:18px;background:#f8fafc;margin:8px 0 16px}.site-title{font-size:1.15rem;font-weight:700}.site-residues{font-size:1.05rem;font-weight:650;word-break:break-word}.small-muted{color:#64748b;font-size:.84rem}
</style>
""", unsafe_allow_html=True)

# M is inside the shield; BRC is outside. The extra SVG viewBox height prevents clipping.
logo = """<svg class='mbrc-logo' viewBox='0 0 470 130' xmlns='http://www.w3.org/2000/svg' role='img' aria-label='MBRC logo'>
<path d='M92 10 L162 28 L162 60 C162 91 140 112 92 124 C44 112 22 91 22 60 L22 28 Z' fill='none' stroke='#111827' stroke-width='5' stroke-linejoin='round'/>
<text x='92' y='82' text-anchor='middle' font-family='Arial,Helvetica,sans-serif' font-size='58' font-weight='800' fill='#111827'>M</text>
<text x='204' y='82' font-family='Arial,Helvetica,sans-serif' font-size='54' font-weight='800' letter-spacing='-3' fill='#111827'>BRC</text>
</svg>"""
st.markdown(f'<div class="mbrc-header">{logo}<div class="mbrc-title">Active Site Finder</div></div><div class="mbrc-subtitle">Structure-first active-site prediction — ranked by the closest RMSD match.</div>', unsafe_allow_html=True)
st.markdown("Enter a **PDB ID**. The app searches experimental PDB structures with Foldseek, ranks them by **RMSD (lower is better)**, and transfers known ligand-binding residues from structural matches to your protein.")

pdb_id = st.text_input("PDB ID", placeholder="e.g. 4HHB", max_chars=4).strip().upper()
run_clicked = st.button("Find active site", type="primary")


def _best_hit(hits: list[fs.Hit]) -> fs.Hit | None:
    valid = [h for h in hits if h.rmsd is not None and h.q_aln and h.t_aln]
    return min(valid, key=lambda h: h.rmsd) if valid else None


def _map_site_pairs(hit: fs.Hit, query_pdb: str, query_chain: str | None) -> list[dict]:
    """Map the experimentally annotated site of one template to the query."""
    target_pdb = asite.fetch_target_structure(hit)
    if not target_pdb or not hit.q_aln or not hit.t_aln:
        return []
    target = asite.parse_ca_residues(target_pdb, hit.chain_id)
    query = asite.parse_ca_residues(query_pdb, query_chain)
    site = asite._template_site_residues(target_pdb, hit.chain_id)
    if not target or not query or not site:
        return []
    site_keys = {(r["resnum"], r["insertion_code"]) for r in site}
    ti, _ = fs._alignment_start_score(hit.t_aln, target, hit.t_start)
    qi, _ = fs._alignment_start_score(hit.q_aln, query, hit.q_start)
    pairs = []
    for qc, tc in zip(hit.q_aln, hit.t_aln):
        if qc != "-" and tc != "-" and 0 <= qi < len(query) and 0 <= ti < len(target):
            tr = target[ti]
            if (tr["resnum"], tr["insertion_code"]) in site_keys:
                qr = query[qi]
                pairs.append({
                    "template_chain": hit.chain_id or "?",
                    "template_resnum": tr["resnum"],
                    "template_resname": tr["resname"],
                    "query_chain": query_chain or "?",
                    "query_resnum": qr["resnum"],
                    "query_resname": qr["resname"],
                    "exact": tr["resname"].upper() == qr["resname"].upper(),
                })
        if qc != "-": qi += 1
        if tc != "-": ti += 1
    return pairs


def _choose_local_triplet(pairs: list[dict]) -> list[dict]:
    """Choose one tight 3-residue local cluster, rather than the whole site."""
    pairs = sorted(pairs, key=lambda p: (p["query_resnum"], p["template_resnum"]))
    if len(pairs) <= 3:
        return pairs
    best = None
    for i in range(len(pairs) - 2):
        group = pairs[i:i + 3]
        span = group[-1]["query_resnum"] - group[0]["query_resnum"]
        exact_count = sum(1 for p in group if p["exact"])
        score = (span, -exact_count, group[0]["query_resnum"])
        if best is None or score < best[0]:
            best = (score, group)
    return best[1] if best else pairs[:3]


def _find_best_site_hit(hits: list[fs.Hit], query_pdb: str, query_chain: str | None) -> tuple[fs.Hit | None, list[dict]]:
    """Find the closest-RMSD homolog that actually has a mappable site.

    A structural hit can be the absolute RMSD winner while having no ligand or
    SITE annotation. That entry cannot produce a SPRITE-style site match, so
    we select the closest usable experimental PDB instead of failing the UI.
    """
    candidates = sorted([h for h in hits if h.rmsd is not None and h.q_aln and h.t_aln], key=lambda h: h.rmsd)
    for hit in candidates[:30]:
        pairs = _map_site_pairs(hit, query_pdb, query_chain)
        if pairs:
            return hit, pairs
    return None, []


def _render_sprite(hit: fs.Hit, pairs: list[dict]) -> None:
    chosen = _choose_local_triplet(pairs)
    rows = []
    for p in chosen:
        check = " ✓" if p["exact"] else ""
        rows.append(
            f'<tr><td>{html.escape(p["template_chain"])}{p["template_resnum"]} {html.escape(p["template_resname"])}</td>'
            f'<td class="sprite-arrow">matches</td>'
            f'<td class="sprite-exact">{html.escape(p["query_chain"])}{p["query_resnum"]} {html.escape(p["query_resname"])}{check}</td></tr>'
        )
    rmsd = f"{hit.rmsd:.2f} Å" if hit.rmsd is not None else "n/a"
    st.markdown(
        f'<div class="sprite-card"><div class="sprite-head">Best local active-site match</div>'
        f'<div class="sprite-meta"><b>{html.escape(hit.pdb_id)}</b> — {html.escape(hit.description or "structural homolog")} · Chain {html.escape(hit.chain_id or "?")} · RMSD {rmsd}</div>'
        f'<table class="sprite-table"><thead><tr><th>Homolog</th><th></th><th>Your protein</th></tr></thead><tbody>{"".join(rows)}</tbody></table>'
        f'<div class="sprite-note">One site-bearing homolog and one compact three-residue local match, in the SPRITE-style format.</div></div>',
        unsafe_allow_html=True,
    )


def _render_viewer(pid: str) -> None:
    safe = urllib.parse.quote(pid.lower(), safe="")
    components.iframe(f"https://molstar.org/viewer/?pdb={safe}", height=760, scrolling=False)


if run_clicked:
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        st.error("Enter a valid 4-character PDB ID, such as 4HHB.")
        st.stop()
    try:
        with st.spinner(f"Fetching {pdb_id} and finding structural matches..."):
            ticket, query_pdb = fs.submit_search_by_pdb_id(pdb_id, mode="tmalign", databases=["pdb100"])
        query_chain = asite.guess_first_chain_id(query_pdb)
        status = st.empty()
        fs.poll_until_complete(ticket, max_wait_seconds=300, on_status=lambda s, e: status.info(f"Foldseek search: {s.lower()} — {e}s elapsed"))
        status.empty()
        with st.spinner("Calculating RMSD for the structural matches..."):
            hits = fs.fetch_results(ticket, databases=["pdb100"])
            fs.populate_missing_rmsd(hits, query_pdb, query_chain_id=query_chain, max_hits=50)
        if not hits:
            st.warning("Foldseek returned no structural matches. Try another PDB entry.")
            st.stop()

        rmsd_hits = sorted([h for h in hits if h.rmsd is not None], key=lambda h: h.rmsd)
        best = _best_hit(hits)
        if best is not None:
            st.subheader("Closest structural match")
            c1, c2, c3, c4 = st.columns(4)
            with c1:
                st.markdown(f'<div class="metric-card"><div class="metric-label">PDB / chain</div><div class="metric-value">{html.escape(best.target_id)}</div><div class="metric-note">lowest RMSD match</div></div>', unsafe_allow_html=True)
            with c2:
                st.markdown(f'<div class="metric-card"><div class="metric-label">RMSD</div><div class="metric-value">{best.rmsd:.2f} Å</div><div class="metric-note">lower is closer</div></div>', unsafe_allow_html=True)
            with c3:
                tm = f"{best.tm_score:.3f}" if best.tm_score is not None else "n/a"
                st.markdown(f'<div class="metric-card"><div class="metric-label">TM-score</div><div class="metric-value">{tm}</div><div class="metric-note">structural similarity</div></div>', unsafe_allow_html=True)
            with c4:
                ident = f"{best.seq_identity * 100:.1f}%" if best.seq_identity is not None else "n/a"
                st.markdown(f'<div class="metric-card"><div class="metric-label">Sequence identity</div><div class="metric-value">{ident}</div><div class="metric-note">aligned residues</div></div>', unsafe_allow_html=True)
        else:
            st.error("The search returned matches, but no valid RMSD could be calculated. No RMSD winner is shown.")

        st.subheader("Structural matches — ranked by RMSD")
        st.dataframe(pd.DataFrame([{
            "Rank": i + 1,
            "PDB / Chain": h.target_id,
            "Description": h.description or "n/a",
            "RMSD (Å)": f"{h.rmsd:.2f}",
            "TM-Score": f"{h.tm_score:.3f}" if h.tm_score is not None else "n/a",
            "Seq. Identity": f"{h.seq_identity * 100:.1f}%" if h.seq_identity is not None else "n/a",
        } for i, h in enumerate(rmsd_hits[:25])]), use_container_width=True, hide_index=True)

        st.subheader("SPRITE-style active-site match")
        with st.spinner("Finding the closest RMSD homolog with an experimentally mappable site..."):
            site_hit, site_pairs = _find_best_site_hit(rmsd_hits, query_pdb, query_chain)
        if site_hit is not None:
            _render_sprite(site_hit, site_pairs)
        else:
            st.info("No structural homolog with a mappable experimental ligand-binding site was found in the returned PDB matches.")

        st.subheader("Predicted active site")
        with st.spinner("Mapping experimentally observed ligand-binding residues..."):
            predicted = asite.predict_active_site(hits, query_pdb_text=query_pdb, query_chain_id=query_chain, top_n_hits=15)
        if predicted:
            confident = [r for r in predicted if r.support_count >= 2]
            site = confident if confident else predicted[:15]
            text = ", ".join(f"{r.query_resname or '?'}{r.display_resnum}" for r in site)
            st.markdown(f'<div class="site-callout"><div class="site-title">Most strongly supported residues</div><div class="site-residues">{html.escape(text)}</div><div class="small-muted">Highest-supported residue has evidence from {site[0].support_count} structural template(s).</div></div>', unsafe_allow_html=True)
        else:
            st.info("No active-site residues could be transferred from the experimental structural matches.")

        st.subheader("3D structure")
        _render_viewer(pdb_id)
    except fs.FoldseekError as exc:
        st.error(f"Structural search failed: {exc}")
    except Exception as exc:
        st.error(f"Unexpected error: {type(exc).__name__}: {exc}")
