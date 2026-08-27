"""MBRC Active Site Finder."""
from __future__ import annotations

import html

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

import active_site as asite
import foldseek_client as fs
from chimera_commands import chimera_select
from site_match import SitePair, choose_local_triplet, correspondence_is_valid
from structure_viewer import viewer_html

st.set_page_config(page_title="MBRC Active Site Finder", page_icon="🧬", layout="wide")

st.markdown("""
<style>
.block-container{max-width:1400px;padding-top:2.7rem;padding-bottom:4rem}
.mbrc-header{display:flex;align-items:center;gap:14px;min-height:104px;overflow:visible}
.mbrc-logo{width:226px;height:104px;display:block;overflow:visible;flex:none}
.mbrc-title{font-size:1.8rem;font-weight:650;line-height:1.05;letter-spacing:-.03em}
.mbrc-subtitle{color:#64748b;font-size:1rem;margin:.25rem 0 1.5rem}
.metric-card{border:1px solid #dbe3ee;border-radius:12px;padding:16px 18px;background:#f8fafc;min-height:105px}
.metric-label{color:#64748b;font-size:.82rem;margin-bottom:5px}.metric-value{font-size:1.55rem;font-weight:700}.metric-note{color:#64748b;font-size:.78rem;margin-top:3px}
.sprite-card{border:1px solid #dbe3ee;border-radius:12px;background:#fff;padding:18px 20px;margin:8px 0 12px;overflow-x:auto}
.sprite-head{font-size:1.1rem;font-weight:700}.sprite-meta{color:#64748b;font-size:.84rem;margin:4px 0 14px}
.sprite-table{border-collapse:collapse;width:100%;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
.sprite-table th{text-align:left;color:#64748b;font-family:Arial,Helvetica,sans-serif;font-size:.78rem;border-bottom:1px solid #e2e8f0;padding:7px 9px}
.sprite-table td{padding:9px;border-bottom:1px solid #eef2f7;font-size:.93rem}.sprite-table tr:last-child td{border-bottom:0}
.sprite-exact{font-weight:700}.sprite-arrow{color:#64748b;padding:0 10px}.sprite-note{color:#64748b;font-family:Arial,Helvetica,sans-serif;font-size:.8rem;margin-top:10px}
.site-callout{border:1px solid #cbd5e1;border-radius:12px;padding:18px;background:#f8fafc;margin:8px 0 12px}.site-title{font-size:1.15rem;font-weight:700}.site-residues{font-size:1.05rem;font-weight:650;word-break:break-word}.small-muted{color:#64748b;font-size:.84rem}
.chimera-box{border:1px solid #dbe3ee;border-radius:10px;padding:12px 14px;background:#fafcff;margin:8px 0 14px}
.viewer-card{border:1px solid #dbe3ee;border-radius:12px;padding:10px;background:#fff}
.viewer-title{font-weight:700;margin:2px 0 8px}
</style>
""", unsafe_allow_html=True)

logo="""<svg class='mbrc-logo' viewBox='0 0 470 140' xmlns='http://www.w3.org/2000/svg' role='img' aria-label='MBRC logo'>
<path d='M92 12 L162 30 L162 62 C162 94 140 119 92 132 C44 119 22 94 22 62 L22 30 Z' fill='none' stroke='#111827' stroke-width='5' stroke-linejoin='round'/>
<text x='92' y='86' text-anchor='middle' font-family='Arial,Helvetica,sans-serif' font-size='60' font-weight='800' fill='#111827'>M</text>
<text x='204' y='86' font-family='Arial,Helvetica,sans-serif' font-size='54' font-weight='800' letter-spacing='-3' fill='#111827'>BRC</text>
</svg>"""
st.markdown(f'<div class="mbrc-header">{logo}<div class="mbrc-title">Active Site Finder</div></div><div class="mbrc-subtitle">Structure-first active-site prediction — ranked by the closest RMSD match.</div>',unsafe_allow_html=True)
st.markdown("Enter a **PDB ID**. The app searches experimental PDB structures with Foldseek, ranks them by **RMSD (lower is better)**, and maps experimentally supported active-site residues to your protein.")

pdb_id=st.text_input("PDB ID",placeholder="e.g. 4HHB",max_chars=4).strip().upper()
run=st.button("Find active site",type="primary")


def valid_hits(hits):
    return [h for h in hits if h.rmsd is not None and h.q_aln and h.t_aln]


def best_hit(hits):
    candidates=valid_hits(hits)
    return min(candidates,key=lambda h:(h.rmsd,-(h.tm_score if h.tm_score is not None else -1.0))) if candidates else None


def _residue_key(r):
    return (r.get("resnum"), r.get("insertion_code") or "")


def map_site(hit,query_pdb,query_chain):
    """Map experimentally annotated homolog residues through the canonical validated alignment."""
    try:
        target_pdb=asite.fetch_target_structure(hit)
        if not target_pdb or not hit.q_aln or not hit.t_aln:
            return []
        mapped=asite.map_binding_site_details(hit,target_pdb,query_pdb,query_chain)
        return [SitePair(
            homolog_chain=p.get("tchain",hit.chain_id or "?"),
            homolog_resnum=p["tn"], homolog_insertion=p.get("ticode","") or "",
            homolog_resname=p["tname"], query_chain=p.get("qchain",query_chain or "?"),
            query_resnum=p["qn"], query_insertion=p.get("qicode","") or "",
            query_resname=p["qname"], exact=bool(p.get("exact")),
        ) for p in mapped]
    except Exception:
        return []


def cached_map_site(hit,query_pdb,query_chain):
    """Cache the same canonical mapping used to determine SPRITE availability."""
    cache=st.session_state.setdefault("site_map_cache",{})
    key=(hit.target_id,query_chain,hit.q_aln,hit.t_aln,hit.q_start,hit.t_start)
    if key not in cache:
        pairs=map_site(hit,query_pdb,query_chain)
        # The availability flag and displayed SPRITE match must use the exact
        # same canonical mapping. This fallback also repairs older cached hit
        # metadata produced by previous app versions.
        if len(choose_local_triplet(pairs)) != 3:
            canonical=asite.get_sprite_match(hit,query_pdb,query_chain)
            pairs=[SitePair(
                homolog_chain=p.get("tchain",hit.chain_id or "?"),
                homolog_resnum=p["tn"], homolog_insertion=p.get("ticode","") or "",
                homolog_resname=p["tname"], query_chain=p.get("qchain",query_chain or "?"),
                query_resnum=p["qn"], query_insertion=p.get("qicode","") or "",
                query_resname=p["qname"], exact=bool(p.get("exact")),
            ) for p in canonical]
        cache[key]=pairs
    return cache[key]


def reveal_button(label,command,key,help_text):
    if st.button(label,key=key,use_container_width=True):
        st.session_state[f"show_{key}"]=True
    if st.session_state.get(f"show_{key}"):
        st.markdown('<div class="chimera-box">',unsafe_allow_html=True)
        st.code(command,language="text")
        st.caption(help_text)
        st.markdown('</div>',unsafe_allow_html=True)


def render_sprite(hit,pairs,query_id):
    chosen=choose_local_triplet(pairs)
    if len(chosen)!=3:
        st.info("This homolog does not have a validated local three-residue experimental site that can be mapped to the original protein. Try another homolog.")
        return None
    rows=[]
    for p in chosen:
        mark=" ✓" if p.exact else ""
        rows.append(
            f'<tr><td>{html.escape(p.homolog_chain)}{p.homolog_resnum}{html.escape(p.homolog_insertion)} {html.escape(p.homolog_resname)}</td>'
            f'<td class="sprite-arrow">matches</td>'
            f'<td class="sprite-exact">{html.escape(p.query_chain)}{p.query_resnum}{html.escape(p.query_insertion)} {html.escape(p.query_resname)}{mark}</td></tr>'
        )
    st.markdown(
        f'<div class="sprite-card"><div class="sprite-head">SPRITE-style local match</div>'
        f'<div class="sprite-meta"><b>{html.escape(hit.pdb_id)}</b> — {html.escape(hit.description or "structural homolog")} · Chain {html.escape(hit.chain_id or "?")} · RMSD {hit.rmsd:.2f} Å</div>'
        f'<table class="sprite-table"><thead><tr><th>Homolog</th><th></th><th>{html.escape(query_id)} (original)</th></tr></thead>'
        f'<tbody>{"".join(rows)}</tbody></table><div class="sprite-note">Exactly one local three-residue experimental-site match is shown. ✓ means the amino-acid identity is conserved.</div></div>',
        unsafe_allow_html=True,
    )
    homolog_cmd=chimera_select("#2",chosen[0].homolog_chain, [p.homolog_resnum for p in chosen])
    original_cmd=chimera_select("#1",chosen[0].query_chain, [p.query_resnum for p in chosen])
    a,b=st.columns(2)
    with a:
        reveal_button("ChimeraX: select homolog match",homolog_cmd,"sprite_homolog",f"Open the selected homolog as model #2. Chain {chosen[0].homolog_chain} is selected.")
    with b:
        reveal_button("ChimeraX: select original match",original_cmd,"sprite_original",f"Open the original protein as model #1. Chain {chosen[0].query_chain} is selected.")
    return chosen


def render_viewer(pid,title,chain=None,residues=None,height=560):
    try:
        markup=viewer_html(pid,chain=chain,residues=residues,height=height)
    except ValueError as exc:
        st.error(f"3D viewer configuration error: {exc}")
        return
    st.markdown(f'<div class="viewer-card"><div class="viewer-title">{html.escape(title)}</div>',unsafe_allow_html=True)
    components.html(markup,height=height,scrolling=False)
    st.markdown('</div>',unsafe_allow_html=True)


if run:
    if len(pdb_id)!=4 or not pdb_id.isalnum():
        st.error("Enter a valid 4-character PDB ID, such as 4HHB.");st.stop()
    try:
        with st.spinner(f"Fetching {pdb_id} and finding structural matches..."):
            ticket,query_pdb=fs.submit_search_by_pdb_id(pdb_id,mode="tmalign",databases=["pdb100"])
        query_chain=asite.guess_first_chain_id(query_pdb)
        if not query_chain:
            st.error("The PDB entry was downloaded, but no protein C-alpha chain could be parsed.");st.stop()
        fs.poll_until_complete(ticket,max_wait_seconds=300)
        with st.spinner("Calculating RMSD for the structural matches..."):
            hits=fs.fetch_results(ticket,databases=["pdb100"])
            fs.populate_missing_rmsd(hits,query_pdb,query_chain_id=query_chain,max_hits=50)
        rmsd_hits=sorted(valid_hits(hits),key=lambda h:(h.rmsd,-(h.tm_score if h.tm_score is not None else -1.0)))
        if not rmsd_hits:
            st.error("Foldseek returned matches, but none had a valid RMSD and alignment. No scientifically meaningful ranking can be shown.");st.stop()
        st.session_state.update(query_id=pdb_id,query_pdb=query_pdb,query_chain=query_chain,hits=hits,rmsd_hits=rmsd_hits)
        st.session_state["selected_hit_index"]=0
        st.session_state["site_map_cache"]={}
        st.session_state.pop("homolog_selector",None)
        for key in ("sprite_homolog","sprite_original","predicted_site"):
            st.session_state.pop(f"show_{key}",None)
        with st.spinner("Mapping experimentally observed ligand-binding residues..."):
            st.session_state["predicted_site_cache"]=asite.predict_active_site(hits,query_pdb_text=query_pdb,query_chain_id=query_chain,top_n_hits=15)
    except fs.FoldseekError as exc:st.error(f"Structural search failed: {exc}")
    except Exception as exc:st.error(f"Unexpected error: {type(exc).__name__}: {exc}")


if "rmsd_hits" in st.session_state:
    query_id=st.session_state["query_id"];query_pdb=st.session_state["query_pdb"];query_chain=st.session_state["query_chain"]
    hits=st.session_state["hits"];rmsd_hits=st.session_state["rmsd_hits"]
    best=best_hit(hits)

    st.subheader("Closest structural match")
    c1,c2,c3,c4=st.columns(4)
    with c1:st.markdown(f'<div class="metric-card"><div class="metric-label">PDB / chain</div><div class="metric-value">{html.escape(best.target_id)}</div><div class="metric-note">lowest valid RMSD</div></div>',unsafe_allow_html=True)
    with c2:st.markdown(f'<div class="metric-card"><div class="metric-label">RMSD</div><div class="metric-value">{best.rmsd:.2f} Å</div><div class="metric-note">lower is closer</div></div>',unsafe_allow_html=True)
    with c3:
        tm=f"{best.tm_score:.3f}" if best.tm_score is not None else "n/a"
        st.markdown(f'<div class="metric-card"><div class="metric-label">TM-score</div><div class="metric-value">{tm}</div><div class="metric-note">structural similarity</div></div>',unsafe_allow_html=True)
    with c4:
        ident=f"{best.seq_identity*100:.1f}%" if best.seq_identity is not None else "n/a"
        st.markdown(f'<div class="metric-card"><div class="metric-label">Sequence identity</div><div class="metric-value">{ident}</div><div class="metric-note">aligned residues</div></div>',unsafe_allow_html=True)

    st.subheader("Structural matches — ranked by RMSD")
    st.caption("Ranking is based only on valid structural RMSD. A low-RMSD structure can still lack an experimentally mappable active-site annotation.")
    st.dataframe(pd.DataFrame([{"Rank":i+1,"PDB / Chain":h.target_id,"Description":h.description or "n/a","RMSD (Å)":f"{h.rmsd:.2f}","TM-Score":f"{h.tm_score:.3f}" if h.tm_score is not None else "n/a","Seq. Identity":f"{h.seq_identity*100:.1f}%" if h.seq_identity is not None else "n/a","SPRITE site":"✓ available" if getattr(h,"sprite_available",False) else "not available"} for i,h in enumerate(rmsd_hits[:50])]),use_container_width=True,hide_index=True)

    st.subheader("SPRITE-style active-site match")
    # Determine availability from the exact mapping that will be displayed.
    # This prevents a stale sprite_available flag from offering a homolog
    # that subsequently renders "no validated local site".
    availability={}
    with st.spinner("Checking local SPRITE matches..."):
        for h in rmsd_hits[:50]:
            availability[h.target_id]=len(choose_local_triplet(cached_map_site(h,query_pdb,query_chain)))==3
    available_hits=[h for h in rmsd_hits if availability.get(h.target_id,False)]
    st.caption(f"{len(available_hits)} of {len(rmsd_hits)} structural matches have a validated three-residue experimental site. Only those are selectable below.")
    if not available_hits:
        st.info("None of the returned structural matches have a validated three-residue experimental site that can be mapped to the original protein.")
        chosen=None
        selected_hit=None
    else:
        labels=[f"{rmsd_hits.index(h)+1}. {h.target_id} · chain {h.chain_id or '?'} · RMSD {h.rmsd:.2f} Å" for h in available_hits]
        previous=st.session_state.get("selected_hit_target")
        default_index=next((i for i,h in enumerate(available_hits) if h.target_id==previous),0)
        selected_label=st.selectbox("Choose the similar protein",labels,index=default_index,key="homolog_selector")
        selected_index=labels.index(selected_label)
        selected_hit=available_hits[selected_index]
        st.session_state["selected_hit_target"]=selected_hit.target_id
        selected_pairs=cached_map_site(selected_hit,query_pdb,query_chain)
        chosen=render_sprite(selected_hit,selected_pairs,query_id)

    st.subheader("Predicted active site")
    predicted=st.session_state.get("predicted_site_cache",[])
    if predicted:
        confident=[r for r in predicted if r.support_count>=2]
        site=confident if confident else predicted[:15]
        text=", ".join(f"{r.query_resname or '?'}{r.display_resnum}" for r in site)
        st.markdown(f'<div class="site-callout"><div class="site-title">Most strongly supported residues</div><div class="site-residues">{html.escape(text)}</div><div class="small-muted">Highest-supported residue has evidence from {site[0].support_count} structural template(s). Experimental PDB SITE records are preferred; ligand-contact inference is used when SITE records are absent.</div></div>',unsafe_allow_html=True)
        pred_cmd=chimera_select("#1",query_chain,[r.query_resnum for r in site])
        reveal_button("ChimeraX: select predicted active site",pred_cmd,"predicted_site",f"Open the original protein as model #1. This selects the predicted residues on chain {query_chain}.")
    else:
        st.info("No experimentally annotated active-site residues could be transferred from the returned structural matches.")

    st.subheader("3D structures")
    v1,v2=st.columns(2)
    chosen_homolog_residues=[p.homolog_resnum for p in chosen] if chosen else []
    chosen_query_residues=[p.query_resnum for p in chosen] if chosen else []
    with v1:
        render_viewer(query_id,"Original protein",chain=query_chain,residues=chosen_query_residues)
    with v2:
        homolog_id=selected_hit.pdb_id if selected_hit else None
        if homolog_id:
            render_viewer(homolog_id,"Selected structural homolog",chain=selected_hit.chain_id,residues=chosen_homolog_residues)
        else:
            st.info("Select a structural homolog to view it.")
