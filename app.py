"""
MBRC Active Site Finder

Primary goal:
    Find the closest structural matches to a query PDB structure by RMSD,
    then transfer experimentally observed ligand-binding residues from the
    closest matching PDB structures onto the query.
"""

from __future__ import annotations

import html
import os

import pandas as pd
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


# ---------------------------------------------------------------------------
# Styling
# ---------------------------------------------------------------------------

st.markdown(
    """
    <style>
    .block-container {
        max-width: 1400px;
        padding-top: 1.5rem;
        padding-bottom: 3rem;
    }
    .mbrc-header {
        display: flex;
        align-items: center;
        gap: 14px;
        margin-bottom: 0.25rem;
    }
    .mbrc-logo {
        width: 58px;
        height: 58px;
        object-fit: contain;
        border-radius: 8px;
    }
    .mbrc-name {
        font-size: 2rem;
        font-weight: 750;
        line-height: 1;
        letter-spacing: -0.03em;
    }
    .mbrc-subtitle {
        color: #64748b;
        font-size: 1rem;
        margin: 0.35rem 0 1.5rem 72px;
    }
    .metric-card {
        border: 1px solid #dbe3ee;
        border-radius: 12px;
        padding: 16px 18px;
        background: #f8fafc;
        min-height: 105px;
    }
    .metric-label {
        color: #64748b;
        font-size: 0.82rem;
        margin-bottom: 5px;
    }
    .metric-value {
        font-size: 1.55rem;
        font-weight: 700;
    }
    .metric-note {
        color: #64748b;
        font-size: 0.78rem;
        margin-top: 3px;
    }
    .site-callout {
        border: 1px solid #cbd5e1;
        border-radius: 12px;
        padding: 18px;
        background: #f8fafc;
        margin: 8px 0 16px 0;
    }
    .site-title {
        font-size: 1.15rem;
        font-weight: 700;
        margin-bottom: 4px;
    }
    .site-residues {
        font-size: 1.05rem;
        font-weight: 650;
        word-break: break-word;
    }
    .small-muted {
        color: #64748b;
        font-size: 0.84rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

logo_path = os.path.join(os.path.dirname(__file__), "assets", "mbrc_logo.png")
if os.path.exists(logo_path):
    with open(logo_path, "rb") as f:
        import base64
        logo_b64 = base64.b64encode(f.read()).decode("ascii")
    logo_html = (
        f'<img class="mbrc-logo" src="data:image/png;base64,{logo_b64}" '
        'alt="MBRC shield logo">'
    )
else:
    logo_html = '<div class="mbrc-logo"></div>'

st.markdown(
    f"""
    <div class="mbrc-header">
        {logo_html}
        <div class="mbrc-name">M<span style="font-weight:500">BRC</span> Active Site Finder</div>
    </div>
    <div class="mbrc-subtitle">
        Structure-first active-site prediction — ranked by the closest RMSD match.
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown(
    "Enter a **PDB ID**. The app searches experimental PDB structures with "
    "Foldseek, ranks them by **RMSD (lower is better)**, and transfers known "
    "ligand-binding residues from the closest structural matches to your protein."
)

pdb_id = st.text_input(
    "PDB ID",
    placeholder="e.g. 4HHB",
    max_chars=4,
    label_visibility="visible",
).strip().upper()

run_clicked = st.button(
    "Find active site",
    type="primary",
    use_container_width=False,
)


# ---------------------------------------------------------------------------
# 3D viewer
# ---------------------------------------------------------------------------

def render_viewer(
    query_pdb: str,
    active_site_resnums: list[int],
    height: int = 620,
) -> None:
    """Render a self-contained 3Dmol viewer with enough iframe height.

    The old viewer was clipped because the iframe was only a few pixels
    taller than the canvas plus caption. This version keeps the viewer and
    legend inside a fixed-height iframe.
    """
    def _escape_js_template(text: str) -> str:
        return (
            text.replace("\\", "\\\\")
            .replace("`", "\\`")
            .replace("</script>", "<\\/script>")
        )

    # 3Dmol accepts a comma-separated residue selector.
    resi_selector = ",".join(str(n) for n in active_site_resnums)
    pdb_js = _escape_js_template(query_pdb)

    html_doc = f"""
    <!doctype html>
    <html>
    <head>
      <meta charset="utf-8">
      <style>
        html, body {{
          margin: 0;
          padding: 0;
          width: 100%;
          height: 100%;
          overflow: hidden;
          background: #ffffff;
          font-family: Arial, sans-serif;
        }}
        #viewer3d {{
          width: 100%;
          height: calc(100% - 38px);
          min-height: 540px;
          position: relative;
        }}
        #legend {{
          height: 38px;
          display: flex;
          align-items: center;
          padding-left: 10px;
          color: #475569;
          font-size: 13px;
          box-sizing: border-box;
        }}
      </style>
      <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    </head>
    <body>
      <div id="viewer3d"></div>
      <div id="legend">
        Gray = protein &nbsp;&nbsp;|&nbsp;&nbsp;
        Orange = predicted active-site residues
      </div>
      <script>
        const el = document.getElementById("viewer3d");
        const viewer = $3Dmol.createViewer(el, {{
          backgroundColor: "white"
        }});

        const pdb = `{pdb_js}`;
        const model = viewer.addModel(pdb, "pdb");

        model.setStyle({{}}, {{
          cartoon: {{color: "lightgray"}}
        }});

        const active = "{resi_selector}";
        if (active) {{
          model.setStyle(
            {{resi: active}},
            {{
              cartoon: {{color: "orange"}},
              stick: {{radius: 0.18, colorscheme: "orangeCarbon"}}
            }}
          );
          model.addStyle(
            {{resi: active}},
            {{
              sphere: {{radius: 0.35, colorscheme: "orangeCarbon"}}
            }}
          );
        }}

        viewer.zoomTo();
        viewer.render();

        window.addEventListener("resize", () => viewer.resize());
      </script>
    </body>
    </html>
    """

    components.html(
        html_doc,
        height=height,
        scrolling=False,
    )


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

if run_clicked:
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        st.error("Enter a valid 4-character PDB ID, such as 4HHB.")
        st.stop()

    try:
        # TM-align mode is intentionally used: it is the Foldseek mode that
        # performs a global structural alignment and gives the most useful
        # RMSD/TM-score comparison for this app's purpose.
        with st.spinner(
            f"Fetching {pdb_id} and finding the closest structural matches..."
        ):
            ticket_id, query_pdb_text = fs.submit_search_by_pdb_id(
                pdb_id,
                mode="tmalign",
                databases=["pdb100"],
            )

        status_box = st.empty()

        def _on_status(status: str, elapsed: int) -> None:
            status_box.info(
                f"Foldseek search: {status.lower()} — {elapsed}s elapsed"
            )

        fs.poll_until_complete(
            ticket_id,
            max_wait_seconds=300,
            on_status=_on_status,
        )
        status_box.empty()

        with st.spinner("Reading structural matches and calculating the RMSD ranking..."):
            hits = fs.fetch_results(ticket_id, databases=["pdb100"])

        if not hits:
            st.warning(
                "Foldseek returned no structural matches. Try a different PDB "
                "entry or a structure with a complete protein chain."
            )
            st.stop()

        # ------------------------------------------------------------------
        # Closest RMSD result
        # ------------------------------------------------------------------
        rmsd_hits = [h for h in hits if h.rmsd is not None]
        if rmsd_hits:
            best = rmsd_hits[0]

            st.subheader("Closest structural match")

            c1, c2, c3, c4 = st.columns(4)
            with c1:
                st.markdown(
                    f"""
                    <div class="metric-card">
                      <div class="metric-label">PDB / chain</div>
                      <div class="metric-value">{html.escape(best.target_id)}</div>
                      <div class="metric-note">lowest RMSD match</div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
            with c2:
                st.markdown(
                    f"""
                    <div class="metric-card">
                      <div class="metric-label">RMSD</div>
                      <div class="metric-value">{best.rmsd:.2f} Å</div>
                      <div class="metric-note">lower is closer</div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
            with c3:
                tm = f"{best.tm_score:.3f}" if best.tm_score is not None else "n/a"
                st.markdown(
                    f"""
                    <div class="metric-card">
                      <div class="metric-label">TM-score</div>
                      <div class="metric-value">{tm}</div>
                      <div class="metric-note">structural similarity</div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
            with c4:
                ident = (
                    f"{best.seq_identity * 100:.1f}%"
                    if best.seq_identity is not None
                    else "n/a"
                )
                st.markdown(
                    f"""
                    <div class="metric-card">
                      <div class="metric-label">Sequence identity</div>
                      <div class="metric-value">{ident}</div>
                      <div class="metric-note">aligned residues</div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )

            if best.description:
                st.caption(best.description)

        else:
            st.warning(
                "Foldseek returned structural matches, but the response did not "
                "include RMSD values. The table is still shown, but an RMSD-based "
                "active-site ranking cannot be guaranteed for this search."
            )

        # ------------------------------------------------------------------
        # Full ranking
        # ------------------------------------------------------------------
        st.subheader("Structural matches — ranked by RMSD")

        if rmsd_hits:
            table_hits = rmsd_hits[:25]
        else:
            table_hits = hits[:25]

        table_df = pd.DataFrame(
            [
                {
                    "Rank": i + 1,
                    "PDB / Chain": h.target_id,
                    "Description": h.description or "n/a",
                    "RMSD (Å)": (
                        f"{h.rmsd:.2f}" if h.rmsd is not None else "n/a"
                    ),
                    "TM-Score": (
                        f"{h.tm_score:.3f}"
                        if h.tm_score is not None
                        else "n/a"
                    ),
                    "Seq. Identity": (
                        f"{h.seq_identity * 100:.1f}%"
                        if h.seq_identity is not None
                        else "n/a"
                    ),
                }
                for i, h in enumerate(table_hits)
            ]
        )
        st.dataframe(
            table_df,
            use_container_width=True,
            hide_index=True,
        )

        # ------------------------------------------------------------------
        # Active-site prediction
        # ------------------------------------------------------------------
        st.subheader("Predicted active site")

        query_chain = asite.guess_first_chain_id(query_pdb_text)

        with st.spinner(
            "Transferring known ligand-binding sites from the closest-RMSD "
            "experimental structures..."
        ):
            predicted = asite.predict_active_site(
                hits,
                query_pdb_text=query_pdb_text,
                query_chain_id=query_chain,
                top_n_hits=15,
            )

        if not predicted:
            st.info(
                "No active-site residues could be transferred. The closest "
                "structures may not have recorded ligand-binding sites, or "
                "Foldseek did not return the alignment strings needed for "
                "residue-by-residue mapping."
            )
        else:
            # If multiple close templates agree, those are the strongest calls.
            confident = [r for r in predicted if r.support_count >= 2]
            shown = confident if confident else predicted[:15]

            residue_text = ", ".join(
                f"{r.query_resname or '?'}{r.display_resnum}"
                for r in shown
            )

            strongest = shown[0]
            st.markdown(
                f"""
                <div class="site-callout">
                  <div class="site-title">Most strongly supported residues</div>
                  <div class="site-residues">{html.escape(residue_text)}</div>
                  <div class="small-muted">
                    The highest-supported residue has evidence from
                    {strongest.support_count} structural template(s).
                    Templates are considered in closest-RMSD order.
                  </div>
                </div>
                """,
                unsafe_allow_html=True,
            )

            site_df = pd.DataFrame(
                [
                    {
                        "Residue": f"{r.query_resname or '?'}{r.display_resnum}",
                        "Query position": r.display_resnum,
                        "Homologs agreeing": r.support_count,
                        "Supporting structures": ", ".join(
                            r.supporting_hits[:6]
                        )
                        + (" ..." if len(r.supporting_hits) > 6 else ""),
                    }
                    for r in shown
                ]
            )
            st.dataframe(
                site_df,
                use_container_width=True,
                hide_index=True,
            )

            st.caption(
                "A residue is more convincing when multiple independent "
                "experimental structures with known ligands map their binding "
                "sites to the same query position."
            )

        # ------------------------------------------------------------------
        # 3D structure
        # ------------------------------------------------------------------
        st.subheader("3D structure")

        active_resnums = [r.query_resnum for r in predicted]
        render_viewer(
            query_pdb_text,
            active_resnums,
            height=650,
        )

    except fs.FoldseekError as exc:
        st.error(f"Structural search failed: {exc}")
    except Exception as exc:  # noqa: BLE001
        st.error(f"Something went wrong: {exc}")
else:
    st.info(
        "Enter a PDB ID above and click **Find active site**. "
        "The closest-RMSD structural match will be shown first."
    )
