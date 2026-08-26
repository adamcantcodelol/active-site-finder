"""
app.py
======
Active Site Finder — type a PDB ID, get back:
  1. A table of structurally similar proteins, with RMSD and TM-score.
  2. A predicted active site for YOUR protein, inferred from what's known
     about its structural relatives (works even on unstudied/"unknown"
     proteins, since the evidence comes from homologs, not the query).
  3. A 3D view of your protein with the predicted active site highlighted.

HOW TO RUN THIS (no coding experience required)
------------------------------------------------
See README.md for full step-by-step setup instructions. Once set up, the
whole app is one command:

    streamlit run app.py

Then a browser tab opens automatically. Type a 4-character PDB ID (like
4HHB) and click the button. That's the entire user interaction.
"""

from __future__ import annotations

import streamlit as st
import pandas as pd
import streamlit.components.v1 as components

import foldseek_client as fs
import active_site as asite

st.set_page_config(page_title="Active Site Finder", layout="wide")

st.title("🔬 Active Site Finder")
st.caption(
    "Type a PDB ID. This tool finds structurally similar proteins, reports "
    "their RMSD/TM-score, and predicts your protein's active site based on "
    "what's known about its relatives — even if your protein itself has "
    "never been studied."
)

pdb_id = st.text_input(
    "PDB ID", placeholder="e.g. 4HHB", max_chars=4,
).strip()

run_clicked = st.button("Analyze structure", type="primary")


# --------------------------------------------------------------------------
# Helper: embed a 3Dmol.js viewer with active-site residues highlighted
# --------------------------------------------------------------------------

def render_viewer(query_pdb: str, active_site_resnums: list[int], height: int = 480) -> None:
    def _escape(text: str) -> str:
        return text.replace("\\", "\\\\").replace("`", "\\`")

    resi_selector = ",".join(str(n) for n in active_site_resnums) or "-9999"  # empty-safe

    html = f"""
    <div id="viewer3d" style="height:{height}px; width:100%;"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
        let el = document.getElementById("viewer3d");
        let viewer = $3Dmol.createViewer(el, {{backgroundColor: "white"}});
        let model = viewer.addModel(`{_escape(query_pdb)}`, "pdb");
        model.setStyle({{}}, {{cartoon: {{color: "lightgray"}}}});
        model.setStyle({{resi: "{resi_selector}"}}, {{stick: {{colorscheme: "yellowCarbon"}}}});
        model.setStyle({{resi: "{resi_selector}"}}, {{cartoon: {{color: "orange"}}}});
        viewer.zoomTo();
        viewer.render();
    </script>
    <p style="font-size:0.85em; color:#555;">
        Gray = overall protein &nbsp;|&nbsp; Orange/yellow sticks = predicted active site residues
    </p>
    """
    components.html(html, height=height + 40, scrolling=False)


# --------------------------------------------------------------------------
# Main flow
# --------------------------------------------------------------------------

if run_clicked:
    if len(pdb_id) != 4:
        st.error("Please enter a valid 4-character PDB ID, e.g. '4HHB'.")
        st.stop()

    try:
        # Step 1: fetch the query structure and submit it to Foldseek.
        # We search only the "pdb100" database (real experimental structures),
        # since active-site evidence has to come from entries with actual
        # bound ligands — AlphaFold models don't have that information.
        with st.spinner(f"Fetching {pdb_id.upper()} and submitting structural search..."):
            ticket_id, query_pdb_text = fs.submit_search_by_pdb_id(
                pdb_id, mode="3diaa", databases=["pdb100"]
            )

        # Step 2: wait for the search to finish.
        status_box = st.empty()

        def _on_status(status: str, elapsed: int) -> None:
            status_box.info(f"Searching for similar structures... ({status}, {elapsed}s)")

        with st.spinner("This usually takes 10-60 seconds..."):
            fs.poll_until_complete(ticket_id, on_status=_on_status)
        status_box.empty()

        # Step 3: get the list of similar proteins.
        with st.spinner("Fetching similar proteins..."):
            hits = fs.fetch_results(ticket_id, databases=["pdb100"])

        if not hits:
            st.warning(
                "No structurally similar proteins were found. This can happen "
                "for very small, very novel, or poorly-formed structures."
            )
            st.stop()

        st.subheader(f"Similar proteins found: {len(hits)}")
        table_df = pd.DataFrame([
            {
                "PDB ID": h.target_id,
                "Description": h.description or "n/a",
                "TM-Score": f"{h.tm_score:.2f}" if h.tm_score is not None else "n/a",
                "RMSD (\u00c5)": f"{h.rmsd:.2f}" if h.rmsd is not None else "n/a",
                "Seq. Identity": f"{h.seq_identity * 100:.1f}%" if h.seq_identity is not None else "n/a",
            }
            for h in hits
        ])
        st.dataframe(table_df, use_container_width=True, hide_index=True)

        # Step 4: predict the active site from homolog evidence.
        with st.spinner("Predicting active site from similar proteins' known binding sites..."):
            predicted = asite.predict_active_site(hits, top_n_hits=15)
            asite.annotate_residue_names(
                predicted, query_pdb_text,
                chain_id=asite.guess_first_chain_id(query_pdb_text),
            )

        st.subheader("Predicted active site")

        if not predicted:
            st.info(
                "No active-site evidence could be transferred from the similar "
                "proteins found (none of them had a recorded ligand-binding site "
                "in the PDB — they may be unbound/apo structures). Try a query "
                "with more homologs, or interpret this protein's function with "
                "other tools."
            )
        else:
            # Keep residues with at least 2 independent homologs agreeing,
            # unless that would leave nothing to show.
            confident = [r for r in predicted if r.support_count >= 2]
            shown = confident if confident else predicted[:10]

            site_df = pd.DataFrame([
                {
                    "Residue #": r.query_resnum,
                    "Residue": r.query_resname or "n/a",
                    "Homologs agreeing": r.support_count,
                    "Evidence from": ", ".join(r.supporting_hits[:5]) + (
                        "..." if len(r.supporting_hits) > 5 else ""
                    ),
                }
                for r in shown
            ])
            st.dataframe(site_df, use_container_width=True, hide_index=True)
            st.caption(
                "Higher 'homologs agreeing' = stronger evidence. Residue numbers "
                "match the PDB file's own numbering."
            )

        # Step 5: 3D view with the active site highlighted.
        st.subheader("3D structure")
        active_resnums = [r.query_resnum for r in (predicted if predicted else [])]
        render_viewer(query_pdb_text, active_resnums)

    except fs.FoldseekError as exc:
        st.error(f"Structural search failed: {exc}")
    except Exception as exc:  # noqa: BLE001
        st.error(f"Something went wrong: {exc}")

else:
    st.info("Enter a PDB ID above and click **Analyze structure** to begin.")
