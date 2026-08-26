# MBRC Active Site Finder

A Streamlit app that finds the closest structural matches for a PDB structure and predicts the query active site by transferring experimentally observed ligand-binding residues from close structural templates.

## What changed

- Searches Foldseek's **PDB100** database in **TM-align mode**.
- Ranks structural matches primarily by **RMSD (lowest first)**.
- Displays the closest-RMSD structure prominently.
- Correctly reads Foldseek alignment fields (`qaln` / `taln`) so residue mapping can work.
- Keeps TM-score separate from Foldseek's homology probability (`prob`).
- Normalizes sequence identity safely.
- Uses actual PDB CA-residue order when transferring binding sites, so PDB numbering gaps/insertion codes are less likely to break the mapping.
- Uses PDBe binding-site annotations from experimental structures.
- Highlights predicted active-site residues in the 3Dmol viewer.
- Fixes the clipped bottom 3D viewer by giving the iframe enough height.
- Adds the supplied shield logo and MBRC-style header.
- Improves error handling and Foldseek polling timeout behavior.

## Run

```bash
pip install -r requirements.txt
streamlit run app.py
```

Then enter a four-character PDB ID such as `4HHB`.

## Important interpretation note

The app is a structure-based prediction, not a proof of catalytic activity. The strongest predictions are residues where multiple close experimental structures with known ligand-binding sites independently map to the same query position. Always verify the result against the structure and biological context.
