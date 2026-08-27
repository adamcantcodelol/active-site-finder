# MBRC Active Site Finder

A Streamlit app that finds the closest structural matches for a PDB structure and predicts the query active site by transferring experimentally observed ligand-binding residues from experimental structural templates.

## What changed

- Searches Foldseek's **PDB100** database in **TM-align mode**.
- Ranks structural matches primarily by **validated RMSD (lowest first)**.
- Displays the closest-RMSD structure separately from the best site-bearing structural match.
- Correctly reads Foldseek alignment fields (`qaln` / `taln`) so residue mapping can work.
- Validates alignment/residue correspondence against the actual deposited C-alpha sequence before accepting a mapping.
- Keeps TM-score separate from Foldseek's homology probability (`prob`).
- Normalizes sequence identity safely.
- Uses actual PDB CA-residue order when transferring binding sites, including PDB numbering gaps and insertion codes.
- Uses experimental PDB/PDBe binding-site annotations with ligand-contact inference as a fallback.
- Provides an interactive **SPRITE-style local match**: select a specific homolog and see one compact three-residue correspondence for that homolog only.
- Generates ChimeraX commands for the original predicted site, the selected homolog's local site, and the original protein's corresponding local site.
- Handles missing/invalid annotations without inventing a residue correspondence.
- Highlights predicted active-site residues in the 3D viewer.
- Fixes the clipped bottom 3D viewer by giving the iframe enough height.
- Uses the supplied shield logo with the **M inside the shield** and BRC outside it.
- Improves error handling and Foldseek polling timeout behavior.
- Includes regression tests for alignment/site selection and ChimeraX command generation.

## Run

```bash
pip install -r requirements.txt
streamlit run app.py
```

Then enter a four-character PDB ID such as `4HHB`.

## ChimeraX workflow

The generated commands assume the **original protein is model #1** and the **selected homolog is model #2**. Open those structures in ChimeraX in that order, then paste the corresponding command into the ChimeraX command line.

The SPRITE-style buttons select only the three residues displayed in the local match. Insertion codes are preserved when present.

## Important interpretation note

The app is a structure-based prediction, not a proof of catalytic activity. A low-RMSD structural match and a site-bearing homolog are not necessarily the same structure. The strongest predictions are residues where multiple experimental structures with known ligand-binding sites independently map to the same query position. Always verify the result against the structure and biological context.
