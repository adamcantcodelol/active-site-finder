"""Small dependency-free HTML renderer for an embedded 3D structure viewer.

Uses 3Dmol.js in the browser so the selected chain/residues can be highlighted
without requiring another Python package. The viewer is deliberately isolated
from the prediction logic so it cannot alter residue mappings.
"""
from __future__ import annotations

import json
import re


_PDB_RE = re.compile(r"^[0-9A-Za-z]{4}$")
_CHAIN_RE = re.compile(r"^[A-Za-z0-9_]+$")


def normalize_pdb_id(pdb_id: str) -> str:
    value = str(pdb_id).strip().lower()
    if not _PDB_RE.fullmatch(value):
        raise ValueError("PDB ID must be exactly four letters/numbers")
    return value


def normalize_residues(residues) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for value in residues or []:
        text = str(value).strip()
        if re.fullmatch(r"[+-]?\d+[A-Za-z]?", text) and text not in seen:
            seen.add(text)
            out.append(text)
    return out


def viewer_html(pdb_id: str, chain: str | None = None, residues=None, height: int = 560) -> str:
    """Return a self-contained 3Dmol.js viewer with optional highlighted residues."""
    pdb = normalize_pdb_id(pdb_id)
    chain_value = None
    if chain and chain != "?":
        chain_value = str(chain).strip()
        if not _CHAIN_RE.fullmatch(chain_value):
            raise ValueError("Invalid chain identifier")
    selected = normalize_residues(residues)
    safe_height = max(300, min(int(height), 900))
    chain_js = json.dumps(chain_value)
    residues_js = json.dumps(selected)
    pdb_js = json.dumps(pdb)
    return f'''<!doctype html>
<html><head>
<meta charset="utf-8">
<script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
<style>html,body,#viewer{{width:100%;height:100%;margin:0;overflow:hidden;background:#fff}}</style>
</head><body><div id="viewer"></div>
<script>
(function() {{
  const pdb = {pdb_js};
  const chain = {chain_js};
  const selected = {residues_js};
  const element = document.getElementById('viewer');
  const viewer = $3Dmol.createViewer(element, {{backgroundColor:'white'}});
  $3Dmol.download('pdb:' + pdb, viewer, {{}}, function() {{
    viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray'}}}});
    if (chain && selected.length) {{
      viewer.setStyle({{chain: chain, resi: selected}}, {{cartoon: {{color: 'orange'}}, stick: {{radius: 0.22}}, sphere: {{radius: 0.32}}}}});
      viewer.zoomTo({{chain: chain, resi: selected}});
    }} else if (chain) {{
      viewer.zoomTo({{chain: chain}});
    }} else {{
      viewer.zoomTo();
    }}
    viewer.render();
  }});
}})();
</script></body></html>'''
