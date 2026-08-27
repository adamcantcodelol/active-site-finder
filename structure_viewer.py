"""Safe HTML renderer for the embedded 3D structure viewer.

The viewer uses 3Dmol.js in the browser. PDB references may arrive from
Foldseek as values such as ``4xz1-assembly1.cif.gz_A``; this module extracts
the actual four-character PDB accession before contacting RCSB.
"""
from __future__ import annotations

import json
import re
from urllib.parse import quote

_PDB_RE = re.compile(r"^[0-9A-Za-z]{4}$")
_CHAIN_RE = re.compile(r"^[A-Za-z0-9_]+$")
_RES_RE = re.compile(r"^[+-]?\d+[A-Za-z]?$" )


def normalize_pdb_id(pdb_id: str) -> str:
    """Normalize either a plain PDB ID or a Foldseek target reference."""
    raw = str(pdb_id).strip()
    if not raw:
        raise ValueError("PDB ID is empty")
    value = raw.split()[0]
    candidate = value[:4]
    if not _PDB_RE.fullmatch(candidate):
        raise ValueError("PDB ID must be exactly four letters/numbers")
    return candidate.lower()


def normalize_residues(residues) -> list[str]:
    """Validate, deduplicate, and preserve residue identifiers."""
    out: list[str] = []
    seen: set[str] = set()
    for value in residues or []:
        text = str(value).strip()
        if _RES_RE.fullmatch(text) and text not in seen:
            seen.add(text)
            out.append(text)
    return out


def viewer_html(pdb_id: str, chain: str | None = None, residues=None, height: int = 560) -> str:
    """Return a 3Dmol.js viewer with an optional local residue highlight."""
    pdb = normalize_pdb_id(pdb_id)
    chain_value = None
    if chain and chain != "?":
        chain_value = str(chain).strip()
        if not _CHAIN_RE.fullmatch(chain_value):
            raise ValueError("Invalid chain identifier")
    selected = normalize_residues(residues)
    safe_height = max(300, min(int(height), 900))
    pdb_url = f"https://files.rcsb.org/download/{quote(pdb.upper())}.pdb"
    return f'''<!doctype html>
<html><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<script src="https://3dmol.csb.pitt.edu/build/3Dmol-min.js"></script>
<style>
html,body,#viewer{{width:100%;height:100%;margin:0;overflow:hidden;background:#fff}}
#status{{font:14px Arial,sans-serif;padding:18px;color:#64748b}}
</style>
</head><body><div id="viewer"><div id="status">Loading structure…</div></div>
<script>
(function() {{
  const chain = {json.dumps(chain_value)};
  const selected = {json.dumps(selected)};
  const url = {json.dumps(pdb_url)};
  const element = document.getElementById('viewer');
  const status = document.getElementById('status');
  if (!window.$3Dmol) {{ status.textContent = '3D viewer library failed to load.'; return; }}
  const viewer = $3Dmol.createViewer(element, {{backgroundColor:'white'}});
  $3Dmol.download(url, viewer, {{}}, function() {{
    if (!viewer.getModel()) {{ status.textContent = 'Structure could not be loaded.'; return; }}
    status.remove();
    viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray'}}}});
    if (chain && selected.length) {{
      const sel = {{chain: chain, resi: selected}};
      viewer.setStyle(sel, {{cartoon: {{color: 'orange'}}, stick: {{radius: 0.22}}, sphere: {{radius: 0.32}}}});
      viewer.zoomTo(sel);
    }} else if (chain) {{
      viewer.zoomTo({{chain: chain}});
    }} else {{
      viewer.zoomTo();
    }}
    viewer.render();
  }});
}})();
</script></body></html>'''
