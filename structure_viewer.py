"""Safe HTML renderer for the embedded 3D structure viewer."""
from __future__ import annotations

import json
import re
from urllib.parse import quote

_PDB_RE = re.compile(r"^[0-9A-Za-z]{4}$")
_CHAIN_RE = re.compile(r"^[A-Za-z0-9_]+$")
_RES_RE = re.compile(r"^[+-]?\d+[A-Za-z]?$")


def normalize_pdb_id(pdb_id: str) -> str:
    raw = str(pdb_id).strip()
    if not raw:
        raise ValueError("PDB ID is empty")
    value = raw.split()[0]
    candidate = value[:4]
    if not _PDB_RE.fullmatch(candidate):
        raise ValueError("PDB ID must be exactly four letters/numbers")
    return candidate.lower()


def normalize_residues(residues) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for value in residues or []:
        text = str(value).strip()
        if _RES_RE.fullmatch(text) and text not in seen:
            seen.add(text)
            out.append(text)
    return out


def viewer_html(pdb_id: str, chain: str | None = None, residues=None, height: int = 560) -> str:
    """Return a self-contained 3Dmol.js viewer.

    We fetch the PDB text ourselves and call ``addModel`` instead of using
    3Dmol's ``download`` helper.  This gives the iframe an explicit HTTP,
    parsing, and rendering error path and avoids the blank iframe seen with
    some browser/CDN combinations.
    """
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
#status{{font:14px Arial,sans-serif;padding:18px;color:#475569}}
</style>
</head><body><div id="viewer"><div id="status">Loading structure…</div></div>
<script>
(async function() {{
  const chain = {json.dumps(chain_value)};
  const selected = {json.dumps(selected)};
  const url = {json.dumps(pdb_url)};
  const element = document.getElementById('viewer');
  const status = document.getElementById('status');
  function fail(message) {{ status.textContent = message; }}

  try {{
    if (!window.$3Dmol) {{ fail('3D viewer library failed to load.'); return; }}
    const response = await fetch(url, {{cache: 'no-store'}});
    if (!response.ok) {{ fail('Could not download the structure (HTTP ' + response.status + ').'); return; }}
    const pdbText = await response.text();
    if (!pdbText || pdbText.indexOf('ATOM') < 0) {{ fail('The downloaded PDB contains no parseable ATOM records.'); return; }}

    const viewer = $3Dmol.createViewer(element, {{backgroundColor: 'white'}});
    const model = viewer.addModel(pdbText, 'pdb');
    if (!model) {{ fail('3Dmol could not parse the PDB structure.'); return; }}

    viewer.setStyle({{}}, {{cartoon: {{color: 'lightgray'}}}});
    let focus = null;
    if (chain && selected.length) {{
      focus = {{chain: chain, resi: selected}};
      viewer.setStyle(focus, {{cartoon: {{color: 'orange'}}, stick: {{radius: 0.22}}, sphere: {{radius: 0.34}}}});
      viewer.addStyle(focus, {{label: {{fontSize: 11, backgroundOpacity: 0.65, inFront: true}}}});
    }} else if (chain) {{
      focus = {{chain: chain}};
    }}

    if (focus) viewer.zoomTo(focus); else viewer.zoomTo();
    viewer.resize();
    viewer.render();
    status.remove();
  }} catch (error) {{
    fail('3D viewer error: ' + (error && error.message ? error.message : 'unknown error'));
  }}
}})();
</script></body></html>'''
