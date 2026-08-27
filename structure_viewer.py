"""3D structure viewer URL/HTML helpers.

The hosted 3Dmol viewer is used here rather than running custom JavaScript in a
Streamlit component iframe.  This avoids browser/CDN/CORS failures inside the
Streamlit sandbox while still allowing exact chain/residue styling through
3Dmol's documented URL selectors.
"""
from __future__ import annotations

import html
import re
from urllib.parse import quote

_PDB_RE = re.compile(r"^[0-9A-Za-z]{4}$")
_CHAIN_RE = re.compile(r"^[A-Za-z0-9_]+$")
_RES_RE = re.compile(r"^[+-]?\d+[A-Za-z]?$")


def normalize_pdb_id(pdb_id: str) -> str:
    raw = str(pdb_id).strip()
    if not raw:
        raise ValueError("PDB ID is empty")
    candidate = raw.split()[0][:4]
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


def viewer_url(pdb_id: str, chain: str | None = None, residues=None) -> str:
    """Build a hosted 3Dmol scene with the selected active-site residues highlighted."""
    pdb = normalize_pdb_id(pdb_id)
    selected = normalize_residues(residues)
    chain_value = None
    if chain and chain != "?":
        chain_value = str(chain).strip()
        if not _CHAIN_RE.fullmatch(chain_value):
            raise ValueError("Invalid chain identifier")

    base = "https://3dmol.org/viewer.html"
    params = [f"pdb={quote(pdb.upper())}", "style=cartoon:color~lightgray"]
    if chain_value and selected:
        selector = f"resi:{','.join(selected)};chain:{chain_value}"
        # Magenta is deliberately distinct from the neutral protein cartoon.
        params.extend([
            f"select={quote(selector, safe=':;,')}",
            "style=stick:radius~0.28,colorscheme~magentaCarbon",
            f"select={quote(selector, safe=':;,')}",
            "style=sphere:radius~0.38,colorscheme~magentaCarbon",
            f"select={quote(selector, safe=':;,')}",
            "labelres=fontSize:12;backgroundOpacity:0.65",
        ])
    elif chain_value:
        params.extend([
            f"select={quote('chain:' + chain_value, safe=':')}",
            "style=cartoon:color~lightgray",
        ])
    return base + "?" + "&".join(params)


def viewer_html(pdb_id: str, chain: str | None = None, residues=None, height: int = 560) -> str:
    """Return a Streamlit-safe iframe embedding the hosted 3Dmol viewer."""
    url = viewer_url(pdb_id, chain=chain, residues=residues)
    safe_height = max(300, min(int(height), 900))
    return (
        f'<iframe src="{html.escape(url, quote=True)}" '
        f'width="100%" height="{safe_height}" '
        f'style="border:0;border-radius:10px;background:#fff" '
        f'loading="lazy" allowfullscreen></iframe>'
    )
