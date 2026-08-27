"""UI helpers for structure viewing and ChimeraX selections."""
from __future__ import annotations

import html
import urllib.parse


def molstar_url(pdb_id: str, chain: str | None = None, residues=None) -> str:
    """Build a Mol* viewer URL with optional residue focus."""
    pdb = str(pdb_id).strip().lower()
    params = {"pdb": pdb}
    if chain and chain != "?":
        params["chain"] = str(chain).strip()
    if residues:
        params["residues"] = ",".join(str(r) for r in residues)
    return "https://molstar.org/viewer/?" + urllib.parse.urlencode(params)


def command_block(command: str) -> str:
    """Render a safe, copy-friendly command block."""
    return (
        '<div style="position:relative">'
        '<pre style="white-space:pre-wrap;margin:0">'
        + html.escape(str(command))
        + '</pre></div>'
    )
