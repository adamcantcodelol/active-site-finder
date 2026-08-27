"""Helpers for generating residue-focused Mol* viewer URLs."""
from __future__ import annotations

import urllib.parse


def _clean(value) -> str:
    return str(value).strip()


def residue_query(residues) -> str:
    """Normalize residue identifiers for a Mol* URL query."""
    out = []
    seen = set()
    for residue in residues or []:
        value = _clean(residue)
        if not value or value in seen:
            continue
        seen.add(value)
        out.append(value)
    return ",".join(out)


def molstar_focus_url(pdb_id: str, chain: str | None = None, residues=None) -> str:
    """Build a stable Mol* URL, including optional residue focus metadata."""
    pdb = _clean(pdb_id).lower()
    if len(pdb) != 4 or not pdb.isalnum():
        raise ValueError("PDB ID must be a 4-character alphanumeric identifier")
    params = {"pdb": pdb}
    if chain:
        params["chain"] = _clean(chain)
    selected = residue_query(residues)
    if selected:
        params["residues"] = selected
    return "https://molstar.org/viewer/?" + urllib.parse.urlencode(params)


def viewer_selection_payload(chain: str, residues) -> dict:
    """Return serializable viewer-selection metadata for UI integrations."""
    selected = residue_query(residues)
    return {"chain": _clean(chain), "residues": selected.split(",") if selected else []}
