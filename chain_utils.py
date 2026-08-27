"""Helpers for selecting and validating protein chains."""
from __future__ import annotations

from collections.abc import Iterable

AA3 = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
}


def normalize_chain(chain) -> str | None:
    if chain is None:
        return None
    value = str(chain).strip()
    return value or None


def choose_chain(chains: Iterable[dict], preferred: str | None = None) -> dict | None:
    """Choose a protein chain deterministically, honoring an explicit choice."""
    items = [c for c in chains if isinstance(c, dict)]
    if preferred:
        preferred = normalize_chain(preferred)
        for item in items:
            if normalize_chain(item.get("chain")) == preferred:
                return item
    protein = [c for c in items if c.get("is_protein", True)]
    pool = protein or items
    if not pool:
        return None
    return max(pool, key=lambda c: (int(c.get("length", 0) or 0), str(c.get("chain", ""))))


def validate_chain_residues(residues: Iterable[dict]) -> bool:
    """Reject obviously non-protein or malformed residue records."""
    seen = False
    for residue in residues:
        if not isinstance(residue, dict):
            return False
        name = str(residue.get("resname", "")).upper()
        number = residue.get("resnum")
        if name not in AA3 or number is None:
            return False
        seen = True
    return seen
