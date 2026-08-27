"""Safe ChimeraX command generation for Active Site Finder."""
from __future__ import annotations

import re
from typing import Iterable

_RESNUM_RE = re.compile(r"^[+-]?\d+[A-Za-z]?$|^[+-]?\d+$")
_MODEL_RE = re.compile(r"^#[0-9]+(?:\.[0-9]+)*$")
_CHAIN_RE = re.compile(r"^[A-Za-z0-9_]+$")


def normalize_residue_number(value) -> str | None:
    """Return a ChimeraX-safe residue number, preserving insertion codes."""
    if value is None:
        return None
    text = str(value).strip()
    if not text or not _RESNUM_RE.fullmatch(text):
        return None
    return text


def residue_sort_key(value: str):
    m = re.fullmatch(r"([+-]?\d+)([A-Za-z]?)", value)
    if not m:
        return (10**18, value)
    return (int(m.group(1)), m.group(2).upper())


def normalize_residues(residues: Iterable) -> list[str]:
    """Deduplicate and sort residue numbers deterministically."""
    values = {r for x in residues if (r := normalize_residue_number(x)) is not None}
    return sorted(values, key=residue_sort_key)


def chimera_select(model: str, chain: str | None, residues: Iterable) -> str:
    """Build a validated ChimeraX select command."""
    model = str(model).strip()
    if not _MODEL_RE.fullmatch(model):
        raise ValueError(f"Invalid ChimeraX model: {model!r}")
    chain_text = ""
    if chain is not None and str(chain).strip():
        chain_text = str(chain).strip()
        if not _CHAIN_RE.fullmatch(chain_text):
            raise ValueError(f"Invalid chain identifier: {chain!r}")
        chain_text = f"/{chain_text}"
    nums = normalize_residues(residues)
    if not nums:
        return "select clear"
    return f"select {model}{chain_text}:{','.join(nums)}"


def chimera_select_pairs(model: str, chain: str | None, pairs: Iterable[dict], number_key: str) -> str:
    """Build a selection directly from validated correspondence records."""
    return chimera_select(model, chain, (p.get(number_key) for p in pairs))
