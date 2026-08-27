"""Helpers for generating safe ChimeraX residue-selection commands.

The UI should pass validated residue numbers from the PDB/alignment layer into
these helpers. Keeping command formatting here prevents the UI from having
slightly different ChimeraX syntax in different places.
"""
from __future__ import annotations

from typing import Iterable, Optional


def _clean_chain(chain: Optional[str]) -> str:
    chain = (chain or "").strip()
    # ChimeraX atom-spec chain identifiers are normally one character. Do not
    # allow whitespace or punctuation to become part of a generated command.
    if len(chain) == 1 and chain.isalnum():
        return chain
    return ""


def _clean_model(model: str) -> str:
    model = str(model).strip()
    if not model:
        return "#1"
    if model.startswith("#"):
        body = model[1:]
    else:
        body = model
    # Permit model paths such as #1.2, but only numeric components.
    parts = body.split(".")
    if not parts or any(not p.isdigit() for p in parts):
        return "#1"
    return "#" + ".".join(parts)


def residue_spec(resnum: int, insertion_code: str = "") -> str:
    """Return ChimeraX's residue-number token, preserving insertion codes."""
    number = int(resnum)
    icode = (insertion_code or "").strip()
    if len(icode) > 1 or (icode and not icode.isalnum()):
        icode = ""
    return f"{number}{icode}"


def select_command(
    model: str,
    chain: Optional[str],
    residues: Iterable[tuple[int, str] | int],
) -> str:
    """Build a deterministic ChimeraX select command.

    Invalid residue values are ignored. An empty residue collection returns
    ``select clear`` so the caller never emits a syntactically incomplete
    command.
    """
    clean: set[str] = set()
    for item in residues:
        try:
            if isinstance(item, tuple):
                number, insertion = item
                clean.add(residue_spec(int(number), str(insertion or "")))
            else:
                clean.add(residue_spec(int(item)))
        except (TypeError, ValueError):
            continue

    if not clean:
        return "select clear"

    ordered = sorted(clean, key=lambda value: (int(''.join(c for c in value if c.isdigit())), value))
    chain_id = _clean_chain(chain)
    chain_part = f"/{chain_id}" if chain_id else ""
    return f"select {_clean_model(model)}{chain_part}:{','.join(ordered)}"
