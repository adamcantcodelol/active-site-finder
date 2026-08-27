"""Canonical, testable logic for SPRITE-style local active-site matches."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

# Maps the standard three-letter amino acid code (as stored in a PDB file,
# e.g. "ARG") to the one-letter code used in Foldseek alignment strings
# (e.g. "R"). Do NOT just take the first letter of the three-letter code --
# for many residues (ARG, ASP, GLU, GLN, LYS, PHE, TRP, TYR, ASN...) the
# first letter is not the correct one-letter code.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y",
}


@dataclass(frozen=True)
class SitePair:
    """One experimentally annotated homolog residue mapped to the query."""

    homolog_chain: str
    homolog_resnum: int
    homolog_insertion: str
    homolog_resname: str
    query_chain: str
    query_resnum: int
    query_insertion: str
    query_resname: str
    exact: bool


def _gap(a: int, b: int) -> int:
    return max(0, int(b) - int(a) - 1)


def choose_local_triplet(pairs: Iterable[SitePair]) -> list[SitePair]:
    """Choose exactly one compact 3-residue cluster from one homolog site.

    The cluster is minimized by query and homolog numbering gaps first, then
    by span, and finally prefers conserved amino-acid identities. This makes
    the result deterministic and avoids displaying an entire active site.
    """
    ordered = sorted(
        list(pairs),
        key=lambda p: (
            p.query_resnum,
            p.query_insertion,
            p.homolog_resnum,
            p.homolog_insertion,
        ),
    )
    if len(ordered) <= 3:
        return ordered

    best: tuple[tuple[int, int, int, int, int], list[SitePair]] | None = None
    for i in range(len(ordered) - 2):
        group = ordered[i : i + 3]
        query_gap = sum(_gap(group[j].query_resnum, group[j + 1].query_resnum) for j in range(2))
        homolog_gap = sum(_gap(group[j].homolog_resnum, group[j + 1].homolog_resnum) for j in range(2))
        query_span = group[-1].query_resnum - group[0].query_resnum
        homolog_span = group[-1].homolog_resnum - group[0].homolog_resnum
        exact = sum(1 for p in group if p.exact)
        score = (query_gap, homolog_gap, query_span, homolog_span, -exact)
        if best is None or score < best[0]:
            best = (score, group)
    return best[1] if best else []


def correspondence_is_valid(
    query_residue: Mapping,
    homolog_residue: Mapping,
    aligned_query: str,
    aligned_homolog: str,
) -> bool:
    """Validate a non-gap alignment column against actual residue records."""
    if aligned_query == "-" or aligned_homolog == "-":
        return False
    qname = str(query_residue.get("resname", "")).upper()
    hname = str(homolog_residue.get("resname", "")).upper()
    q_one = THREE_TO_ONE.get(qname)
    h_one = THREE_TO_ONE.get(hname)
    aligned_query = aligned_query.upper()
    aligned_homolog = aligned_homolog.upper()
    q_ok = q_one == aligned_query or aligned_query == "X"
    h_ok = h_one == aligned_homolog or aligned_homolog == "X"
    return bool(q_one and h_one and q_ok and h_ok)


def format_pair(pair: SitePair) -> str:
    """Human-readable SPRITE-style correspondence."""
    h_ins = pair.homolog_insertion or ""
    q_ins = pair.query_insertion or ""
    check = " ✓" if pair.exact else ""
    return (
        f"{pair.homolog_chain}{pair.homolog_resnum}{h_ins} {pair.homolog_resname} "
        f"matches {pair.query_chain}{pair.query_resnum}{q_ins} {pair.query_resname}{check}"
    )
