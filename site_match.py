"""Canonical, testable logic for SPRITE-style local active-site matches."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence


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
    return bool(qname and hname and qname[0] == aligned_query.upper() and hname[0] == aligned_homolog.upper())


def format_pair(pair: SitePair) -> str:
    """Human-readable SPRITE-style correspondence."""
    h_ins = pair.homolog_insertion or ""
    q_ins = pair.query_insertion or ""
    check = " ✓" if pair.exact else ""
    return (
        f"{pair.homolog_chain}{pair.homolog_resnum}{h_ins} {pair.homolog_resname} "
        f"matches {pair.query_chain}{pair.query_resnum}{q_ins} {pair.query_resname}{check}"
    )
