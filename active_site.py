"""Active-site prediction by structural homology.

The prediction is based on the closest experimental PDB structures by RMSD.
For each template, experimentally annotated PDB SITE records are mapped through
Foldseek's structural alignment onto the query structure. Multiple templates
supporting the same query residue increase confidence.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from foldseek_client import Hit, fetch_target_structure


@dataclass
class ActiveSiteResidue:
    query_resnum: int
    query_resname: Optional[str]
    support_count: int
    supporting_hits: list[str] = field(default_factory=list)
    insertion_code: str = ""

    @property
    def display_resnum(self) -> str:
        return f"{self.query_resnum}{self.insertion_code}".strip()


def parse_site_residues(
    pdb_text: str,
    chain_id: Optional[str] = None,
) -> list[dict]:
    """Read binding-site residues from standard PDB SITE records.

    SITE records are deposited annotations describing residues that make up a
    named structural site. Using the downloaded PDB file avoids the retired
    PDBe whole-entry binding-sites endpoint and preserves author numbering.
    """
    residues: list[dict] = []
    seen: set[tuple[str, int, str]] = set()

    for line in pdb_text.splitlines():
        if not line.startswith("SITE") or len(line) < 27:
            continue

        # A SITE record contains four residue slots. Each slot is:
        # residue name 18-20, chain 22, residue number 23-26, insertion 27.
        for offset in (18, 29, 40, 51):
            if len(line) < offset + 9:
                continue
            resname = line[offset:offset + 3].strip()
            chain = line[offset + 3:offset + 4].strip()
            resnum_text = line[offset + 4:offset + 8].strip()
            insertion = line[offset + 8:offset + 9].strip()
            if not resname or not resnum_text:
                continue
            if chain_id is not None and chain != chain_id:
                continue
            try:
                resnum = int(resnum_text)
            except ValueError:
                continue
            key = (chain, resnum, insertion)
            if key in seen:
                continue
            seen.add(key)
            residues.append({
                "chain_id": chain,
                "resnum": resnum,
                "insertion_code": insertion,
                "resname": resname,
            })

    return residues


def parse_ca_residues(
    pdb_text: str,
    chain_id: Optional[str] = None,
) -> list[dict]:
    """Return ordered CA residues as they occur in a PDB chain."""
    residues: list[dict] = []
    seen: set[tuple[int, str]] = set()

    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[12:16].strip() != "CA":
            continue
        chain = line[21].strip()
        if chain_id is not None and chain != chain_id:
            continue
        try:
            resnum = int(line[22:26].strip())
        except ValueError:
            continue
        insertion = line[26].strip()
        key = (resnum, insertion)
        if key in seen:
            continue
        seen.add(key)
        residues.append({
            "resnum": resnum,
            "insertion_code": insertion,
            "resname": line[17:20].strip(),
        })
    return residues


def parse_ca_residue_names(
    pdb_text: str,
    chain_id: Optional[str] = None,
) -> dict[int, str]:
    return {
        r["resnum"]: r["resname"]
        for r in parse_ca_residues(pdb_text, chain_id)
    }


def guess_first_chain_id(pdb_text: str) -> Optional[str]:
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            return line[21].strip()
    return None


def _map_binding_site_to_query(
    hit: Hit,
    target_structure: str,
    query_structure: str,
    query_chain_id: Optional[str],
) -> list[tuple[int, str, str]]:
    """Map template SITE residues through a Foldseek alignment."""
    if not hit.q_aln or not hit.t_aln:
        return []
    if hit.q_start is None or hit.t_start is None:
        return []

    target_residues = parse_ca_residues(target_structure, hit.chain_id)
    query_residues = parse_ca_residues(query_structure, query_chain_id)
    if not target_residues or not query_residues:
        return []

    site_residues = parse_site_residues(target_structure, hit.chain_id)
    site_keys = {
        (r["resnum"], r["insertion_code"])
        for r in site_residues
    }
    if not site_keys:
        return []

    # Foldseek sequence positions are 1-indexed. Advance each index only when
    # that side of the alignment contains a residue.
    q_index = hit.q_start - 1
    t_index = hit.t_start - 1
    mapped: list[tuple[int, str, str]] = []

    for q_char, t_char in zip(hit.q_aln, hit.t_aln):
        q_present = q_char != "-"
        t_present = t_char != "-"

        if q_present and t_present:
            if 0 <= t_index < len(target_residues) and 0 <= q_index < len(query_residues):
                target_res = target_residues[t_index]
                target_key = (target_res["resnum"], target_res["insertion_code"])
                if target_key in site_keys:
                    query_res = query_residues[q_index]
                    mapped.append((
                        query_res["resnum"],
                        query_res["insertion_code"],
                        query_res["resname"],
                    ))

        if q_present:
            q_index += 1
        if t_present:
            t_index += 1

    return mapped


def predict_active_site(
    hits: list[Hit],
    query_pdb_text: str,
    query_chain_id: Optional[str] = None,
    top_n_hits: int = 15,
) -> list[ActiveSiteResidue]:
    """Predict the query active site from the closest experimental templates.

    Lower RMSD is the primary ranking criterion. TM-score and E-value are only
    tie-breakers. Templates with no RMSD are considered after known-RMSD hits.
    """
    candidates = [
        h for h in hits
        if "pdb" in (h.database or "").lower() and h.q_aln and h.t_aln
    ]
    candidates.sort(key=lambda h: (
        h.rmsd is None,
        h.rmsd if h.rmsd is not None else float("inf"),
        -(h.tm_score if h.tm_score is not None else -1.0),
        h.e_value if h.e_value is not None else float("inf"),
    ))
    candidates = candidates[:max(1, top_n_hits)]

    votes: dict[tuple[int, str], list[str]] = {}
    names: dict[tuple[int, str], str] = {}

    for hit in candidates:
        target_structure = fetch_target_structure(hit)
        if not target_structure:
            continue
        mapped = _map_binding_site_to_query(
            hit, target_structure, query_pdb_text, query_chain_id
        )
        for resnum, insertion, resname in mapped:
            key = (resnum, insertion)
            supporters = votes.setdefault(key, [])
            if hit.target_id not in supporters:
                supporters.append(hit.target_id)
            names[key] = resname

    results = [
        ActiveSiteResidue(
            query_resnum=resnum,
            query_resname=names.get((resnum, insertion)),
            support_count=len(supporters),
            supporting_hits=supporters,
            insertion_code=insertion,
        )
        for (resnum, insertion), supporters in votes.items()
    ]
    results.sort(key=lambda r: (-r.support_count, r.query_resnum, r.insertion_code))
    return results


def annotate_residue_names(
    residues: list[ActiveSiteResidue],
    query_pdb_text: str,
    chain_id: Optional[str] = None,
) -> None:
    names = parse_ca_residue_names(query_pdb_text, chain_id)
    for residue in residues:
        residue.query_resname = names.get(residue.query_resnum)
