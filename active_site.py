"""
Active-site prediction by structural homology.

The prediction is intentionally based on the closest-RMSD experimental
structures first. For each template with a known ligand-binding site, the
Foldseek alignment is transferred onto the query structure. Multiple
independent templates voting for the same query residue increase confidence.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import requests

from foldseek_client import Hit, fetch_target_structure

PDBE_BINDING_SITES_URL = (
    "https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/{pdb_id}"
)


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


def _first_present(d: dict, *keys: str):
    for key in keys:
        if d.get(key) is not None:
            return d[key]
    return None


def get_binding_site_residues(
    pdb_id: str,
    chain_id: Optional[str] = None,
    timeout: int = 15,
) -> list[dict]:
    """Return PDBe ligand-binding residues, preserving chain/insertion data."""
    try:
        resp = requests.get(
            PDBE_BINDING_SITES_URL.format(pdb_id=pdb_id.lower()),
            timeout=timeout,
        )
        if resp.status_code != 200:
            return []
        payload = resp.json()
    except (requests.RequestException, ValueError):
        return []

    records = payload.get(pdb_id.lower(), [])
    residues: list[dict] = []

    for record in records:
        for res in record.get("site_residues", []):
            resnum = _first_present(
                res, "author_residue_number", "residue_number"
            )
            chain = _first_present(res, "chain_id", "auth_asym_id")
            insertion = str(
                _first_present(
                    res,
                    "author_insertion_code",
                    "insertion_code",
                    "pdbx_PDB_ins_code",
                )
                or ""
            ).strip()

            if resnum is None:
                continue
            if chain_id and chain and str(chain) != str(chain_id):
                continue

            try:
                resnum = int(resnum)
            except (TypeError, ValueError):
                continue

            residues.append(
                {
                    "chain_id": str(chain or ""),
                    "author_residue_number": resnum,
                    "insertion_code": insertion,
                }
            )

    return residues


def parse_ca_residues(
    pdb_text: str,
    chain_id: Optional[str] = None,
) -> list[dict]:
    """Return ordered CA residues as they occur in a PDB chain.

    Foldseek's qstart/tstart are sequence positions. Mapping those positions
    through the actual PDB CA list avoids assuming PDB author numbering is
    contiguous (important for insertion codes and numbering gaps).
    """
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

        residues.append(
            {
                "resnum": resnum,
                "insertion_code": insertion,
                "resname": line[17:20].strip(),
            }
        )

    return residues


def parse_ca_residue_names(
    pdb_text: str, chain_id: Optional[str] = None
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


def _site_keys(site_residues: list[dict]) -> set[tuple[int, str]]:
    return {
        (r["author_residue_number"], r.get("insertion_code", ""))
        for r in site_residues
    }


def _map_binding_site_to_query(
    hit: Hit,
    target_structure: str,
    query_structure: str,
    query_chain_id: Optional[str],
) -> list[tuple[int, str, str]]:
    """Map template binding residues through the Foldseek alignment.

    Returns (query residue number, insertion code, residue name).
    """
    if not hit.q_aln or not hit.t_aln:
        return []
    if hit.q_start is None or hit.t_start is None:
        return []

    target_residues = parse_ca_residues(target_structure, hit.chain_id)
    query_residues = parse_ca_residues(query_structure, query_chain_id)

    if not target_residues or not query_residues:
        return []

    site_residues = get_binding_site_residues(
        hit.pdb_id, chain_id=hit.chain_id
    )
    site_keys = _site_keys(site_residues)
    if not site_keys:
        return []

    # Foldseek positions are 1-indexed sequence positions.
    q_index = hit.q_start - 1
    t_index = hit.t_start - 1
    mapped: list[tuple[int, str, str]] = []

    for q_char, t_char in zip(hit.q_aln, hit.t_aln):
        q_present = q_char != "-"
        t_present = t_char != "-"

        if q_present and t_present:
            if 0 <= t_index < len(target_residues):
                target_res = target_residues[t_index]
                target_key = (
                    target_res["resnum"],
                    target_res["insertion_code"],
                )
                if target_key in site_keys and 0 <= q_index < len(query_residues):
                    query_res = query_residues[q_index]
                    mapped.append(
                        (
                            query_res["resnum"],
                            query_res["insertion_code"],
                            query_res["resname"],
                        )
                    )

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
    """Predict the query active site using the closest RMSD PDB templates.

    Templates without RMSD are placed after templates with known RMSD.
    Only experimental PDB hits are used because their binding-site evidence
    comes from experimentally observed ligands.
    """
    candidates = [
        h
        for h in hits
        if "pdb" in (h.database or "").lower()
        and h.q_aln
        and h.t_aln
    ]
    candidates.sort(
        key=lambda h: (
            h.rmsd is None,
            h.rmsd if h.rmsd is not None else float("inf"),
            -(h.tm_score if h.tm_score is not None else -1.0),
        )
    )
    candidates = candidates[:top_n_hits]

    votes: dict[tuple[int, str], list[str]] = {}
    names: dict[tuple[int, str], str] = {}

    for hit in candidates:
        target_structure = fetch_target_structure(hit)
        if not target_structure:
            continue

        mapped = _map_binding_site_to_query(
            hit,
            target_structure,
            query_pdb_text,
            query_chain_id,
        )

        for resnum, insertion, resname in mapped:
            key = (resnum, insertion)
            if hit.target_id not in votes.setdefault(key, []):
                votes[key].append(hit.target_id)
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

    results.sort(
        key=lambda r: (-r.support_count, r.query_resnum, r.insertion_code)
    )
    return results


# Backward-compatible helper for callers that already have a predicted list.
def annotate_residue_names(
    residues: list[ActiveSiteResidue],
    query_pdb_text: str,
    chain_id: Optional[str] = None,
) -> None:
    names = parse_ca_residue_names(query_pdb_text, chain_id=chain_id)
    for residue in residues:
        residue.query_resname = names.get(residue.query_resnum)
