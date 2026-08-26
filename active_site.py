"""
active_site.py
===============
Predicts the active site of a protein structure by homology transfer:

  1. Take the structural homologs Foldseek found (with their alignments
     to the query).
  2. For each homolog that is a real experimental PDB entry, look up
     which of ITS residues touch a bound ligand/cofactor (its known
     "binding site") using the PDBe REST API.
  3. Walk the Foldseek alignment to translate those residue numbers from
     the homolog's numbering into the QUERY's numbering.
  4. Count how many independent homologs point at each query residue.
     Residues flagged by several unrelated homologs are a strong,
     evidence-based prediction of the active site — this works even if
     the query protein itself has never been experimentally studied,
     because the evidence comes entirely from its structural relatives.

This is the same general strategy used by established homology-based
function-annotation tools (e.g. COFACTOR, ProFunc): transfer known
functional-site annotations across structural alignments and use
agreement across multiple templates as a confidence signal.

DATA SOURCE
-----------
PDBe REST API: https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/{pdb_id}
Confirmed response shape (from PDBe's own official tutorial notebook):
    data[pdb_id] -> list of binding-site records, each with a
    'site_residues' list of dicts containing at least
    'author_residue_number' (an int). Chain and residue-name keys are
    not guaranteed to use one exact spelling across all PDBe API
    versions, so those are looked up defensively below.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional
import requests

from foldseek_client import Hit

PDBE_BINDING_SITES_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/{pdb_id}"


@dataclass
class ActiveSiteResidue:
    query_resnum: int
    query_resname: Optional[str]       # 3-letter code, e.g. "HIS" (from the query PDB file)
    support_count: int                 # how many independent homologs flagged this position
    supporting_hits: list[str] = field(default_factory=list)  # e.g. ["1ABC_A", "2XYZ_B"]


def _first_present(d: dict, *keys: str):
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return None


def get_binding_site_residues(pdb_id: str, chain_id: Optional[str] = None, timeout: int = 15) -> list[dict]:
    """Fetch known ligand-binding-site residues for a real PDB entry.

    Returns a list of {"chain_id": str, "author_residue_number": int} for
    every residue PDBe reports as contacting a bound ligand/cofactor,
    optionally restricted to one chain. Returns [] if the entry has no
    recorded binding sites (e.g. an apo structure) or the lookup fails —
    both are normal, non-fatal outcomes, not errors.
    """
    try:
        resp = requests.get(
            PDBE_BINDING_SITES_URL.format(pdb_id=pdb_id.lower()), timeout=timeout
        )
    except requests.exceptions.RequestException:
        return []
    if resp.status_code != 200:
        return []

    try:
        payload = resp.json()
        records = payload.get(pdb_id.lower(), [])
    except ValueError:
        return []

    residues = []
    for record in records:
        for res in record.get("site_residues", []):
            resnum = _first_present(res, "author_residue_number", "residue_number")
            chain = _first_present(res, "chain_id", "auth_asym_id")
            if resnum is None:
                continue
            if chain_id and chain and str(chain) != str(chain_id):
                continue
            residues.append({"chain_id": chain, "author_residue_number": int(resnum)})
    return residues


def _map_target_residues_to_query(hit: Hit, target_binding_residues: set[int]) -> set[int]:
    """Walk a Foldseek alignment (q_aln vs t_aln) and translate a set of
    target residue numbers into the corresponding query residue numbers.

    ASSUMPTION / KNOWN LIMITATION: this assumes the target's numbering in
    the alignment (t_start, incrementing by 1 per non-gap character) lines
    up with the PDB "author" residue numbering used by the binding-site
    API. This holds for the common case of a single, contiguous,
    standard-numbered chain. It can drift for entries with insertion
    codes, numbering gaps, or multiple segments — treat results as a
    strong hint to verify visually, not ground truth.
    """
    if not (hit.q_aln and hit.t_aln and hit.q_start and hit.t_start):
        return set()

    mapped_query_positions = set()
    q_pos = hit.q_start
    t_pos = hit.t_start

    for q_char, t_char in zip(hit.q_aln, hit.t_aln):
        q_is_residue = q_char != "-"
        t_is_residue = t_char != "-"

        if q_is_residue and t_is_residue and t_pos in target_binding_residues:
            mapped_query_positions.add(q_pos)

        if q_is_residue:
            q_pos += 1
        if t_is_residue:
            t_pos += 1

    return mapped_query_positions


def predict_active_site(hits: list[Hit], top_n_hits: int = 15) -> list[ActiveSiteResidue]:
    """Aggregate active-site evidence across the top structural homologs.

    Only homologs from real experimental PDB entries (database "pdb100")
    are used as evidence sources, since AlphaFold-model hits have no
    experimentally observed ligands to transfer.
    """
    votes: dict[int, list[str]] = {}

    candidates = [h for h in hits if "pdb" in (h.database or "").lower()][:top_n_hits]

    for hit in candidates:
        pdb_id = hit.target_id.split("_")[0]
        chain = hit.target_id.split("_")[1] if "_" in hit.target_id else None

        site_residues = get_binding_site_residues(pdb_id, chain_id=chain)
        if not site_residues:
            continue

        target_resnums = {r["author_residue_number"] for r in site_residues}
        mapped = _map_target_residues_to_query(hit, target_resnums)

        for q_resnum in mapped:
            votes.setdefault(q_resnum, []).append(hit.target_id)

    results = [
        ActiveSiteResidue(
            query_resnum=resnum,
            query_resname=None,  # filled in later once we have the query PDB text
            support_count=len(supporters),
            supporting_hits=supporters,
        )
        for resnum, supporters in votes.items()
    ]
    results.sort(key=lambda r: -r.support_count)
    return results


def annotate_residue_names(residues: list[ActiveSiteResidue], query_pdb_text: str,
                            chain_id: Optional[str] = None) -> None:
    """Fill in `query_resname` (e.g. "HIS") for each predicted residue by
    reading the query's own PDB file. Modifies the list in place.
    """
    resname_by_num = parse_ca_residue_names(query_pdb_text, chain_id=chain_id)
    for r in residues:
        r.query_resname = resname_by_num.get(r.query_resnum)


def parse_ca_residue_names(pdb_text: str, chain_id: Optional[str] = None) -> dict[int, str]:
    """Minimal PDB ATOM-record parser: returns {residue_number: 3-letter
    residue name} for every alpha-carbon (CA) atom, optionally restricted
    to one chain. No external dependency (e.g. BioPython) required.
    """
    result: dict[int, str] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        chain = line[21].strip()
        if chain_id and chain != chain_id:
            continue
        resname = line[17:20].strip()
        try:
            resnum = int(line[22:26].strip())
        except ValueError:
            continue
        result[resnum] = resname
    return result


def guess_first_chain_id(pdb_text: str) -> Optional[str]:
    """Return the chain ID of the first ATOM record found, as a sane
    default for a single-chain query structure."""
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            return line[21].strip()
    return None
