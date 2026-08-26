"""Active-site prediction by structural homology.

The closest structural matches are ranked by RMSD. Experimentally annotated
PDB SITE residues are preferred; when a PDB has no SITE records, the code falls
back to residues making atomic contact with a non-water HETATM ligand. Those
template residues are then mapped through the Foldseek alignment onto the
query structure.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from foldseek_client import Hit, fetch_target_structure


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y",
}


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


def _safe_int(value: str) -> Optional[int]:
    try:
        return int(value.strip())
    except (TypeError, ValueError):
        return None


def parse_site_residues(pdb_text: str, chain_id: Optional[str] = None) -> list[dict]:
    """Parse standard PDB SITE records using the official column layout."""
    residues: list[dict] = []
    seen: set[tuple[str, int, str]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith("SITE"):
            continue
        for offset in (18, 29, 40, 51):
            if len(line) < offset + 10:
                continue
            resname = line[offset:offset + 3].strip()
            chain = line[offset + 4:offset + 5].strip()
            resnum = _safe_int(line[offset + 5:offset + 9])
            insertion = line[offset + 9:offset + 10].strip()
            if not resname or resnum is None:
                continue
            if chain_id is not None and chain != chain_id:
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


def parse_ca_residues(pdb_text: str, chain_id: Optional[str] = None) -> list[dict]:
    """Return ordered protein CA residues for a chain."""
    residues: list[dict] = []
    seen: set[tuple[int, str]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        chain = line[21].strip()
        if chain_id is not None and chain != chain_id:
            continue
        resnum = _safe_int(line[22:26])
        if resnum is None:
            continue
        insertion = line[26].strip()
        resname = line[17:20].strip().upper()
        key = (resnum, insertion)
        if key in seen:
            continue
        seen.add(key)
        residues.append({
            "resnum": resnum,
            "insertion_code": insertion,
            "resname": resname,
            # foldseek_client._alignment_start_score uses this one-letter code
            # to anchor an alignment to the actual modeled CA residues.
            "one": THREE_TO_ONE.get(resname, "X"),
        })
    return residues


def parse_ca_residue_names(pdb_text: str, chain_id: Optional[str] = None) -> dict[int, str]:
    return {r["resnum"]: r["resname"] for r in parse_ca_residues(pdb_text, chain_id)}


def guess_first_chain_id(pdb_text: str) -> Optional[str]:
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            return line[21].strip()
    return None


def _parse_atom(line: str) -> Optional[dict]:
    if len(line) < 54:
        return None
    try:
        return {
            "atom": line[12:16].strip(),
            "resname": line[17:20].strip(),
            "chain": line[21].strip(),
            "resnum": int(line[22:26].strip()),
            "icode": line[26].strip(),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
        }
    except (TypeError, ValueError):
        return None


def parse_ligand_contact_residues(
    pdb_text: str,
    chain_id: Optional[str] = None,
    cutoff_angstrom: float = 4.0,
) -> list[dict]:
    """Infer binding residues from atomic contacts with non-water HETATM ligands."""
    excluded = {
        "HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "ACT",
        "CL", "NA", "K", "CA", "MG", "MN", "ZN", "FE", "CO", "NI",
    }
    ligand_atoms: list[dict] = []
    protein_atoms: list[dict] = []
    for line in pdb_text.splitlines():
        if line.startswith("HETATM"):
            atom = _parse_atom(line)
            if atom and atom["resname"].upper() not in excluded:
                ligand_atoms.append(atom)
        elif line.startswith("ATOM"):
            atom = _parse_atom(line)
            if atom and (chain_id is None or atom["chain"] == chain_id):
                protein_atoms.append(atom)
    if not ligand_atoms or not protein_atoms:
        return []
    cutoff2 = cutoff_angstrom * cutoff_angstrom
    contacts: dict[tuple[str, int, str], dict] = {}
    for atom in protein_atoms:
        key = (atom["chain"], atom["resnum"], atom["icode"])
        if atom["atom"].upper().startswith("H"):
            continue
        for lig in ligand_atoms:
            dx = atom["x"] - lig["x"]
            dy = atom["y"] - lig["y"]
            dz = atom["z"] - lig["z"]
            if dx * dx + dy * dy + dz * dz <= cutoff2:
                contacts[key] = {
                    "chain_id": atom["chain"],
                    "resnum": atom["resnum"],
                    "insertion_code": atom["icode"],
                    "resname": atom["resname"],
                }
                break
    return list(contacts.values())


def _template_site_residues(pdb_text: str, chain_id: Optional[str]) -> list[dict]:
    site = parse_site_residues(pdb_text, chain_id)
    if site:
        return site
    return parse_ligand_contact_residues(pdb_text, chain_id)


def _map_binding_site_to_query(
    hit: Hit,
    target_structure: str,
    query_structure: str,
    query_chain_id: Optional[str],
) -> list[tuple[int, str, str]]:
    """Map template binding residues through a Foldseek alignment."""
    if not hit.q_aln or not hit.t_aln or hit.q_start is None or hit.t_start is None:
        return []
    target_residues = parse_ca_residues(target_structure, hit.chain_id)
    query_residues = parse_ca_residues(query_structure, query_chain_id)
    site_residues = _template_site_residues(target_structure, hit.chain_id)
    if not target_residues or not query_residues or not site_residues:
        return []
    site_keys = {(r["resnum"], r["insertion_code"]) for r in site_residues}
    q_index = hit.q_start - 1
    t_index = hit.t_start - 1
    mapped: list[tuple[int, str, str]] = []
    for q_char, t_char in zip(hit.q_aln, hit.t_aln):
        q_present = q_char != "-"
        t_present = t_char != "-"
        if q_present and t_present and 0 <= t_index < len(target_residues) and 0 <= q_index < len(query_residues):
            target_res = target_residues[t_index]
            if (target_res["resnum"], target_res["insertion_code"]) in site_keys:
                query_res = query_residues[q_index]
                mapped.append((query_res["resnum"], query_res["insertion_code"], query_res["resname"]))
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
    """Predict the active site using the closest experimental PDB templates."""
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
        mapped = _map_binding_site_to_query(hit, target_structure, query_pdb_text, query_chain_id)
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
