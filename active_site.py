"""Active-site prediction by structural homology.

Experimental PDB SITE records are preferred. If a structure has no SITE
records, ligand-contact residues are used as a conservative fallback. All
residue mappings are validated against the actual deposited C-alpha sequence
before they are exposed to the UI.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from foldseek_client import Hit, fetch_target_structure

THREE_TO_ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q",
    "GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K",
    "MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
    "TYR":"Y","VAL":"V","MSE":"M","SEP":"S","TPO":"T","PTR":"Y",
}

@dataclass
class ActiveSiteResidue:
    query_resnum: int
    query_resname: Optional[str]
    support_count: int
    supporting_hits: list[str] = field(default_factory=list)
    insertion_code: str = ""
    evidence_type: str = "SITE"

    @property
    def display_resnum(self) -> str:
        return f"{self.query_resnum}{self.insertion_code}".strip()


def _safe_int(value: str) -> Optional[int]:
    try:
        return int(value.strip())
    except (TypeError, ValueError):
        return None


def parse_site_residues(pdb_text: str, chain_id: Optional[str] = None) -> list[dict]:
    residues: list[dict] = []
    seen: set[tuple[str, int, str]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith("SITE"):
            continue
        for offset in (18, 29, 40, 51):
            if len(line) < offset + 9:
                continue
            resname = line[offset:offset + 3].strip().upper()
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
            residues.append({"chain_id": chain, "resnum": resnum, "insertion_code": insertion, "resname": resname})
    return residues


def parse_ca_residues(pdb_text: str, chain_id: Optional[str] = None) -> list[dict]:
    residues: list[dict] = []
    seen: set[tuple[int, str]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
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
        residues.append({"resnum": resnum, "insertion_code": insertion, "resname": resname, "one": THREE_TO_ONE.get(resname, "X")})
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
        return {"atom": line[12:16].strip(), "resname": line[17:20].strip(), "chain": line[21].strip(), "resnum": int(line[22:26].strip()), "icode": line[26].strip(), "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54])}
    except (TypeError, ValueError):
        return None


def parse_ligand_contact_residues(pdb_text: str, chain_id: Optional[str] = None, cutoff_angstrom: float = 4.0) -> list[dict]:
    excluded = {"HOH","WAT","DOD","SO4","PO4","GOL","EDO","PEG","ACT","CL","NA","K","CA","MG","MN","ZN","FE","CO","NI"}
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
        if atom["atom"].upper().startswith("H"):
            continue
        key = (atom["chain"], atom["resnum"], atom["icode"])
        for lig in ligand_atoms:
            dx, dy, dz = atom["x"] - lig["x"], atom["y"] - lig["y"], atom["z"] - lig["z"]
            if dx * dx + dy * dy + dz * dz <= cutoff2:
                contacts[key] = {"chain_id": atom["chain"], "resnum": atom["resnum"], "insertion_code": atom["icode"], "resname": atom["resname"]}
                break
    return list(contacts.values())


def _template_site_residues(pdb_text: str, chain_id: Optional[str]) -> tuple[list[dict], str]:
    site = parse_site_residues(pdb_text, chain_id)
    if site:
        return site, "SITE"
    contacts = parse_ligand_contact_residues(pdb_text, chain_id)
    return contacts, "LIGAND_CONTACT"


def _anchor_alignment(alignment: str, residues: list[dict], start_hint: Optional[int]) -> tuple[int, float]:
    if not alignment or not residues:
        return -1, 0.0
    ungapped = [c.upper() for c in alignment if c != "-"]
    if not ungapped:
        return -1, 0.0
    window = min(len(ungapped), 80)
    probe = ungapped[:window]
    hint = max(0, (start_hint or 1) - 1)
    best_idx, best_score = -1, -1e9
    for idx in range(len(residues)):
        if idx + len(probe) > len(residues):
            break
        matches = sum(residues[idx + j]["one"] == aa or aa == "X" for j, aa in enumerate(probe))
        distance_penalty = min(abs(idx - hint), 200) * 0.002
        score = matches - distance_penalty
        if score > best_score:
            best_idx, best_score = idx, score
    identity = max(0.0, best_score / max(1, window))
    if identity < 0.60:
        return -1, identity
    return best_idx, identity


def map_binding_site_details(hit: Hit, target_structure: str, query_structure: str, query_chain_id: Optional[str]) -> list[dict]:
    if not hit.q_aln or not hit.t_aln:
        return []
    target = parse_ca_residues(target_structure, hit.chain_id)
    query = parse_ca_residues(query_structure, query_chain_id)
    site, evidence_type = _template_site_residues(target_structure, hit.chain_id)
    if not target or not query or not site:
        return []
    ti, t_conf = _anchor_alignment(hit.t_aln, target, hit.t_start)
    qi, q_conf = _anchor_alignment(hit.q_aln, query, hit.q_start)
    if ti < 0 or qi < 0 or min(t_conf, q_conf) < 0.60:
        return []
    site_keys = {(r["resnum"], r["insertion_code"]) for r in site}
    out: list[dict] = []
    for qc, tc in zip(hit.q_aln, hit.t_aln):
        q_present, t_present = qc != "-", tc != "-"
        if q_present and t_present and 0 <= qi < len(query) and 0 <= ti < len(target):
            qr, tr = query[qi], target[ti]
            q_ok = qr["one"] == qc.upper() or qc.upper() == "X"
            t_ok = tr["one"] == tc.upper() or tc.upper() == "X"
            if q_ok and t_ok and (tr["resnum"], tr["insertion_code"]) in site_keys:
                out.append({"tchain": hit.chain_id or "?", "tn": tr["resnum"], "ticode": tr["insertion_code"], "tname": tr["resname"], "qchain": query_chain_id or "?", "qn": qr["resnum"], "qicode": qr["insertion_code"], "qname": qr["resname"], "exact": tr["one"] == qr["one"], "evidence_type": evidence_type})
        if q_present:
            qi += 1
        if t_present:
            ti += 1
    return out


def select_local_triplet(pairs: list[dict]) -> list[dict]:
    """Select one compact, alignment-contiguous three-residue local cluster."""
    ordered = sorted(pairs, key=lambda p: (p["qn"], p.get("qicode", ""), p["tn"], p.get("ticode", "")))
    if len(ordered) <= 3:
        return ordered
    best = None
    for i in range(len(ordered) - 2):
        group = ordered[i:i + 3]
        q_gaps = sum(max(0, group[j + 1]["qn"] - group[j]["qn"] - 1) for j in range(2))
        t_gaps = sum(max(0, group[j + 1]["tn"] - group[j]["tn"] - 1) for j in range(2))
        span = (group[-1]["qn"] - group[0]["qn"]) + (group[-1]["tn"] - group[0]["tn"])
        exact = sum(1 for p in group if p.get("exact"))
        score = (q_gaps + t_gaps, span, -exact, group[0]["qn"], group[0]["tn"])
        if best is None or score < best[0]:
            best = (score, group)
    return best[1] if best else ordered[:3]


def get_sprite_match(hit: Hit, query_structure: str, query_chain_id: Optional[str]) -> list[dict]:
    """Return the single validated three-residue SPRITE-style match for one hit."""
    target_structure = fetch_target_structure(hit)
    if not target_structure:
        return []
    mapped = map_binding_site_details(hit, target_structure, query_structure, query_chain_id)
    return select_local_triplet(mapped)


def best_site_bearing_hit(hits: list[Hit], query_structure: str, query_chain_id: Optional[str], max_hits: int = 50) -> tuple[Optional[Hit], list[dict]]:
    """Find the lowest-RMSD hit with a validated, displayable local site."""
    candidates = [h for h in hits if h.rmsd is not None and h.q_aln and h.t_aln]
    candidates.sort(key=lambda h: (h.rmsd, -(h.tm_score if h.tm_score is not None else -1.0), h.e_value if h.e_value is not None else float("inf")))
    for hit in candidates[:max(1, max_hits)]:
        match = get_sprite_match(hit, query_structure, query_chain_id)
        if len(match) == 3:
            return hit, match
    return None, []


def _map_binding_site_to_query(hit: Hit, target_structure: str, query_structure: str, query_chain_id: Optional[str]) -> list[tuple[int, str, str]]:
    return [(p["qn"], p["qicode"], p["qname"]) for p in map_binding_site_details(hit, target_structure, query_structure, query_chain_id)]


def _set_sprite_status(hit: Hit, available: bool) -> None:
    """Attach UI-safe status to a hit without changing its scientific identity."""
    setattr(hit, "sprite_available", bool(available))
    base = (hit.description or "structural homolog").strip()
    marker = " · SPRITE: available" if available else " · SPRITE: no validated local site"
    for old in (" · SPRITE: available", " · SPRITE: no validated local site"):
        if base.endswith(old):
            base = base[:-len(old)]
    hit.description = base + marker


def predict_active_site(hits: list[Hit], query_pdb_text: str, query_chain_id: Optional[str] = None, top_n_hits: int = 15) -> list[ActiveSiteResidue]:
    candidates = [h for h in hits if "pdb" in (h.database or "").lower() and h.q_aln and h.t_aln and h.rmsd is not None]
    candidates.sort(key=lambda h: (h.rmsd, -(h.tm_score if h.tm_score is not None else -1.0), h.e_value if h.e_value is not None else float("inf")))

    # Evaluate every returned structural hit once so the UI can tell the user
    # which choices actually have a valid SPRITE-style local match.
    for hit in candidates:
        try:
            target = fetch_target_structure(hit)
            mapped = map_binding_site_details(hit, target, query_pdb_text, query_chain_id) if target else []
            _set_sprite_status(hit, len(select_local_triplet(mapped)) == 3)
        except Exception:
            _set_sprite_status(hit, False)

    votes: dict[tuple[int, str], list[str]] = {}
    names: dict[tuple[int, str], str] = {}
    for hit in candidates[:max(1, top_n_hits)]:
        try:
            target = fetch_target_structure(hit)
            if not target:
                continue
            for p in map_binding_site_details(hit, target, query_pdb_text, query_chain_id):
                key = (p["qn"], p["qicode"])
                supporters = votes.setdefault(key, [])
                if hit.target_id not in supporters:
                    supporters.append(hit.target_id)
                names[key] = p["qname"]
        except Exception:
            continue
    results = [ActiveSiteResidue(key[0], names.get(key), len(supporters), supporters, key[1]) for key, supporters in votes.items()]
    results.sort(key=lambda r: (-r.support_count, r.query_resnum, r.insertion_code))
    return results


def annotate_residue_names(residues: list[ActiveSiteResidue], query_pdb_text: str, chain_id: Optional[str] = None) -> None:
    names = parse_ca_residue_names(query_pdb_text, chain_id)
    for residue in residues:
        residue.query_resname = names.get(residue.query_resnum)
