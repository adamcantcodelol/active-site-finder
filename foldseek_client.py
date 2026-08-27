"""Foldseek client for the MBRC Active Site Finder."""
from __future__ import annotations

import dataclasses
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Optional, Sequence

import numpy as np
import requests

FOLDSEEK_API_BASE = "https://search.foldseek.com/api"
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
SEARCH_MODES = ("3diaa", "tmalign")


class FoldseekError(RuntimeError):
    pass


@dataclass
class Hit:
    target_id: str
    description: str
    database: str
    e_value: Optional[float]
    seq_identity: Optional[float]
    score: Optional[float]
    tm_score: Optional[float]
    rmsd: Optional[float]
    q_start: Optional[int]
    q_end: Optional[int]
    t_start: Optional[int]
    t_end: Optional[int]
    q_aln: Optional[str]
    t_aln: Optional[str]
    raw: dict = field(default_factory=dict, repr=False)

    @property
    def pdb_id(self) -> str:
        token = self.target_id.split()[0]
        return token.split("-", 1)[0][:4].upper()

    @property
    def chain_id(self) -> Optional[str]:
        token = self.target_id.split()[0]
        if "_" not in token:
            return None
        return token.rsplit("_", 1)[1] or None

    def as_table_row(self) -> dict:
        return {
            "PDB / Chain": self.target_id,
            "Description": self.description or "n/a",
            "RMSD (Å)": _fmt(self.rmsd),
            "TM-Score": _fmt(self.tm_score, digits=3),
            "Seq. Identity": _fmt(self.seq_identity * 100 if self.seq_identity is not None else None, suffix="%"),
            "E-value": _fmt(self.e_value, sci=True),
        }


def _fmt(value: Optional[float], sci: bool = False, suffix: str = "", digits: int = 2) -> str:
    if value is None:
        return "n/a"
    return f"{value:.2e}" if sci else f"{value:.{digits}f}{suffix}"


def _first_present(d: dict, *keys: str) -> Any:
    for key in keys:
        if d.get(key) is not None:
            return d[key]
    return None


def _to_float(value: Any) -> Optional[float]:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


def _to_int(value: Any) -> Optional[int]:
    try:
        return None if value is None else int(value)
    except (TypeError, ValueError):
        return None


def _normalize_seq_identity(value: Any) -> Optional[float]:
    value = _to_float(value)
    if value is None:
        return None
    if value > 1:
        value /= 100.0
    return max(0.0, min(1.0, value))


def _hit_from_raw(raw: dict, database: str) -> Hit:
    return Hit(
        target_id=str(_first_present(raw, "target", "id") or "unknown").strip(),
        description=str(_first_present(raw, "tDescription", "description", "theader", "header") or ""),
        database=database,
        e_value=_to_float(_first_present(raw, "eval", "evalue", "e_value")),
        seq_identity=_normalize_seq_identity(_first_present(raw, "fident", "seqId", "seqid")),
        score=_to_float(_first_present(raw, "score", "bits")),
        tm_score=_to_float(_first_present(raw, "alntmscore", "tmScore", "tmscore", "tm_score")),
        rmsd=_to_float(_first_present(raw, "rmsd", "rmsD", "RMSD")),
        q_start=_to_int(_first_present(raw, "qStartPos", "qstart", "qStart")),
        q_end=_to_int(_first_present(raw, "qEndPos", "qend", "qEnd")),
        t_start=_to_int(_first_present(raw, "dbStartPos", "tstart", "tStart")),
        t_end=_to_int(_first_present(raw, "dbEndPos", "tend", "tEnd")),
        q_aln=_first_present(raw, "qaln", "qAln", "qAlignedSeq"),
        t_aln=_first_present(raw, "taln", "tAln", "dbAln", "tAlignedSeq"),
        raw=raw,
    )


def submit_search(structure_text: str, filename: str = "query.pdb", mode: str = "tmalign", databases: Sequence[str] = ("pdb100",), timeout: int = 30) -> str:
    if mode not in SEARCH_MODES:
        raise ValueError(f"mode must be one of {SEARCH_MODES}, got {mode!r}")
    if not databases:
        raise ValueError("At least one database must be selected.")
    try:
        resp = requests.post(f"{FOLDSEEK_API_BASE}/ticket", files={"q": (filename, structure_text, "application/octet-stream")}, data=[("mode", mode), *[("database[]", db) for db in databases]], timeout=timeout)
        resp.raise_for_status()
        payload = resp.json()
    except requests.exceptions.RequestException as exc:
        raise FoldseekError(f"Could not reach Foldseek API: {exc}") from exc
    except ValueError as exc:
        raise FoldseekError("Foldseek returned invalid JSON while creating the job.") from exc
    ticket_id = payload.get("id")
    if not ticket_id:
        raise FoldseekError(f"Foldseek did not return a ticket id ({payload.get('reason', 'unknown error')}).")
    return str(ticket_id)


def submit_search_by_pdb_id(pdb_id: str, mode: str = "tmalign", databases: Sequence[str] = ("pdb100",)) -> tuple[str, str]:
    pdb_id = pdb_id.strip().upper()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        raise ValueError("PDB IDs must be exactly 4 letters/numbers, e.g. 4HHB.")
    try:
        resp = requests.get(RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id), timeout=15)
        resp.raise_for_status()
    except requests.exceptions.RequestException as exc:
        raise FoldseekError(f"Could not fetch PDB '{pdb_id}' from RCSB: {exc}") from exc
    structure_text = resp.text
    if not any(line.startswith(("ATOM", "HETATM", "MODEL", "HEADER")) for line in structure_text.splitlines()):
        raise FoldseekError(f"RCSB returned an unexpected file for PDB ID '{pdb_id}'.")
    return submit_search(structure_text, filename=f"{pdb_id}.pdb", mode=mode, databases=databases), structure_text


def get_ticket_status(ticket_id: str, timeout: int = 15) -> str:
    try:
        resp = requests.get(f"{FOLDSEEK_API_BASE}/ticket/{ticket_id}", timeout=timeout)
        resp.raise_for_status()
        return str(resp.json().get("status", "UNKNOWN")).upper()
    except requests.exceptions.RequestException as exc:
        raise FoldseekError(f"Could not reach Foldseek API: {exc}") from exc
    except ValueError as exc:
        raise FoldseekError("Foldseek returned invalid JSON for job status.") from exc


def poll_until_complete(ticket_id: str, max_wait_seconds: int = 300, poll_interval_seconds: int = 3, on_status: Optional[Callable[[str, int], None]] = None) -> None:
    if max_wait_seconds <= 0 or poll_interval_seconds <= 0:
        raise ValueError("Timeout and poll interval must be positive.")
    started = time.monotonic()
    while True:
        status = get_ticket_status(ticket_id)
        elapsed = int(time.monotonic() - started)
        if on_status:
            on_status(status, elapsed)
        if status == "COMPLETE":
            return
        if status == "ERROR":
            raise FoldseekError(f"Foldseek job {ticket_id} failed on the server.")
        if status not in {"PENDING", "RUNNING", "STARTED"}:
            raise FoldseekError(f"Unexpected Foldseek job status: {status!r}.")
        if elapsed >= max_wait_seconds:
            raise FoldseekError(f"Timed out after {elapsed}s waiting for ticket {ticket_id}.")
        time.sleep(min(poll_interval_seconds, max_wait_seconds - elapsed))


def fetch_results(ticket_id: str, databases: Sequence[str], timeout: int = 30) -> list[Hit]:
    all_hits: list[Hit] = []
    for db_index, db_name in enumerate(databases):
        try:
            resp = requests.get(f"{FOLDSEEK_API_BASE}/result/{ticket_id}/{db_index}", timeout=timeout)
            if resp.status_code == 404:
                continue
            resp.raise_for_status()
            payload = resp.json()
        except requests.exceptions.RequestException as exc:
            raise FoldseekError(f"Could not fetch results for '{db_name}': {exc}") from exc
        except ValueError as exc:
            raise FoldseekError(f"Foldseek returned invalid JSON for '{db_name}'.") from exc
        for db_result in payload.get("results", []):
            if not isinstance(db_result, dict):
                continue
            reported_db = str(db_result.get("db") or db_name)
            for alignment_group in db_result.get("alignments", []):
                records = alignment_group if isinstance(alignment_group, list) else [alignment_group]
                for raw in records:
                    if isinstance(raw, dict):
                        all_hits.append(_hit_from_raw(raw, reported_db))
    _sort_hits(all_hits)
    return all_hits


def _sort_hits(hits: list[Hit]) -> None:
    hits.sort(key=lambda h: (h.rmsd is None, h.rmsd if h.rmsd is not None else float("inf"), -(h.tm_score if h.tm_score is not None else -1.0), h.e_value if h.e_value is not None else float("inf")))


_THREE_TO_ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M","SEP":"S","TPO":"T","PTR":"Y"}


def _parse_ca_by_chain(pdb_text: str, chain_id: Optional[str]) -> list[dict]:
    residues: list[dict] = []
    seen: set[tuple[str, int, str]] = set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        chain = line[21].strip()
        if chain_id is not None and chain != chain_id:
            continue
        try:
            resnum = int(line[22:26].strip())
            icode = line[26].strip()
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)
        except ValueError:
            continue
        key = (chain, resnum, icode)
        if key in seen:
            continue
        seen.add(key)
        resname = line[17:20].strip().upper()
        residues.append({"resnum": resnum, "icode": icode, "resname": resname, "one": _THREE_TO_ONE.get(resname, "X"), "xyz": xyz})
    return residues


def _alignment_start_score(alignment: str, residues: list[dict], start_hint: Optional[int]) -> tuple[int, float]:
    """Anchor a local Foldseek alignment to the deposited, modeled C-alpha sequence.

    PDB files can omit unresolved residues, so a Foldseek sequence position is
    not necessarily a Python index. Search the modeled chain using the actual
    alignment letters. A weak anchor is rejected instead of producing a bogus
    residue correspondence.
    """
    if not alignment or not residues:
        return -1, 0.0
    probe = [c.upper() for c in alignment if c != "-"]
    if not probe:
        return -1, 0.0
    window = min(len(probe), 80)
    probe = probe[:window]
    hint = max(0, (start_hint or 1) - 1)
    best_idx, best_identity, best_score = -1, -1.0, float("-inf")
    for idx in range(len(residues) - len(probe) + 1):
        matches = sum(residues[idx + j]["one"] == aa or aa == "X" for j, aa in enumerate(probe))
        identity = matches / max(1, window)
        distance_penalty = min(abs(idx - hint), 200) * 0.002
        score = matches - distance_penalty
        if score > best_score:
            best_idx, best_identity, best_score = idx, identity, score
    if best_idx < 0 or best_identity < 0.60:
        return -1, best_identity
    return best_idx, best_identity


def _aligned_ca_pairs(hit: Hit, query_pdb: str, target_pdb: str, query_chain_id: Optional[str]) -> tuple[np.ndarray, np.ndarray]:
    if not hit.q_aln or not hit.t_aln:
        return np.empty((0, 3)), np.empty((0, 3))
    qres = _parse_ca_by_chain(query_pdb, query_chain_id)
    tres = _parse_ca_by_chain(target_pdb, hit.chain_id)
    if not qres or not tres:
        return np.empty((0, 3)), np.empty((0, 3))
    qi, q_conf = _alignment_start_score(hit.q_aln, qres, hit.q_start)
    ti, t_conf = _alignment_start_score(hit.t_aln, tres, hit.t_start)
    if qi < 0 or ti < 0 or min(q_conf, t_conf) < 0.60:
        return np.empty((0, 3)), np.empty((0, 3))
    qpts: list[np.ndarray] = []
    tpts: list[np.ndarray] = []
    for qc, tc in zip(hit.q_aln, hit.t_aln):
        q_present, t_present = qc != "-", tc != "-"
        if q_present and t_present and 0 <= qi < len(qres) and 0 <= ti < len(tres):
            qr, tr = qres[qi], tres[ti]
            if (qc.upper() == "X" or qr["one"] == qc.upper()) and (tc.upper() == "X" or tr["one"] == tc.upper()):
                qpts.append(qr["xyz"])
                tpts.append(tr["xyz"])
        if q_present:
            qi += 1
        if t_present:
            ti += 1
        if (qi >= len(qres) or ti >= len(tres)) and len(qpts) >= 3:
            break
    if len(qpts) < 3:
        return np.empty((0, 3)), np.empty((0, 3))
    return np.asarray(qpts), np.asarray(tpts)


def _kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    a0 = a - a.mean(axis=0)
    b0 = b - b.mean(axis=0)
    h = a0.T @ b0
    u, _, vt = np.linalg.svd(h)
    d = 1.0 if np.linalg.det(vt.T @ u.T) >= 0 else -1.0
    rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    aligned = a0 @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - b0) ** 2, axis=1))))


def populate_missing_rmsd(hits: list[Hit], query_pdb: str, query_chain_id: Optional[str] = None, max_hits: int = 50) -> list[Hit]:
    """Calculate C-alpha RMSD only from residue-validated alignment pairs."""
    attempted = 0
    for hit in hits:
        if hit.rmsd is not None:
            continue
        if attempted >= max_hits:
            break
        attempted += 1
        target = fetch_target_structure(hit)
        if not target:
            continue
        qpts, tpts = _aligned_ca_pairs(hit, query_pdb, target, query_chain_id)
        if len(qpts) < 3:
            continue
        try:
            hit.rmsd = _kabsch_rmsd(qpts, tpts)
        except np.linalg.LinAlgError:
            continue
    _sort_hits(hits)
    return hits


def fetch_target_structure(hit: Hit, timeout: int = 15) -> Optional[str]:
    try:
        resp = requests.get(RCSB_DOWNLOAD_URL.format(pdb_id=hit.pdb_id), timeout=timeout)
        return resp.text if resp.status_code == 200 else None
    except requests.exceptions.RequestException:
        return None


def hits_to_dicts(hits: Sequence[Hit]) -> list[dict]:
    return [dataclasses.asdict(hit) for hit in hits]
