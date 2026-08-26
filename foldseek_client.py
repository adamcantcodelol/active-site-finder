"""
Foldseek client for the Active Site Finder.

The app deliberately uses the Foldseek web server in TM-align mode because
the primary ranking criterion is the closest structural match by RMSD.
"""
from __future__ import annotations

import dataclasses
import time
from dataclasses import dataclass, field
from typing import Any, Callable, Optional, Sequence

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
    seq_identity: Optional[float]  # normalized fraction, 0..1
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
        return self.target_id.split("_", 1)[0][:4].upper()

    @property
    def chain_id(self) -> Optional[str]:
        return self.target_id.split("_", 1)[1] if "_" in self.target_id else None

    def as_table_row(self) -> dict:
        return {
            "PDB / Chain": self.target_id,
            "Description": self.description or "n/a",
            "RMSD (Å)": _fmt(self.rmsd),
            "TM-Score": _fmt(self.tm_score),
            "Seq. Identity": _fmt(
                self.seq_identity * 100 if self.seq_identity is not None else None,
                suffix="%",
            ),
            "E-value": _fmt(self.e_value, sci=True),
        }


def _fmt(value: Optional[float], sci: bool = False, suffix: str = "") -> str:
    if value is None:
        return "n/a"
    return f"{value:.2e}" if sci else f"{value:.2f}{suffix}"


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


def _normalize_seq_identity(value: Any) -> Optional[float]:
    value = _to_float(value)
    if value is None:
        return None
    # Foldseek's fident is a fraction; tolerate percentage-style API values too.
    if value > 1:
        value /= 100.0
    return max(0.0, min(1.0, value))


def _hit_from_raw(raw: dict, database: str) -> Hit:
    target = str(_first_present(raw, "target", "id") or "unknown").strip()
    return Hit(
        target_id=target,
        description=str(
            _first_present(raw, "tDescription", "description", "theader", "header") or ""
        ),
        database=database,
        e_value=_to_float(_first_present(raw, "eval", "evalue", "e_value")),
        seq_identity=_normalize_seq_identity(
            _first_present(raw, "fident", "seqId", "seqid")
        ),
        score=_to_float(_first_present(raw, "score", "bits")),
        # Foldseek distinguishes TM-score from `prob` (homology probability).
        tm_score=_to_float(
            _first_present(raw, "alntmscore", "tmScore", "tmscore", "tm_score")
        ),
        rmsd=_to_float(_first_present(raw, "rmsd", "rmsD", "RMSD")),
        q_start=_to_int(_first_present(raw, "qStartPos", "qstart", "qStart")),
        q_end=_to_int(_first_present(raw, "qEndPos", "qend", "qEnd")),
        t_start=_to_int(_first_present(raw, "dbStartPos", "tstart", "tStart")),
        t_end=_to_int(_first_present(raw, "dbEndPos", "tend", "tEnd")),
        # Official Foldseek alignment field names are qaln and taln.
        q_aln=_first_present(raw, "qaln", "qAln", "qAlignedSeq"),
        t_aln=_first_present(raw, "taln", "tAln", "dbAln", "tAlignedSeq"),
        raw=raw,
    )


def _to_int(value: Any) -> Optional[int]:
    try:
        return None if value is None else int(value)
    except (TypeError, ValueError):
        return None


def submit_search(
    structure_text: str,
    filename: str = "query.pdb",
    mode: str = "tmalign",
    databases: Sequence[str] = ("pdb100",),
    timeout: int = 30,
) -> str:
    if mode not in SEARCH_MODES:
        raise ValueError(f"mode must be one of {SEARCH_MODES}, got {mode!r}")
    if not databases:
        raise ValueError("At least one database must be selected.")

    try:
        resp = requests.post(
            f"{FOLDSEEK_API_BASE}/ticket",
            files={"q": (filename, structure_text)},
            data={"mode": mode, "database[]": list(databases)},
            timeout=timeout,
        )
        resp.raise_for_status()
        payload = resp.json()
    except requests.exceptions.RequestException as exc:
        raise FoldseekError(f"Could not reach Foldseek API: {exc}") from exc
    except ValueError as exc:
        raise FoldseekError("Foldseek returned invalid JSON while creating the job.") from exc

    ticket_id = payload.get("id")
    if not ticket_id:
        raise FoldseekError(
            f"Foldseek did not return a ticket id ({payload.get('reason', 'unknown error')})."
        )
    return str(ticket_id)


def submit_search_by_pdb_id(
    pdb_id: str,
    mode: str = "tmalign",
    databases: Sequence[str] = ("pdb100",),
) -> tuple[str, str]:
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
    ticket_id = submit_search(
        structure_text,
        filename=f"{pdb_id}.pdb",
        mode=mode,
        databases=databases,
    )
    return ticket_id, structure_text


def get_ticket_status(ticket_id: str, timeout: int = 15) -> str:
    try:
        resp = requests.get(f"{FOLDSEEK_API_BASE}/ticket/{ticket_id}", timeout=timeout)
        resp.raise_for_status()
        payload = resp.json()
    except requests.exceptions.RequestException as exc:
        raise FoldseekError(f"Could not reach Foldseek API: {exc}") from exc
    except ValueError as exc:
        raise FoldseekError("Foldseek returned invalid JSON for job status.") from exc
    return str(payload.get("status", "UNKNOWN")).upper()


def poll_until_complete(
    ticket_id: str,
    max_wait_seconds: int = 300,
    poll_interval_seconds: int = 3,
    on_status: Optional[Callable[[str, int], None]] = None,
) -> None:
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
            raise FoldseekError(
                f"Timed out after {elapsed}s waiting for ticket {ticket_id}."
            )

        time.sleep(min(poll_interval_seconds, max_wait_seconds - elapsed))


def fetch_results(
    ticket_id: str,
    databases: Sequence[str],
    timeout: int = 30,
) -> list[Hit]:
    all_hits: list[Hit] = []

    for db_index, db_name in enumerate(databases):
        try:
            resp = requests.get(
                f"{FOLDSEEK_API_BASE}/result/{ticket_id}/{db_index}",
                timeout=timeout,
            )
            if resp.status_code == 404:
                continue
            resp.raise_for_status()
            payload = resp.json()
        except requests.exceptions.RequestException as exc:
            raise FoldseekError(
                f"Could not fetch results for '{db_name}': {exc}"
            ) from exc
        except ValueError as exc:
            raise FoldseekError(
                f"Foldseek returned invalid JSON for '{db_name}'."
            ) from exc

        for db_result in payload.get("results", []):
            if not isinstance(db_result, dict):
                continue
            reported_db = str(db_result.get("db") or db_name)
            for alignment_group in db_result.get("alignments", []):
                records = (
                    alignment_group
                    if isinstance(alignment_group, list)
                    else [alignment_group]
                )
                for raw in records:
                    if isinstance(raw, dict):
                        all_hits.append(_hit_from_raw(raw, reported_db))

    # The app's primary goal is the closest structural match by RMSD.
    # Put hits with a known RMSD first; lower RMSD is better.
    # TM-score is a useful tie-breaker.
    all_hits.sort(
        key=lambda h: (
            h.rmsd is None,
            h.rmsd if h.rmsd is not None else float("inf"),
            -(h.tm_score if h.tm_score is not None else -1.0),
            h.e_value if h.e_value is not None else float("inf"),
        )
    )
    return all_hits


def fetch_target_structure(hit: Hit, timeout: int = 15) -> Optional[str]:
    try:
        resp = requests.get(
            RCSB_DOWNLOAD_URL.format(pdb_id=hit.pdb_id),
            timeout=timeout,
        )
        if resp.status_code == 200:
            return resp.text
    except requests.exceptions.RequestException:
        pass
    return None


def hits_to_dicts(hits: Sequence[Hit]) -> list[dict]:
    return [dataclasses.asdict(hit) for hit in hits]
