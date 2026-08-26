"""
foldseek_client.py
===================
Framework-agnostic Python client for the public Foldseek structure search
web service (https://search.foldseek.com/api).

This module owns the full asynchronous "ticket" lifecycle:

    1. submit_search()        -> POST /api/ticket            (create a job)
    2. poll_until_complete()  -> GET  /api/ticket/{id}        (wait for it)
    3. fetch_results()        -> GET  /api/result/{id}/{n}    (read hits)

It is intentionally UI-agnostic (no Streamlit / FastAPI imports here) so it
can be reused from a Streamlit script, a FastAPI route, a CLI, or a test
suite.

-----------------------------------------------------------------------
SOURCING NOTE - please read before trusting field names blindly
-----------------------------------------------------------------------
Foldseek does not publish a formal OpenAPI/JSON schema for its web API.
The endpoints and request shape below are confirmed against:
  - the official MMseqs2-App API example referenced by Foldseek's own docs
  - the Arcadia-Science `foldseek_apiquery.py` community client
  - the ToolUniverse `FoldseekTool` implementation (which uses the exact
    `/api/result/{ticket}/{db_index}` JSON endpoint used below)

The *result* JSON's field names (`seqId`, `eval`, `qStartPos`, ...) are
confirmed. The names for TM-score / RMSD / aligned-sequence strings are
NOT formally documented anywhere public, and may differ depending on
search `mode` ("3diaa" vs "tmalign") or may be entirely absent for a given
hit. To avoid silently fabricating data, every accessor below tries a
short list of plausible key names and falls back to `None` rather than
guessing a single name and failing/crashing. If you find the real key
names by inspecting a live response (see `Hit.raw`), tighten
`_first_present()` accordingly.
-----------------------------------------------------------------------
"""

from __future__ import annotations

import time
import dataclasses
from dataclasses import dataclass, field
from typing import Optional, Sequence, Any

import requests

FOLDSEEK_API_BASE = "https://search.foldseek.com/api"
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
ALPHAFOLD_MODEL_URL = "https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"

# Databases Foldseek's public server supports as of this writing.
# (Kept as a plain list so the UI layer can render checkboxes from it.)
AVAILABLE_DATABASES = [
    "pdb100",
    "afdb50",
    "afdb-swissprot",
    "afdb-proteome",
    "mgnify_esm30",
    "gmgcl_id",
]

SEARCH_MODES = ["3diaa", "tmalign"]


class FoldseekError(RuntimeError):
    """Raised for any unrecoverable Foldseek API failure (submission,
    polling timeout, server-side ERROR status, malformed response, ...)."""


# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------

@dataclass
class Hit:
    """One structural homolog reported by Foldseek, normalized into a
    stable shape regardless of which raw JSON keys the server used."""

    target_id: str
    description: str
    database: str
    e_value: Optional[float]
    seq_identity: Optional[float]      # fraction, 0-1
    score: Optional[float]             # bit score (3diaa) or raw score
    tm_score: Optional[float]          # None if server didn't provide one
    rmsd: Optional[float]              # None if server didn't provide one
    q_start: Optional[int]
    q_end: Optional[int]
    t_start: Optional[int]
    t_end: Optional[int]
    q_aln: Optional[str]               # aligned query sequence block (with gaps)
    t_aln: Optional[str]               # aligned target sequence block (with gaps)
    raw: dict = field(default_factory=dict, repr=False)  # original record, for debugging

    def as_table_row(self) -> dict:
        """Flat dict matching the columns the UI needs to display."""
        return {
            "Target ID": self.target_id,
            "Description": self.description or "n/a",
            "Database": self.database,
            "TM-Score": _fmt(self.tm_score),
            "RMSD (\u00c5)": _fmt(self.rmsd),
            "E-value": _fmt(self.e_value, sci=True),
            "Seq. Identity": _fmt(
                self.seq_identity * 100 if self.seq_identity is not None else None,
                suffix="%",
            ),
        }


def _fmt(value: Optional[float], sci: bool = False, suffix: str = "") -> str:
    if value is None:
        return "n/a"
    if sci:
        return f"{value:.2e}"
    return f"{value:.2f}{suffix}"


def _first_present(d: dict, *keys: str) -> Any:
    """Return the first non-None value found under any of `keys`."""
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return None


def _hit_from_raw(raw: dict, database: str) -> Hit:
    """Normalize one raw alignment record from /api/result into a Hit.

    Tries several plausible key spellings for the less-well-documented
    fields (tm_score, rmsd, aligned sequence blocks) — see module
    docstring for why.
    """
    q_start = _first_present(raw, "qStartPos", "qstart", "qStart")
    q_end = _first_present(raw, "qEndPos", "qend", "qEnd")
    t_start = _first_present(raw, "dbStartPos", "tstart", "tStart")
    t_end = _first_present(raw, "dbEndPos", "tend", "tEnd")

    return Hit(
        target_id=str(_first_present(raw, "target", "id") or "unknown"),
        description=str(_first_present(raw, "tDescription", "description", "header") or ""),
        database=database,
        e_value=_to_float(_first_present(raw, "eval", "evalue", "e_value")),
        seq_identity=_to_float(_first_present(raw, "seqId", "seqid", "fident")),
        score=_to_float(_first_present(raw, "score", "bits")),
        tm_score=_to_float(
            _first_present(raw, "tmScore", "tmscore", "alntmscore", "prob")
        ),
        rmsd=_to_float(_first_present(raw, "rmsd", "rmsD", "RMSD")),
        q_start=int(q_start) if q_start is not None else None,
        q_end=int(q_end) if q_end is not None else None,
        t_start=int(t_start) if t_start is not None else None,
        t_end=int(t_end) if t_end is not None else None,
        q_aln=_first_present(raw, "qAln", "qaln", "qAlignedSeq"),
        t_aln=_first_present(raw, "dbAln", "taln", "tAlignedSeq"),
        raw=raw,
    )


def _to_float(v: Any) -> Optional[float]:
    if v is None:
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------
# 1. Submit
# --------------------------------------------------------------------------

def submit_search(
    structure_text: str,
    filename: str = "query.pdb",
    mode: str = "3diaa",
    databases: Sequence[str] = ("pdb100", "afdb50", "afdb-swissprot"),
    timeout: int = 30,
) -> str:
    """POST a structure (PDB or mmCIF text) to Foldseek and return the
    resulting ticket id.

    Parameters
    ----------
    structure_text: raw contents of a .pdb / .cif file (as text)
    filename: filename to report to the server (extension matters for
        format auto-detection; keep .pdb or .cif)
    mode: "3diaa" (fast structural-alphabet search, default) or
        "tmalign" (slower, TM-align based re-scoring)
    databases: which Foldseek-hosted databases to search
    """
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
    except requests.exceptions.RequestException as exc:
        raise FoldseekError(f"Could not reach Foldseek API: {exc}") from exc

    if resp.status_code != 200:
        raise FoldseekError(
            f"Foldseek submission failed: HTTP {resp.status_code} - {resp.text[:300]}"
        )

    payload = resp.json()
    ticket_id = payload.get("id")
    if not ticket_id:
        # Foldseek returns {"status": "...", "reason": "..."} on rate limits
        reason = payload.get("reason", "no reason given")
        raise FoldseekError(f"Foldseek did not return a ticket id ({reason}).")
    return ticket_id


def submit_search_by_pdb_id(
    pdb_id: str,
    mode: str = "3diaa",
    databases: Sequence[str] = ("pdb100", "afdb50", "afdb-swissprot"),
) -> tuple[str, str]:
    """Fetch a structure from RCSB by 4-character PDB ID, then submit it.

    Returns (ticket_id, structure_text) — the structure text is returned
    too so the caller can render the *query* in the 3D viewer without a
    second network round trip.
    """
    pdb_id = pdb_id.strip().upper()
    if len(pdb_id) != 4:
        raise ValueError("PDB IDs are 4 characters, e.g. '4HHB'.")

    try:
        resp = requests.get(RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id), timeout=15)
    except requests.exceptions.RequestException as exc:
        raise FoldseekError(f"Could not reach RCSB PDB: {exc}") from exc
    if resp.status_code != 200:
        raise FoldseekError(f"PDB ID '{pdb_id}' not found on RCSB (HTTP {resp.status_code}).")

    structure_text = resp.text
    ticket_id = submit_search(
        structure_text, filename=f"{pdb_id}.pdb", mode=mode, databases=databases
    )
    return ticket_id, structure_text


# --------------------------------------------------------------------------
# 2. Poll
# --------------------------------------------------------------------------

def get_ticket_status(ticket_id: str, timeout: int = 15) -> str:
    """Single status check. Returns one of PENDING / RUNNING / COMPLETE / ERROR."""
    try:
        resp = requests.get(f"{FOLDSEEK_API_BASE}/ticket/{ticket_id}", timeout=timeout)
    except requests.exceptions.RequestException as exc:
        raise FoldseekError(f"Could not reach Foldseek API: {exc}") from exc
    if resp.status_code != 200:
        raise FoldseekError(f"Status check failed: HTTP {resp.status_code}")
    return resp.json().get("status", "UNKNOWN")


def poll_until_complete(
    ticket_id: str,
    max_wait_seconds: int = 180,
    poll_interval_seconds: int = 3,
    on_status: Optional[callable] = None,
) -> None:
    """Block until the ticket reaches COMPLETE, raising FoldseekError on
    ERROR or timeout.

    `on_status(status: str, elapsed: int)` is called after every poll if
    provided, so a UI layer can update a progress indicator.
    """
    elapsed = 0
    while elapsed <= max_wait_seconds:
        status = get_ticket_status(ticket_id)
        if on_status:
            on_status(status, elapsed)

        if status == "COMPLETE":
            return
        if status == "ERROR":
            raise FoldseekError(
                f"Foldseek job {ticket_id} failed on the server (status=ERROR)."
            )

        time.sleep(poll_interval_seconds)
        elapsed += poll_interval_seconds

    raise FoldseekError(
        f"Timed out after {max_wait_seconds}s waiting for ticket {ticket_id}."
    )


# --------------------------------------------------------------------------
# 3. Fetch results
# --------------------------------------------------------------------------

def fetch_results(
    ticket_id: str,
    databases: Sequence[str],
    timeout: int = 30,
) -> list[Hit]:
    """Fetch and normalize alignment results for every searched database.

    Foldseek's `/api/result/{ticket}/{n}` endpoint is indexed by the
    *position* of the database in the original submission (0, 1, 2, ...),
    not by database name — hence `databases` must be passed in the same
    order used in `submit_search`.
    """
    all_hits: list[Hit] = []

    for db_index, db_name in enumerate(databases):
        try:
            resp = requests.get(
                f"{FOLDSEEK_API_BASE}/result/{ticket_id}/{db_index}", timeout=timeout
            )
        except requests.exceptions.RequestException as exc:
            raise FoldseekError(f"Could not reach Foldseek API: {exc}") from exc

        if resp.status_code == 404:
            # This db index may not have produced hits / doesn't exist — skip.
            continue
        if resp.status_code != 200:
            raise FoldseekError(
                f"Result fetch failed for db '{db_name}': HTTP {resp.status_code}"
            )

        payload = resp.json()
        for db_result in payload.get("results", []):
            reported_db = db_result.get("db", db_name)
            for alignment_group in db_result.get("alignments", []):
                # The API sometimes nests one list per query chain; flatten.
                records = alignment_group if isinstance(alignment_group, list) else [alignment_group]
                for raw in records:
                    if isinstance(raw, dict):
                        all_hits.append(_hit_from_raw(raw, reported_db))

    # Best hits first: prefer TM-score when we have it, else lower e-value.
    all_hits.sort(
        key=lambda h: (
            -(h.tm_score if h.tm_score is not None else -1),
            h.e_value if h.e_value is not None else float("inf"),
        )
    )
    return all_hits


# --------------------------------------------------------------------------
# Helper: fetch a hit's *actual* 3D coordinates for the viewer
# --------------------------------------------------------------------------

def fetch_target_structure(hit: Hit, timeout: int = 15) -> Optional[str]:
    """Best-effort fetch of the real PDB-format coordinates for a Hit, by
    going straight to the canonical source database rather than relying
    on any undocumented Foldseek structure-serving route.

    - pdb100 hits look like "1ABC_A" (pdb id + chain)      -> RCSB
    - afdb*  hits look like "AF-P12345-F1-model_v4" (or similar UniProt-
      derived ids)                                          -> AlphaFold DB
    Returns None if the source can't be determined or the fetch fails;
    callers should treat that as "no structure available to render".
    """
    target = hit.target_id
    db = (hit.database or "").lower()

    try:
        if db.startswith("afdb") or target.upper().startswith("AF-"):
            uniprot = _extract_uniprot_accession(target)
            if not uniprot:
                return None
            resp = requests.get(ALPHAFOLD_MODEL_URL.format(uniprot=uniprot), timeout=timeout)
        else:
            pdb_id = target.split("_")[0][:4]
            resp = requests.get(RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id), timeout=timeout)

        if resp.status_code == 200:
            return resp.text
    except requests.exceptions.RequestException:
        pass
    return None


def _extract_uniprot_accession(target_id: str) -> Optional[str]:
    """'AF-P12345-F1-model_v4' -> 'P12345'; falls back to the raw id."""
    parts = target_id.split("-")
    if len(parts) >= 2 and parts[0].upper() == "AF":
        return parts[1]
    return target_id or None


def hits_to_dicts(hits: Sequence[Hit]) -> list[dict]:
    """JSON-serializable representation, e.g. for a FastAPI response."""
    return [dataclasses.asdict(h) for h in hits]
