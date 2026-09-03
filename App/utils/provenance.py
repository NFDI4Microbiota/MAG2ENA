"""Unified provenance events for MAG2ENA metadata and submission workflows."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from typing import Iterable
from uuid import uuid4

import pandas as pd


PROVENANCE_COLUMNS = (
    "event_id",
    "timestamp",
    "stage",
    "action",
    "entity_type",
    "sample_name",
    "field",
    "source",
    "before",
    "after",
    "method",
    "method_version",
    "decision",
    "details",
)

_SECRET_KEY_RE = re.compile(
    r"(?:password|passwd|credential|secret|access[_ -]?token|api[_ -]?key)",
    re.IGNORECASE,
)
_EMPTY_TEXT = {"", "<na>", "nan", "none"}


def utc_now() -> str:
    """Return a timezone-aware ISO 8601 timestamp."""
    return datetime.now(timezone.utc).isoformat()


def _text(value) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    text = str(value)
    return "" if text.strip().casefold() in _EMPTY_TEXT else text


def _safe_details(value, *, key: str = ""):
    """Return JSON-safe details while redacting credential-like keys."""
    if _SECRET_KEY_RE.search(str(key)):
        return "[REDACTED]"
    if isinstance(value, dict):
        return {
            str(child_key): _safe_details(child_value, key=str(child_key))
            for child_key, child_value in value.items()
        }
    if isinstance(value, (list, tuple, set)):
        return [_safe_details(item) for item in value]
    if isinstance(value, Path):
        return value.name
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def make_event(
    *,
    stage: str,
    action: str,
    entity_type: str,
    sample_name: str = "",
    field: str = "",
    source: str = "",
    before="",
    after="",
    method: str = "",
    method_version: str = "",
    decision: str = "",
    details: dict | None = None,
    timestamp: str | None = None,
    event_id: str | None = None,
) -> dict:
    """Create a serializable, credential-safe provenance event."""
    if any(_SECRET_KEY_RE.search(str(value)) for value in (field, source, action)):
        before = after = "[REDACTED]"
    return {
        "event_id": event_id or str(uuid4()),
        "timestamp": timestamp or utc_now(),
        "stage": _text(stage),
        "action": _text(action),
        "entity_type": _text(entity_type),
        "sample_name": _text(sample_name),
        "field": _text(field),
        "source": _text(source),
        "before": _text(before),
        "after": _text(after),
        "method": _text(method),
        "method_version": _text(method_version),
        "decision": _text(decision),
        "details": _safe_details(details or {}),
    }


def dataframe_change_events(
    before: pd.DataFrame,
    after: pd.DataFrame,
    *,
    stage: str,
    action: str,
    source: str,
    method: str = "",
    method_version: str = "",
    decision: str = "",
    details: dict | None = None,
    key_column: str = "sample_name",
) -> list[dict]:
    """Return one provenance event per changed dataframe cell."""
    before = before.copy()
    after = after.copy()
    indices = list(dict.fromkeys([*before.index.tolist(), *after.index.tolist()]))
    columns = list(dict.fromkeys([*before.columns.tolist(), *after.columns.tolist()]))
    events = []

    for index in indices:
        sample_name = ""
        if index in after.index and key_column in after.columns:
            sample_name = _text(after.at[index, key_column])
        if not sample_name and index in before.index and key_column in before.columns:
            sample_name = _text(before.at[index, key_column])

        for column in columns:
            old_value = (
                before.at[index, column]
                if index in before.index and column in before.columns
                else ""
            )
            new_value = (
                after.at[index, column]
                if index in after.index and column in after.columns
                else ""
            )
            if _text(old_value) == _text(new_value):
                continue
            event_details = {"row_index": str(index)}
            event_details.update(details or {})
            events.append(
                make_event(
                    stage=stage,
                    action=action,
                    entity_type="metadata_cell",
                    sample_name=sample_name,
                    field=str(column),
                    source=source,
                    before=old_value,
                    after=new_value,
                    method=method,
                    method_version=method_version,
                    decision=decision,
                    details=event_details,
                )
            )
    return events


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def sha256_path(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def dataframe_sha256(df: pd.DataFrame) -> str:
    content = df.fillna("").astype(str).to_csv(index=False, lineterminator="\n")
    return sha256_bytes(content.encode("utf-8"))


def events_to_dataframe(events: Iterable[dict]) -> pd.DataFrame:
    """Flatten provenance events into a stable, human-readable table."""
    rows = []
    for raw_event in events:
        event = make_event(
            stage=raw_event.get("stage", ""),
            action=raw_event.get("action", ""),
            entity_type=raw_event.get("entity_type", ""),
            sample_name=raw_event.get("sample_name", ""),
            field=raw_event.get("field", ""),
            source=raw_event.get("source", ""),
            before=raw_event.get("before", ""),
            after=raw_event.get("after", ""),
            method=raw_event.get("method", ""),
            method_version=raw_event.get("method_version", ""),
            decision=raw_event.get("decision", ""),
            details=raw_event.get("details", {}),
            timestamp=raw_event.get("timestamp") or utc_now(),
            event_id=raw_event.get("event_id") or str(uuid4()),
        )
        event["details"] = json.dumps(
            event["details"], ensure_ascii=False, sort_keys=True
        )
        rows.append(event)
    return pd.DataFrame(rows, columns=PROVENANCE_COLUMNS)


def events_to_tsv(events: Iterable[dict]) -> str:
    return events_to_dataframe(events).to_csv(sep="\t", index=False, lineterminator="\n")


def events_to_json(events: Iterable[dict]) -> str:
    safe_records = events_to_dataframe(events).to_dict(orient="records")
    for record in safe_records:
        record["details"] = json.loads(record["details"] or "{}")
    return json.dumps(safe_records, ensure_ascii=False, indent=2)


def write_provenance_files(events: Iterable[dict], directory: str | Path) -> tuple[Path, Path]:
    """Write the unified human- and machine-readable ledgers."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    events = list(events)
    tsv_path = directory / "provenance.tsv"
    json_path = directory / "provenance.json"
    tsv_path.write_text(events_to_tsv(events), encoding="utf-8")
    json_path.write_text(events_to_json(events), encoding="utf-8")
    return tsv_path, json_path
