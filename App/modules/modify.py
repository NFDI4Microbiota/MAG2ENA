"""Streamlit workflow for metadata-only updates of existing ENA MAG samples."""

from __future__ import annotations

import json
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from pathlib import Path

import pandas as pd
import streamlit as st
from modules.submit import (
    CHECKLIST_XML,
    EXAMPLES_DIR,
    MAG2ENA_VERSION,
    _ena_taxonomy_search,
    _merge_into_metadata,
    _parse_checkm_file,
    _parse_gtdbtk_file,
    build_column_config,
    job_submitted_dialog,
    load_fields_from_xml,
    validate_dataframe,
)
from rq import get_current_job
from utils.db import TaskStatus
from utils.ena_metadata import attrs_to_mag_metadata, fetch_ena_sample_metadata
from utils.ena_sample_update import (
    BATCH_SIZE,
    EnaSampleUpdateError,
    build_modify_xml,
    build_rollback_xml,
    build_sample_set_xml,
    discover_mag_samples_by_study,
    document_to_row,
    fetch_sample_document,
    receipt_success,
    rows_diff,
    sample_candidates_from_accessions,
    submit_modify_xml,
)
from utils.encryption import encrypt_job_folder
from utils.environment_resolver import (
    OlsEnvoClient,
    apply_environment_proposals,
    proposals_to_review_dataframe,
    resolve_environment_batch,
)
from utils.geography import render_geographic_overview
from utils.provenance import (
    dataframe_change_events,
    dataframe_sha256,
    events_to_dataframe,
    events_to_json,
    events_to_tsv,
    make_event,
    sha256_bytes,
    sha256_path,
    write_provenance_files,
)
from utils.tasks import enqueue_task, manager, scrub_job_credentials

_HERE = Path(__file__).parent
_APP_DIR = _HERE.parent
SUBMISSIONS_DIR = _APP_DIR / "jobs"
_MAX_FETCH_WORKERS = 4


# Background job helpers
def clear_session_state() -> None:
    """Discard all Update-page discovery, credentials, editor, and review state."""
    for key in list(st.session_state):
        if key.startswith("modify_"):
            st.session_state.pop(key, None)


def _safe_error(error: Exception, *credentials: str) -> str:
    text = str(error)
    for credential in credentials:
        if credential:
            text = text.replace(str(credential), "[REDACTED]")
    return text


def _append_job_event(ledger: list[dict], directory: Path, **kwargs) -> dict:
    event = make_event(**kwargs)
    ledger.append(event)
    write_provenance_files(ledger, directory)
    return event


def _fetch_documents_parallel(
    candidates: list[dict],
    *,
    portal: str,
    username: str,
    password: str,
    progress=None,
) -> tuple[list[dict], list[dict]]:
    """Fetch at most four authenticated sample XML records concurrently."""
    documents = []
    failures = []
    with ThreadPoolExecutor(max_workers=_MAX_FETCH_WORKERS) as executor:
        future_to_candidate = {
            executor.submit(
                fetch_sample_document,
                candidate["sample_accession"],
                portal=portal,
                username=username,
                password=password,
            ): candidate
            for candidate in candidates
        }
        total = len(future_to_candidate)
        for completed, future in enumerate(as_completed(future_to_candidate), start=1):
            candidate = future_to_candidate[future]
            try:
                document = future.result()
                document["discovery"] = candidate
                documents.append(document)
            except Exception as exc:
                failures.append(
                    {
                        "sample_accession": candidate["sample_accession"],
                        "error": _safe_error(exc, username, password),
                    }
                )
            if progress:
                progress(completed / max(total, 1))
    documents.sort(key=lambda item: item.get("accession", ""))
    return documents, failures


def _modify_task_impl(
    documents,
    before_rows,
    updated_rows,
    editable_fields,
    submission,
    validated_modify_sha256,
    validation_receipt="",
    provenance_events=None,
):
    job = get_current_job()
    job_id = job.id if job else datetime.now(UTC).strftime("%Y%m%d_%H%M%S")
    scrub_job_credentials(job)
    manager.store_start(job_id, TaskStatus.RUNNING)

    output_dir = SUBMISSIONS_DIR / job_id
    output_dir.mkdir(parents=True, exist_ok=True)
    ledger = list(provenance_events or [])
    portal = submission["portal"]
    username = submission["ena_user"]
    password = submission["ena_password"]

    _append_job_event(
        ledger,
        output_dir,
        stage="submission",
        action="metadata_update_started",
        entity_type="sample_batch",
        source=portal,
        after="running",
        method="MAG2ENA metadata-only RQ worker",
        method_version=MAG2ENA_VERSION,
        details={"job_id": job_id, "sample_count": len(updated_rows)},
    )

    if validation_receipt:
        (output_dir / "validation_receipt.xml").write_text(
            validation_receipt, encoding="utf-8"
        )

    # Refetch immediately before MODIFY.  A changed hash means another actor
    # edited ENA after the user reviewed the diff, so the batch is rejected.
    candidates = [
        {"sample_accession": document["accession"]} for document in documents
    ]
    current_documents, fetch_failures = _fetch_documents_parallel(
        candidates,
        portal=portal,
        username=username,
        password=password,
    )
    if fetch_failures:
        raise EnaSampleUpdateError(
            "Could not recheck all samples before update: "
            + "; ".join(
                f"{item['sample_accession']}: {item['error']}" for item in fetch_failures
            )
        )
    original_hashes = {
        document["accession"]: document["canonical_hash"] for document in documents
    }
    changed_remotely = [
        document["accession"]
        for document in current_documents
        if original_hashes.get(document["accession"]) != document["canonical_hash"]
    ]
    if changed_remotely:
        raise EnaSampleUpdateError(
            "ENA metadata changed after review for: " + ", ".join(changed_remotely)
        )
    _append_job_event(
        ledger,
        output_dir,
        stage="validation",
        action="optimistic_concurrency_check",
        entity_type="sample_batch",
        source="authenticated ENA sample XML",
        after="passed",
        method="canonical XML SHA-256 comparison",
        method_version=MAG2ENA_VERSION,
        decision="continue",
        details={"sample_count": len(current_documents)},
    )

    field_defs = load_fields_from_xml(CHECKLIST_XML, include_assembly_fields=False)
    worker_errors = validate_dataframe(
        pd.DataFrame(updated_rows, dtype="string"), field_defs
    )
    if not worker_errors.empty:
        raise EnaSampleUpdateError(
            "Metadata no longer passes worker-side checklist/taxonomy validation: "
            + worker_errors.to_json(orient="records")
        )
    changed_by_sample: dict[str, set[str]] = {}
    for change in rows_diff(before_rows, updated_rows):
        changed_by_sample.setdefault(change["sample_accession"], set()).add(
            change["field"]
        )
    for document in documents:
        conflicting = set(document.get("conflicting_duplicate_tags", {})) & (
            changed_by_sample.get(document["accession"], set())
        )
        if conflicting:
            raise EnaSampleUpdateError(
                f"Cannot edit conflicting repeated attributes for {document['accession']}: "
                + ", ".join(sorted(conflicting))
            )

    before_xml = build_sample_set_xml(documents)
    modify_xml = build_modify_xml(
        documents, updated_rows, editable_fields, validate_only=False
    )
    rollback_xml = build_rollback_xml(documents)
    modify_sha256 = sha256_bytes(modify_xml)
    if modify_sha256 != validated_modify_sha256:
        raise EnaSampleUpdateError(
            "The queued metadata differs from the exact XML accepted by ENA dry-run validation"
        )
    (output_dir / "before.xml").write_bytes(before_xml)
    (output_dir / "modify.xml").write_bytes(modify_xml)
    (output_dir / "rollback.xml").write_bytes(rollback_xml)

    diff = pd.DataFrame(rows_diff(before_rows, updated_rows))
    diff.to_csv(output_dir / "metadata_diff.tsv", sep="\t", index=False)
    _append_job_event(
        ledger,
        output_dir,
        stage="submission",
        action="modify_xml_generated",
        entity_type="file",
        source="authoritative ENA XML and curator-approved diff",
        after=modify_sha256,
        method="preserving ENA SAMPLE XML merge",
        method_version=MAG2ENA_VERSION,
        details={
            "sample_count": len(updated_rows),
            "changed_cell_count": len(diff),
            "modify_xml_sha256": modify_sha256,
            "before_xml_sha256": sha256_bytes(before_xml),
            "rollback_xml_sha256": sha256_bytes(rollback_xml),
        },
    )

    receipt = submit_modify_xml(
        modify_xml,
        portal=portal,
        username=username,
        password=password,
        asynchronous=True,
    )
    (output_dir / "webin_log.xml").write_text(receipt, encoding="utf-8")
    success, messages = receipt_success(receipt)
    _append_job_event(
        ledger,
        output_dir,
        stage="submission",
        action="ena_metadata_update_receipt_received",
        entity_type="receipt",
        source=portal,
        after="success" if success else "failure",
        method="ENA Webin REST v2 MODIFY",
        method_version=MAG2ENA_VERSION,
        decision="accepted" if success else "rejected",
        details={
            "job_id": job_id,
            "sample_count": len(updated_rows),
            "messages": messages,
            "receipt_sha256": sha256_bytes(receipt.encode("utf-8")),
        },
    )
    if not success:
        error_text = "ENA rejected the metadata update. " + " | ".join(messages)
        if any("No new BioSample was created" in message for message in messages):
            error_text += (
                " The authenticated Webin account could not update the existing "
                "BioSample through the selected portal; verify whether the sample belongs "
                "to Production or Testing and that it appears in the account's Samples Report."
            )
        (output_dir / "error.txt").write_text(error_text, encoding="utf-8")
        raise EnaSampleUpdateError(error_text)

    result = {
        "job_type": "metadata_update",
        "submitted": len(updated_rows),
        "modified": len(updated_rows),
        "errors": 0,
        "changed_cells": len(diff),
        "portal": portal,
    }
    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    _append_job_event(
        ledger,
        output_dir,
        stage="submission",
        action="metadata_update_completed",
        entity_type="sample_batch",
        source=portal,
        after="completed",
        method="MAG2ENA metadata-only RQ worker",
        method_version=MAG2ENA_VERSION,
        decision="completed",
        details={
            "job_id": job_id,
            "modified": len(updated_rows),
            "changed_cells": len(diff),
            "result_sha256": sha256_path(result_path),
        },
    )
    return result


def modify_task(
    documents,
    before_rows,
    updated_rows,
    editable_fields,
    submission,
    validated_modify_sha256,
    validation_receipt="",
    provenance_events=None,
    email=None,
    encryption_password=None,
):
    """RQ entry point for metadata-only ENA updates."""
    job = get_current_job()
    job_id = job.id if job else datetime.now(UTC).strftime("%Y%m%d_%H%M%S")
    scrub_job_credentials(job)
    try:
        return _modify_task_impl(
            documents,
            before_rows,
            updated_rows,
            editable_fields,
            submission,
            validated_modify_sha256,
            validation_receipt=validation_receipt,
            provenance_events=provenance_events,
        )
    except Exception as exc:
        output_dir = SUBMISSIONS_DIR / job_id
        output_dir.mkdir(parents=True, exist_ok=True)
        safe = _safe_error(
            exc, submission.get("ena_user", ""), submission.get("ena_password", "")
        )
        (output_dir / "error.txt").write_text(safe, encoding="utf-8")
        ledger = list(provenance_events or [])
        json_path = output_dir / "provenance.json"
        try:
            if json_path.exists():
                loaded = json.loads(json_path.read_text(encoding="utf-8"))
                if isinstance(loaded, list):
                    ledger = loaded
        except Exception:
            pass
        try:
            _append_job_event(
                ledger,
                output_dir,
                stage="submission",
                action="metadata_update_failed",
                entity_type="sample_batch",
                source=submission.get("portal", ""),
                after="failed",
                method="MAG2ENA metadata-only RQ worker",
                method_version=MAG2ENA_VERSION,
                decision="failed",
                details={"job_id": job_id, "error_type": type(exc).__name__, "error": safe},
            )
        except Exception:
            pass
        raise
    finally:
        if encryption_password:
            output_dir = SUBMISSIONS_DIR / job_id
            if output_dir.is_dir() and any(output_dir.iterdir()):
                encrypt_job_folder(
                    output_dir,
                    encryption_password,
                    allowed_root=SUBMISSIONS_DIR,
                )


# Page state and provenance
def _events() -> list[dict]:
    if "modify_provenance_events" not in st.session_state:
        st.session_state.modify_provenance_events = [
            make_event(
                stage="workflow",
                action="metadata_update_workflow_started",
                entity_type="workflow",
                source="Update page",
                method="MAG2ENA web workflow",
                method_version=MAG2ENA_VERSION,
            )
        ]
    return st.session_state.modify_provenance_events


def _append_events(events) -> None:
    ledger = _events()
    ledger.extend(list(events))
    st.session_state.modify_provenance_events = ledger


def _invalidate_validation() -> None:
    for key in (
        "modify_validated_hash",
        "modify_validation_receipt",
        "modify_validation_messages",
    ):
        st.session_state.pop(key, None)


def _clear_loaded_batch() -> None:
    for key in list(st.session_state.keys()):
        if key.startswith("modify_") and key not in {
            "modify_candidates",
            "modify_batch_number",
            "modify_portal",
            "modify_user",
            "modify_password",
            "modify_source_mode",
            "modify_study_accession",
            "modify_sample_accessions",
        }:
            st.session_state.pop(key, None)


def _reset_modify_workflow() -> None:
    for key in list(st.session_state.keys()):
        if key.startswith("modify_"):
            st.session_state.pop(key, None)


def _record_dataframe_changes(before: pd.DataFrame, after: pd.DataFrame, **details) -> None:
    events = dataframe_change_events(
        before,
        after,
        stage="metadata",
        action=details.pop("action", "metadata_edited"),
        source=details.pop("source", "Update metadata editor"),
        method=details.pop("method", "curator edit"),
        method_version=MAG2ENA_VERSION,
        decision=details.pop("decision", "user edited"),
        details=details,
        key_column="sample_name",
    )
    _append_events(events)


def _batch_candidates(candidates: list[dict], number: int) -> list[dict]:
    start = (number - 1) * BATCH_SIZE
    return candidates[start : start + BATCH_SIZE]


def _candidate_preview(candidates: list[dict]) -> pd.DataFrame:
    rows = []
    for item in candidates:
        rows.append(
            {
                "sample accession": item.get("sample_accession", ""),
                "assembly accession(s)": ", ".join(item.get("assembly_accessions", [])),
                "assembly name(s)": ", ".join(item.get("assembly_names", [])),
                "scientific name": item.get("scientific_name", ""),
            }
        )
    return pd.DataFrame(rows)


# Load a batch and make the same assistance used on Home available here
def _load_selected_batch(
    batch: list[dict], portal: str, username: str, password: str, field_defs: dict
) -> None:
    progress = st.progress(0, text="Retrieving authoritative ENA sample XML…")
    documents, failures = _fetch_documents_parallel(
        batch,
        portal=portal,
        username=username,
        password=password,
        progress=lambda value: progress.progress(
            value, text="Retrieving authoritative ENA sample XML…"
        ),
    )
    progress.empty()
    if failures:
        st.warning(f"Could not retrieve {len(failures)} sample(s).")
        st.dataframe(pd.DataFrame(failures), hide_index=True, width="stretch")
    if not documents:
        st.error("No editable sample records were retrieved.")
        return

    direct_non_mimag = []
    accepted_documents = []
    for document in documents:
        discovery = document.get("discovery", {})
        discovered_from_mag_assembly = bool(discovery.get("assembly_accessions"))
        checklist_values = {
            str(item.get("value", "") or "").strip()
            for item in document.get("attributes", {}).get("ENA-CHECKLIST", [])
        }
        if discovered_from_mag_assembly or "ERC000047" in checklist_values:
            accepted_documents.append(document)
        else:
            direct_non_mimag.append(document["accession"])
    documents = accepted_documents
    if direct_non_mimag:
        st.warning(
            "Skipped sample accession(s) that were not discovered through a MAG assembly "
            "and do not declare the ERC000047 MIMAG checklist: "
            + ", ".join(direct_non_mimag)
        )
    if not documents:
        st.error("The selected records do not contain editable ERC000047 MAG samples.")
        return

    # A canonical accession may be reached through both a BioSample and ERS ID.
    unique_documents = {}
    for document in documents:
        unique_documents[document["accession"]] = document
    documents = list(unique_documents.values())
    field_names = ["sample_accession", *field_defs.keys()]
    rows = [document_to_row(document, field_names) for document in documents]
    metadata = pd.DataFrame(rows, columns=field_names, dtype="string")
    st.session_state.modify_documents = documents
    st.session_state.modify_before_df = metadata.copy()
    st.session_state.modify_metadata_df = metadata.copy()
    st.session_state.modify_editor_key = st.session_state.get("modify_editor_key", 0) + 1
    ledger = list(_events())
    ledger.append(
        make_event(
            stage="workflow",
            action="metadata_update_batch_loaded",
            entity_type="workflow",
            source="Update page",
            method="MAG2ENA web workflow",
            method_version=MAG2ENA_VERSION,
            details={"portal": portal, "batch_size": len(documents)},
        )
    )
    st.session_state.modify_provenance_events = ledger
    _append_events(
        make_event(
            stage="metadata",
            action="ena_sample_loaded",
            entity_type="sample",
            sample_name=document.get("alias", ""),
            source=document.get("source_url", "authenticated ENA sample XML"),
            after=document["canonical_hash"],
            method="authenticated ENA sample retrieval",
            method_version=MAG2ENA_VERSION,
            details={
                "sample_accession": document["accession"],
                "xml_sha256": document["canonical_hash"],
            },
        )
        for document in documents
    )
    _invalidate_validation()
    st.rerun()


def _parent_metadata_for_dataframe(df: pd.DataFrame) -> dict[str, dict]:
    accessions = []
    for value in df.get("sample derived from", pd.Series(dtype="string")):
        text = str(value or "").strip()
        if text and text not in {"<NA>", "nan", "None"} and text not in accessions:
            accessions.append(text)

    results = {}
    with ThreadPoolExecutor(max_workers=_MAX_FETCH_WORKERS) as executor:
        future_to_accession = {
            executor.submit(fetch_ena_sample_metadata, accession): accession
            for accession in accessions
        }
        for future in as_completed(future_to_accession):
            accession = future_to_accession[future]
            try:
                record = future.result()
            except Exception as exc:
                record = {"_error": str(exc)}
            if not record.get("_error"):
                results[accession] = record
    return results


def _clear_modify_envo_state(*, keep_undo: bool = True) -> None:
    for key in (
        "modify_envo_proposals",
        "modify_envo_review",
        "modify_envo_fingerprint",
        "modify_envo_review_key",
        "modify_envo_scope",
    ):
        st.session_state.pop(key, None)
    if not keep_undo:
        st.session_state.pop("modify_envo_undo_df", None)


def _store_assisted_update(
    current: pd.DataFrame,
    updated: pd.DataFrame,
    *,
    action: str,
    source: str,
    method: str,
    decision: str = "user applied",
    details: dict | None = None,
    clear_envo: bool = True,
) -> None:
    _record_dataframe_changes(
        current,
        updated,
        action=action,
        source=source,
        method=method,
        decision=decision,
        **(details or {}),
    )
    st.session_state.modify_metadata_df = updated.astype("string")
    st.session_state.modify_editor_key = st.session_state.get("modify_editor_key", 0) + 1
    if clear_envo:
        _clear_modify_envo_state(keep_undo=False)
    _invalidate_validation()


def _render_parent_autofill() -> None:
    current = st.session_state.modify_metadata_df
    accessions = sorted(
        {
            str(value).strip()
            for value in current.get("sample derived from", pd.Series(dtype="string")).dropna()
            if str(value).strip() not in {"", "<NA>", "None", "nan"}
        }
    )
    st.caption(
        "Reuse collection, location, isolation-source, metagenomic-source, and environmental "
        "metadata from `sample derived from`. Existing values are preserved."
    )
    if accessions:
        st.info(
            f"Found {len(accessions)} unique parent accession(s): "
            + ", ".join(f"`{value}`" for value in accessions[:6])
            + ("…" if len(accessions) > 6 else "")
        )
    else:
        st.warning("No accessions are present in `sample derived from`.")

    preview_col, button_col = st.columns([3, 1])
    with preview_col:
        preview_accession = st.text_input(
            "Preview a parent accession",
            placeholder="ERS…, SRS…, or DRS…",
            key="modify_parent_preview_accession",
        )
    with button_col:
        st.write("")
        if st.button("Preview", key="modify_parent_preview_button", width="stretch"):
            if preview_accession.strip():
                with st.spinner("Retrieving parent metadata from ENA…"):
                    raw = fetch_ena_sample_metadata(preview_accession.strip())
                if raw.get("_error"):
                    st.error(raw["_error"])
                    st.session_state.pop("modify_parent_preview", None)
                else:
                    st.session_state.modify_parent_preview = {
                        "mapped": attrs_to_mag_metadata(raw),
                        "source": raw.get("_metadata_source", "ENA"),
                        "assemblies": raw.get("_fallback_assembly_accessions", []),
                        "conflicts": raw.get("_fallback_conflicts", {}),
                    }
            else:
                st.warning("Enter an accession first.")

    preview = st.session_state.get("modify_parent_preview")
    if preview:
        st.caption(f"Metadata source: {preview['source']}")
        rows = [
            {"field": field, "value": value}
            for field, value in preview["mapped"].items()
            if value
        ]
        if rows:
            st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")
        if preview["assemblies"]:
            st.caption("Linked assembly fallback: " + ", ".join(preview["assemblies"]))
        if preview["conflicts"]:
            st.warning(
                "Conflicting fallback values were omitted for: "
                + ", ".join(sorted(preview["conflicts"]))
            )

    if st.button(
        "Autofill from all parent accessions",
        key="modify_parent_autofill_button",
        type="primary",
        width="stretch",
        disabled=not accessions,
    ):
        with st.spinner(f"Retrieving {len(accessions)} parent record(s)…"):
            parent_records = _parent_metadata_for_dataframe(current)
        updated = current.copy()
        fallback_accessions = []
        for index, row in updated.iterrows():
            parent = str(row.get("sample derived from", "") or "").strip()
            raw = parent_records.get(parent)
            if not raw:
                continue
            if raw.get("_fallback_assembly_accessions"):
                fallback_accessions.append(parent)
            for field, value in attrs_to_mag_metadata(raw).items():
                if field not in updated or not str(value or "").strip():
                    continue
                existing = str(updated.at[index, field] or "").strip()
                if existing in {"", "<NA>", "nan", "None"}:
                    updated.at[index, field] = str(value)
        _append_events(
            [
                make_event(
                    stage="metadata",
                    action="parent_metadata_retrieved",
                    entity_type="parent_sample",
                    source=raw.get("_metadata_source", "ENA parent sample"),
                    after=accession,
                    method="parent-record retrieval with exact assembly fallback",
                    method_version=MAG2ENA_VERSION,
                    details={
                        "parent_accession": accession,
                        "source_url": raw.get("_metadata_source_url", ""),
                        "fallback_assembly_accessions": raw.get(
                            "_fallback_assembly_accessions", []
                        ),
                        "omitted_conflicts": raw.get("_fallback_conflicts", {}),
                    },
                )
                for accession, raw in parent_records.items()
            ]
        )
        _store_assisted_update(
            current,
            updated,
            action="parent_metadata_autofilled",
            source="ENA parent sample records",
            method="MAG2ENA parent-record autofill",
            details={"fallback_accessions": sorted(set(fallback_accessions))},
        )
        st.session_state.modify_parent_metadata = parent_records
        st.toast(
            f"Filled {len(rows_diff(current.to_dict('records'), updated.to_dict('records')))} cell(s).",
            icon="✅",
        )
        st.rerun()


def _render_taxonomy_assistance() -> None:
    st.caption(
        "Search ENA taxonomy or resolve missing tax IDs in batch. The final update "
        "validation still requires species-rank, non-metagenome taxa."
    )
    query = st.text_input(
        "Organism name",
        placeholder="uncultured Firmicutes bacterium",
        key="modify_taxonomy_query",
    )
    search_col, batch_col = st.columns(2)
    with search_col:
        if st.button("Search taxonomy", key="modify_taxonomy_search", width="stretch"):
            st.session_state.modify_taxonomy_results = _ena_taxonomy_search(query)
    with batch_col:
        if st.button(
            "Batch resolve missing tax IDs",
            key="modify_taxonomy_batch",
            width="stretch",
        ):
            current = st.session_state.modify_metadata_df
            updated = current.copy()
            evidence = {}
            unresolved = []
            with st.spinner("Resolving organism names through ENA taxonomy…"):
                for index, row in updated.iterrows():
                    organism = str(row.get("organism", "") or "").strip()
                    tax_id = str(row.get("tax_id", "") or "").strip()
                    if not organism or tax_id not in {"", "<NA>", "None", "nan"}:
                        continue
                    hits = _ena_taxonomy_search(organism)
                    if not hits:
                        unresolved.append(organism)
                        continue
                    chosen = next(
                        (
                            hit
                            for hit in hits
                            if str(hit.get("scientificName", "")).casefold()
                            == organism.casefold()
                        ),
                        hits[0],
                    )
                    updated.at[index, "tax_id"] = str(chosen.get("taxId", ""))
                    evidence[str(index)] = chosen
            _store_assisted_update(
                current,
                updated,
                action="taxonomy_resolved",
                source="ENA taxonomy API",
                method="batch taxonomy suggestion",
                details={"resolved_candidates": evidence, "unresolved": unresolved},
            )
            if unresolved:
                st.warning("Unresolved: " + ", ".join(sorted(set(unresolved))[:8]))
            st.rerun()

    results = st.session_state.get("modify_taxonomy_results")
    if results is not None:
        if not results:
            st.info("No taxonomy results found.")
            return
        shown = results[:10]
        st.dataframe(
            pd.DataFrame(
                [
                    {
                        "tax_id": item.get("taxId", ""),
                        "scientific name": item.get("scientificName", ""),
                        "rank": item.get("rank", ""),
                    }
                    for item in shown
                ]
            ),
            hide_index=True,
            width="stretch",
        )
        selected = st.selectbox(
            "Taxonomy candidate",
            options=range(len(shown)),
            format_func=lambda index: (
                f"{shown[index].get('scientificName', '')} "
                f"(taxId: {shown[index].get('taxId', '')})"
            ),
            key="modify_taxonomy_candidate",
        )
        sample_names = st.multiselect(
            "Apply to loaded MAGs",
            options=st.session_state.modify_metadata_df["sample_name"].tolist(),
            key="modify_taxonomy_samples",
        )
        if st.button(
            f"Apply candidate to {len(sample_names)} MAG(s)",
            key="modify_taxonomy_apply",
            type="primary",
            width="stretch",
            disabled=not sample_names,
        ):
            current = st.session_state.modify_metadata_df
            updated = current.copy()
            chosen = shown[selected]
            mask = updated["sample_name"].isin(sample_names)
            updated.loc[mask, "organism"] = str(chosen.get("scientificName", ""))
            updated.loc[mask, "tax_id"] = str(chosen.get("taxId", ""))
            _store_assisted_update(
                current,
                updated,
                action="taxonomy_selected",
                source="ENA taxonomy API",
                method="curator-selected taxonomy candidate",
                decision="curator approved",
                details={
                    "query": query,
                    "selected_tax_id": chosen.get("taxId", ""),
                    "selected_scientific_name": chosen.get("scientificName", ""),
                    "rank": chosen.get("rank", ""),
                },
            )
            st.session_state.pop("modify_taxonomy_results", None)
            st.rerun()


def _render_environment_assistance() -> None:
    st.caption(
        "Resolve broad-scale context, local context, and environmental medium with "
        "field-aware ENVO ranking. Nothing is applied without curator acceptance."
    )
    scope = st.radio(
        "Resolution scope",
        ["Fill empty fields", "Review and normalize all fields"],
        horizontal=True,
        key="modify_resolver_scope",
    )
    if (
        st.session_state.get("modify_envo_proposals")
        and st.session_state.get("modify_envo_scope") != scope
    ):
        _clear_modify_envo_state(keep_undo=True)
    generate_col, undo_col = st.columns([3, 1])
    with generate_col:
        generate = st.button(
            "Generate batch proposals",
            key="modify_envo_generate",
            type="primary",
            width="stretch",
        )
    with undo_col:
        if st.button(
            "Undo last apply",
            key="modify_envo_undo",
            width="stretch",
            disabled="modify_envo_undo_df" not in st.session_state,
        ):
            current = st.session_state.modify_metadata_df
            restored = st.session_state.pop("modify_envo_undo_df")
            _store_assisted_update(
                current,
                restored,
                action="environment_resolution_undone",
                source="environmental resolver",
                method="undo last curator-approved apply",
                decision="curator reversed",
            )
            st.rerun()

    if generate:
        current = st.session_state.modify_metadata_df
        parent_records = _parent_metadata_for_dataframe(current)
        st.session_state.modify_parent_metadata = parent_records
        with st.spinner("Querying OLS4 and ranking field-aware ENVO candidates…"):
            proposals = resolve_environment_batch(
                current,
                OlsEnvoClient(),
                empty_only=scope == "Fill empty fields",
                parent_metadata_by_accession=parent_records,
            )
        st.session_state.modify_envo_proposals = proposals
        st.session_state.modify_envo_review = proposals_to_review_dataframe(proposals)
        st.session_state.modify_envo_fingerprint = dataframe_sha256(current)
        st.session_state.modify_envo_scope = scope
        st.session_state.modify_envo_review_key = (
            st.session_state.get("modify_envo_review_key", 0) + 1
        )
        _append_events(
            make_event(
                stage="metadata",
                action="environment_term_proposed",
                entity_type="ontology_proposal",
                field=proposal.get("field", ""),
                source="OLS4 / ENVO",
                after=proposal.get("selected_value", ""),
                method="field-aware environmental resolver",
                method_version=proposal.get("ontology_version", MAG2ENA_VERSION),
                decision="pending curator review",
                details={
                    "proposal_id": proposal.get("proposal_id", ""),
                    "sample_names": proposal.get("sample_names", []),
                    "source_evidence": proposal.get("evidence", ""),
                    "confidence": proposal.get("confidence", ""),
                    "score": proposal.get("score", 0),
                    "reason": proposal.get("reason", ""),
                    "alternatives": proposal.get("candidates", []),
                },
            )
            for proposal in proposals
        )
        st.rerun()

    proposals = st.session_state.get("modify_envo_proposals")
    if not proposals:
        return
    review = st.data_editor(
        st.session_state.modify_envo_review,
        hide_index=True,
        width="stretch",
        height=min(520, 90 + 35 * len(proposals)),
        disabled=[
            "proposal_id",
            "rows",
            "field",
            "source evidence",
            "proposed value",
            "confidence",
            "score",
            "reason",
        ],
        column_config={"accept": st.column_config.CheckboxColumn("Accept", default=False)},
        key=f"modify_envo_review_editor_{st.session_state.get('modify_envo_review_key', 0)}",
    )
    st.session_state.modify_envo_review = review
    select_col, clear_col = st.columns(2)
    with select_col:
        if st.button(
            "Select verified and high-confidence",
            key="modify_envo_select_safe",
            width="stretch",
        ):
            review.loc[
                review["confidence"].isin(["Verified", "High"])
                & review["proposed value"].astype(str).str.len().gt(0),
                "accept",
            ] = True
            st.session_state.modify_envo_review = review
            st.session_state.modify_envo_review_key += 1
            st.rerun()
    with clear_col:
        if st.button("Clear selections", key="modify_envo_clear", width="stretch"):
            review["accept"] = False
            st.session_state.modify_envo_review = review
            st.session_state.modify_envo_review_key += 1
            st.rerun()

    with_candidates = [proposal for proposal in proposals if proposal.get("candidates")]
    if with_candidates:
        with st.expander("Inspect alternatives or change a proposal"):
            lookup = {proposal["proposal_id"]: proposal for proposal in with_candidates}
            proposal_id = st.selectbox(
                "Proposal",
                options=list(lookup),
                format_func=lambda value: (
                    f"{value} · {lookup[value]['field']} · {lookup[value]['row_count']} row(s)"
                ),
                key="modify_envo_alternative_proposal",
            )
            proposal = lookup[proposal_id]
            candidate_index = st.selectbox(
                "Candidate",
                options=range(len(proposal["candidates"])),
                format_func=lambda index: (
                    f"{proposal['candidates'][index]['value']} · "
                    f"{proposal['candidates'][index]['score']:.3f}"
                ),
                key=f"modify_envo_alternative_{proposal_id}",
            )
            candidate = proposal["candidates"][candidate_index]
            if candidate.get("description"):
                st.caption(candidate["description"][:500])
            if st.button(
                "Use this candidate",
                key="modify_envo_use_alternative",
                width="stretch",
            ):
                previous = proposal.get("selected_value", "")
                proposal.update(
                    {
                        "selected_value": candidate["value"],
                        "selected_label": candidate["label"],
                        "selected_id": candidate["obo_id"],
                        "score": candidate["score"],
                        "confidence": "User selected",
                        "reason": candidate["reason"],
                    }
                )
                _append_events(
                    [
                        make_event(
                            stage="metadata",
                            action="environment_alternative_selected",
                            entity_type="ontology_proposal",
                            field=proposal["field"],
                            source="OLS4 / ENVO",
                            before=previous,
                            after=candidate["value"],
                            method="curator-selected resolver alternative",
                            method_version=proposal.get(
                                "ontology_version", MAG2ENA_VERSION
                            ),
                            decision="curator selected candidate",
                            details={"proposal_id": proposal_id},
                        )
                    ]
                )
                previous_decisions = dict(zip(review["proposal_id"], review["accept"]))
                replacement = proposals_to_review_dataframe(proposals)
                replacement["accept"] = replacement["proposal_id"].map(
                    previous_decisions
                ).fillna(False)
                st.session_state.modify_envo_proposals = proposals
                st.session_state.modify_envo_review = replacement
                st.session_state.modify_envo_review_key += 1
                st.rerun()

    accepted = review.loc[review["accept"], "proposal_id"].tolist()
    if st.button(
        f"Apply {len(accepted)} accepted proposal(s)",
        key="modify_envo_apply",
        disabled=not accepted,
        type="primary",
        width="stretch",
    ):
        current = st.session_state.modify_metadata_df
        if dataframe_sha256(current) != st.session_state.modify_envo_fingerprint:
            st.error("The metadata changed. Generate fresh resolver proposals.")
            return
        updated, audit = apply_environment_proposals(
            current,
            proposals,
            accepted,
            overwrite=scope == "Review and normalize all fields",
        )
        _append_events(
            make_event(
                stage="metadata",
                action="environment_proposal_reviewed",
                entity_type="ontology_proposal",
                field=proposal.get("field", ""),
                source="OLS4 / ENVO",
                after=proposal.get("selected_value", ""),
                method="field-aware environmental resolver",
                method_version=proposal.get("ontology_version", MAG2ENA_VERSION),
                decision=(
                    "curator approved"
                    if proposal.get("proposal_id") in accepted
                    else "not selected"
                ),
                details={"proposal_id": proposal.get("proposal_id", "")},
            )
            for proposal in proposals
        )
        st.session_state.modify_envo_undo_df = current.copy()
        _store_assisted_update(
            current,
            updated,
            action="environment_term_applied",
            source="OLS4 / ENVO",
            method="field-aware environmental resolver",
            decision="curator approved",
            details={"audit": audit},
            clear_envo=False,
        )
        _clear_modify_envo_state(keep_undo=True)
        st.rerun()


def _render_fill_column(field_defs: dict) -> None:
    st.caption(
        "Apply a value across a metadata field. Immutable sample accessions and aliases "
        "are excluded; keep `Fill empty cells only` selected unless replacement is intended."
    )
    columns = [
        column
        for column in st.session_state.modify_metadata_df.columns
        if column not in {"sample_accession", "sample_name"}
    ]
    field = st.selectbox("Column", columns, key="modify_fill_column")
    empty_only = st.checkbox(
        "Fill empty cells only", value=True, key="modify_fill_empty_only"
    )
    definition = field_defs.get(field, {})
    if definition.get("type") == "enum":
        value = st.selectbox(
            "Value", definition.get("enum", []), key="modify_fill_enum_value"
        )
    else:
        value = st.text_input("Value", key="modify_fill_text_value")
    if st.button("Apply to column", key="modify_fill_apply", type="primary"):
        if not str(value).strip():
            st.warning("Enter a value first.")
            return
        current = st.session_state.modify_metadata_df
        updated = current.copy()
        if empty_only:
            mask = updated[field].isna() | updated[field].astype(str).str.strip().isin(
                ["", "<NA>", "None", "nan"]
            )
            updated.loc[mask, field] = value
        else:
            updated[field] = value
        _store_assisted_update(
            current,
            updated,
            action="column_filled",
            source="Update fill-column tool",
            method="bulk metadata update",
            details={"field": field, "empty_cells_only": empty_only},
        )
        st.rerun()


def _render_quality_imports() -> None:
    st.caption(
        "Import CheckM/CheckM2 quality metrics or GTDB-Tk organism names by original "
        "sample alias. Existing values are preserved."
    )
    checkm_col, gtdb_col = st.columns(2)
    with checkm_col:
        st.markdown("**CheckM / CheckM2**")
        example_columns = st.columns(2)
        for container, filename, label in (
            (example_columns[0], "checkm_example.tsv", "CheckM example"),
            (example_columns[1], "checkm2_example.tsv", "CheckM2 example"),
        ):
            path = EXAMPLES_DIR / filename
            if path.exists():
                with container:
                    st.download_button(
                        label,
                        data=path.read_bytes(),
                        file_name=filename,
                        mime="text/tab-separated-values",
                        width="stretch",
                    )
        checkm = st.file_uploader(
            "quality_report.tsv or storage.tsv",
            type=["tsv", "txt"],
            key="modify_checkm_upload",
        )
        if checkm:
            parsed = _parse_checkm_file(checkm)
            if parsed is None:
                st.error("Unrecognised CheckM/CheckM2 format.")
            elif st.button("Apply CheckM data", key="modify_checkm_apply"):
                current = st.session_state.modify_metadata_df
                updated = _merge_into_metadata(current, parsed)
                content = checkm.getvalue()
                _append_events(
                    [
                        make_event(
                            stage="metadata",
                            action="quality_file_imported",
                            entity_type="file",
                            source=checkm.name,
                            after=sha256_bytes(content),
                            method="CheckM/CheckM2 parser",
                            method_version=MAG2ENA_VERSION,
                            details={"size_bytes": len(content)},
                        )
                    ]
                )
                _store_assisted_update(
                    current,
                    updated,
                    action="quality_metrics_imported",
                    source=checkm.name,
                    method="CheckM/CheckM2 parser",
                )
                st.rerun()

    with gtdb_col:
        st.markdown("**GTDB-Tk**")
        example = EXAMPLES_DIR / "gtdbtk_example.tsv"
        if example.exists():
            st.download_button(
                "GTDB-Tk example",
                data=example.read_bytes(),
                file_name=example.name,
                mime="text/tab-separated-values",
                width="stretch",
            )
        gtdbtk = st.file_uploader(
            "gtdbtk.bac120.summary.tsv or gtdbtk.ar53.summary.tsv",
            type=["tsv"],
            key="modify_gtdtk_upload",
        )
        if gtdbtk:
            parsed = _parse_gtdbtk_file(gtdbtk)
            if parsed is None:
                st.error("Unrecognised GTDB-Tk summary format.")
            elif st.button("Apply GTDB-Tk data", key="modify_gtdbtk_apply"):
                current = st.session_state.modify_metadata_df
                updated = _merge_into_metadata(current, parsed)
                content = gtdbtk.getvalue()
                _append_events(
                    [
                        make_event(
                            stage="metadata",
                            action="taxonomy_file_imported",
                            entity_type="file",
                            source=gtdbtk.name,
                            after=sha256_bytes(content),
                            method="GTDB-Tk rank-aware parser",
                            method_version=MAG2ENA_VERSION,
                            details={"size_bytes": len(content)},
                        )
                    ]
                )
                _store_assisted_update(
                    current,
                    updated,
                    action="organism_names_imported",
                    source=gtdbtk.name,
                    method="GTDB-Tk rank-aware parser",
                )
                st.rerun()


def _render_assistance(field_defs: dict) -> None:
    with st.expander("Metadata assistance", expanded=False):
        st.caption(
            "The same assistance used during initial submission is available for existing "
            "MAG samples. Every applied value remains editable and is recorded in provenance."
        )
        parent_tab, tax_tab, envo_tab, fill_tab, import_tab = st.tabs(
            [
                "ENA autofill",
                "Taxonomy resolver",
                "Environmental resolver",
                "Fill column",
                "Import quality files",
            ]
        )
        with parent_tab:
            _render_parent_autofill()
        with tax_tab:
            _render_taxonomy_assistance()
        with envo_tab:
            _render_environment_assistance()
        with fill_tab:
            _render_fill_column(field_defs)
        with import_tab:
            _render_quality_imports()


# Review and submission
def _changed_payload() -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    before_rows = st.session_state.modify_before_df.fillna("").to_dict(orient="records")
    updated_rows = st.session_state.modify_metadata_df.fillna("").to_dict(orient="records")
    changes = rows_diff(before_rows, updated_rows)
    changed_accessions = {item["sample_accession"] for item in changes}
    documents = [
        document
        for document in st.session_state.modify_documents
        if document["accession"] in changed_accessions
    ]
    before_rows = [
        row for row in before_rows if row["sample_accession"] in changed_accessions
    ]
    updated_rows = [
        row for row in updated_rows if row["sample_accession"] in changed_accessions
    ]
    return documents, before_rows, updated_rows, changes


def _render_provenance() -> None:
    with st.expander("Provenance and history", expanded=False):
        events = _events()
        history = events_to_dataframe(events)
        st.caption(
            "Append-only history of retrieval, metadata evidence, edits, resolver decisions, "
            "validation, confirmation, and ENA outcomes. Credentials are never recorded."
        )
        if history.empty:
            st.info("No provenance events have been recorded yet.")
            return
        columns = [
            "timestamp",
            "stage",
            "action",
            "sample_name",
            "field",
            "source",
            "before",
            "after",
            "decision",
        ]
        st.dataframe(
            history[columns],
            hide_index=True,
            width="stretch",
            height=min(480, 80 + 35 * len(history)),
        )
        tsv_col, json_col = st.columns(2)
        with tsv_col:
            st.download_button(
                "Download provenance.tsv",
                data=events_to_tsv(events).encode("utf-8"),
                file_name="modify_provenance.tsv",
                mime="text/tab-separated-values",
                width="stretch",
                key="download_modify_provenance_tsv",
            )
        with json_col:
            st.download_button(
                "Download provenance.json",
                data=events_to_json(events).encode("utf-8"),
                file_name="modify_provenance.json",
                mime="application/json",
                width="stretch",
                key="download_modify_provenance_json",
            )


def _current_xml(field_defs: dict, *, validate_only: bool = False) -> bytes:
    documents, _, updated_rows, _ = _changed_payload()
    return build_modify_xml(
        documents,
        updated_rows,
        field_defs.keys(),
        validate_only=validate_only,
    )


def _local_update_errors(field_defs: dict, diff: pd.DataFrame) -> pd.DataFrame:
    changed_accessions = set(diff["sample_accession"]) if not diff.empty else set()
    changed_rows = st.session_state.modify_metadata_df[
        st.session_state.modify_metadata_df["sample_accession"].isin(changed_accessions)
    ].copy()
    errors = validate_dataframe(changed_rows.reset_index(drop=True), field_defs)
    extra = []
    changed_by_sample: dict[str, set[str]] = {}
    if not diff.empty:
        for _, change in diff.iterrows():
            changed_by_sample.setdefault(change["sample_accession"], set()).add(
                change["field"]
            )
    for document in st.session_state.modify_documents:
        changed_fields = changed_by_sample.get(document["accession"], set())
        conflicting = set(document.get("conflicting_duplicate_tags", {})) & changed_fields
        if conflicting:
            extra.append(
                {
                    "row": document["accession"],
                    "field": ", ".join(sorted(conflicting)),
                    "value": "conflicting repeated ENA attributes",
                    "expected": "Resolve conflicting duplicate attributes in ENA before editing them here",
                }
            )
    if extra:
        errors = pd.concat([errors, pd.DataFrame(extra)], ignore_index=True)
    return errors


def runUI():
    st.markdown("## Update submitted MAG metadata")
    st.caption(
        "Update existing ENA MAG sample records without resubmitting assemblies. "
        "MAG2ENA preserves the authoritative sample XML, validates a field-level diff, "
        "and records retrieval, edits, validation, decisions, and ENA outcomes in provenance."
    )
    st.info(
        "Sample accessions, aliases, centre names, assembly files, assembly versions, and "
        "study/sample relationships are not changed by this workflow. Select the portal "
        "where the samples were submitted and use their owning Webin account; public XML "
        "from the other portal or another account may be readable but cannot be modified."
    )

    loaded = bool(st.session_state.get("modify_documents"))
    field_defs = load_fields_from_xml(CHECKLIST_XML, include_assembly_fields=False)

    portal_col, user_col, password_col = st.columns(3)
    with portal_col:
        portal = st.selectbox(
            "ENA portal",
            ["Testing", "Production"],
            key="modify_portal",
            disabled=loaded,
            help=(
                "Choose the portal containing the records. Testing is temporary and "
                "cannot modify the authoritative Production record."
            ),
        )
    with user_col:
        username = st.text_input(
            "ENA Webin username",
            key="modify_user",
            help="Use the Webin account that owns the samples on the selected portal.",
        )
    with password_col:
        password = st.text_input(
            "ENA Webin password",
            type="password",
            key="modify_password",
            help=(
                "Required for authenticated retrieval, validation, and modification; "
                "it is not stored in provenance."
            ),
        )
    has_credentials = bool(username.strip() and password.strip())

    if not loaded:
        mode = st.radio(
            "How should MAG2ENA find the MAG samples?",
            ["Study accession", "Sample accessions"],
            horizontal=True,
            key="modify_source_mode",
            help=(
                "Use a study to discover linked MAG assemblies, or enter one or more "
                "known sample accessions directly."
            ),
        )
        if mode == "Study accession":
            study_accession = st.text_input(
                "Study accession",
                placeholder="PRJEB…, PRJNA…, PRJDB…, ERP…, SRP…, or DRP…",
                key="modify_study_accession",
                help=(
                    "MAG2ENA loads only samples linked to assemblies classified by ENA "
                    "as MAGs."
                ),
            )
            if st.button(
                "Discover MAG samples",
                type="primary",
                width="stretch",
                disabled=not (study_accession.strip() and has_credentials),
                help="Enter a study accession and Webin credentials first.",
            ):
                try:
                    with st.spinner("Finding MAG assemblies and their samples in ENA…"):
                        candidates = discover_mag_samples_by_study(
                            study_accession,
                            portal=portal,
                            username=username,
                            password=password,
                        )
                    st.session_state.modify_candidates = candidates
                    st.session_state.modify_batch_number = 1
                    _events().append(
                        make_event(
                            stage="discovery",
                            action="mag_samples_discovered",
                            entity_type="study",
                            source="ENA Portal API assembly index",
                            after=str(len(candidates)),
                            method="study accession and exact MAG assembly-type query",
                            method_version=MAG2ENA_VERSION,
                            details={"study_accession": study_accession},
                        )
                    )
                except Exception as exc:
                    st.error(_safe_error(exc, username, password))
        else:
            accessions_text = st.text_area(
                "Sample accessions",
                placeholder="ERS… or SAMEA…, separated by lines, spaces, or commas",
                key="modify_sample_accessions",
                help=(
                    "Use ENA sample accessions (ERS/SRS/DRS) or BioSample accessions "
                    "(SAMEA/SAMN/SAMD); duplicates are removed."
                ),
            )
            if st.button(
                "Prepare sample list",
                type="primary",
                width="stretch",
                disabled=not accessions_text.strip(),
                help="Enter at least one sample accession.",
            ):
                try:
                    candidates = sample_candidates_from_accessions(accessions_text)
                    if not candidates:
                        raise ValueError("Enter at least one sample accession")
                    st.session_state.modify_candidates = candidates
                    st.session_state.modify_batch_number = 1
                    _events().append(
                        make_event(
                            stage="discovery",
                            action="sample_accession_list_prepared",
                            entity_type="sample_batch",
                            source="curator-provided sample accessions",
                            after=str(len(candidates)),
                            method="validated accession parsing and deduplication",
                            method_version=MAG2ENA_VERSION,
                            details={"sample_count": len(candidates)},
                        )
                    )
                except ValueError as exc:
                    st.error(str(exc))

        candidates = st.session_state.get("modify_candidates", [])
        if candidates:
            batch_count = math.ceil(len(candidates) / BATCH_SIZE)
            st.success(
                f"Found {len(candidates):,} distinct MAG sample(s), divided into "
                f"{batch_count:,} batch(es) of at most {BATCH_SIZE:,}."
            )
            batch_number = st.number_input(
                "Batch to review",
                min_value=1,
                max_value=batch_count,
                step=1,
                key="modify_batch_number",
            )
            batch = _batch_candidates(candidates, int(batch_number))
            st.dataframe(
                _candidate_preview(batch),
                hide_index=True,
                width="stretch",
                height=min(500, 80 + 35 * len(batch)),
            )
            st.caption(
                f"Showing batch {int(batch_number):,}/{batch_count:,}: "
                f"{len(batch):,} sample(s). Only this batch will be loaded or updated."
            )
            if st.button(
                f"Load batch {int(batch_number):,}",
                type="primary",
                width="stretch",
                disabled=not has_credentials,
                help=(
                    "Loads only this batch after checking that it is available to the "
                    "account on the selected portal."
                ),
            ):
                _load_selected_batch(batch, portal, username, password, field_defs)
        return

    header_col, reset_col = st.columns([4, 1])
    with header_col:
        st.success(
            f"Loaded {len(st.session_state.modify_documents):,} authenticated ENA sample "
            "record(s) in the current batch."
        )
    with reset_col:
        if st.button("Choose another batch", width="stretch"):
            _clear_loaded_batch()
            st.rerun()

    metadata_before_edit = st.session_state.modify_metadata_df.copy()
    column_config = build_column_config(field_defs)
    column_config["sample_accession"] = st.column_config.TextColumn(
        "ENA sample accession", help="Immutable ENA object identifier"
    )
    column_config["sample_name"] = st.column_config.TextColumn(
        "Original alias", help="Immutable original ENA sample alias"
    )
    edited = st.data_editor(
        metadata_before_edit,
        hide_index=True,
        width="stretch",
        height=650,
        disabled=["sample_accession", "sample_name"],
        column_config=column_config,
        key=f"modify_metadata_editor_{st.session_state.get('modify_editor_key', 0)}",
    ).astype("string")
    if not edited.equals(metadata_before_edit):
        _record_dataframe_changes(metadata_before_edit, edited)
        st.session_state.modify_metadata_df = edited
        _invalidate_validation()

    render_geographic_overview(
        st.session_state.modify_metadata_df,
        key="update_geographic_overview",
    )

    _render_assistance(field_defs)
    _render_provenance()

    documents, before_rows, updated_rows, changes = _changed_payload()
    diff = pd.DataFrame(changes)
    st.markdown("### Review and validate")
    if diff.empty:
        st.info("No metadata values have changed.")
        return

    metric_samples, metric_cells = st.columns(2)
    metric_samples.metric("Samples to update", len(updated_rows))
    metric_cells.metric("Changed metadata cells", len(diff))
    st.dataframe(diff, hide_index=True, width="stretch", height=min(500, 80 + 35 * len(diff)))
    st.download_button(
        "Download proposed diff",
        data=diff.to_csv(sep="\t", index=False),
        file_name="metadata_diff.tsv",
        mime="text/tab-separated-values",
        width="stretch",
    )
    if "tax_id" in set(diff["field"]):
        st.warning(
            "A sample taxonomy change is included. ENA advises contacting its helpdesk "
            "when an associated assembly's displayed taxonomy must also be corrected."
        )

    local_errors = _local_update_errors(field_defs, diff)
    if not local_errors.empty:
        st.error(f"Local validation found {len(local_errors)} issue(s).")
        st.dataframe(local_errors, hide_index=True, width="stretch")
        return
    st.success("Local ERC000047 and species-level taxonomy validation passed.")

    actual_xml = _current_xml(field_defs)
    actual_hash = sha256_bytes(actual_xml)
    if st.button(
        "Validate proposed update with ENA",
        type="primary",
        width="stretch",
        disabled=not has_credentials,
        help="Runs ENA MODIFY + VALIDATE without changing the sample records.",
    ):
        if has_credentials:
            try:
                with st.spinner("Running ENA VALIDATE + MODIFY without changing records…"):
                    validation_xml = _current_xml(field_defs, validate_only=True)
                    receipt = submit_modify_xml(
                        validation_xml,
                        portal=portal,
                        username=username,
                        password=password,
                        asynchronous=False,
                    )
                # Webin receipts do not echo HTTP Basic credentials. Preserve
                # the XML exactly; blind string replacement can corrupt tags
                # when a credential happens to be a short/common substring.
                safe_receipt = receipt
                success, messages = receipt_success(safe_receipt)
                st.session_state.modify_validation_receipt = safe_receipt
                st.session_state.modify_validation_messages = messages
                if success:
                    st.session_state.modify_validated_hash = actual_hash
                    _append_events(
                        [
                            make_event(
                                stage="validation",
                                action="ena_modify_dry_run",
                                entity_type="sample_batch",
                                source=portal,
                                after="passed",
                                method="ENA Webin REST v2 MODIFY + VALIDATE",
                                method_version=MAG2ENA_VERSION,
                                decision="eligible for update",
                                details={
                                    "sample_count": len(updated_rows),
                                    "changed_cells": len(diff),
                                    "modify_xml_sha256": actual_hash,
                                    "messages": messages,
                                },
                            )
                        ]
                    )
                    st.success("ENA accepted the dry-run validation. No records were changed.")
                else:
                    st.session_state.pop("modify_validated_hash", None)
                    st.error("ENA rejected the proposed update: " + " | ".join(messages))
            except Exception as exc:
                st.session_state.pop("modify_validated_hash", None)
                st.error(_safe_error(exc, username, password))

    validated = st.session_state.get("modify_validated_hash") == actual_hash
    if st.session_state.get("modify_validated_hash") and not validated:
        st.warning("Metadata changed after ENA validation. Validate the new diff again.")
    if not validated:
        return

    st.success("The exact current diff has passed ENA dry-run validation.")
    if portal == "Production":
        st.warning(
            "Production updates modify the existing public/private ENA sample records. "
            "They do not create new sample accessions and are not an assembly resubmission."
        )
    notification_col, encryption_col = st.columns(2)
    with notification_col:
        notification_email = st.text_input(
            "Notification email (optional)", key="modify_notification_email"
        )
    with encryption_col:
        encryption_password = st.text_input(
            "Results password (optional)",
            type="password",
            key="modify_encryption_password",
            help="Encrypts XML, receipts, diff, rollback file, and provenance at rest.",
        )
    confirmed = st.checkbox(
        f"I approve updating {len(updated_rows):,} existing ENA sample record(s) "
        f"with {len(diff):,} metadata change(s).",
        key="modify_final_confirmation",
        help="Required because Production updates alter existing ENA sample metadata.",
    )
    if st.button(
        "Queue metadata update",
        type="primary",
        width="stretch",
        disabled=not (confirmed and has_credentials and documents and updated_rows),
        help="Available only after validation, explicit approval, and credential entry.",
    ):
        if not documents or not updated_rows or not changes:
            st.error("The metadata update is empty; change at least one value first.")
            return
        _append_events(
            [
                make_event(
                    stage="submission",
                    action="metadata_update_confirmed",
                    entity_type="sample_batch",
                    source="Update page",
                    after="approved",
                    method="explicit curator confirmation",
                    method_version=MAG2ENA_VERSION,
                    decision="user approved",
                    details={
                        "portal": portal,
                        "sample_count": len(updated_rows),
                        "changed_cells": len(diff),
                        "modify_xml_sha256": actual_hash,
                    },
                )
            ]
        )
        submission = {
            "portal": portal,
            "ena_user": username,
            "ena_password": password,
        }
        kwargs = {
            "documents": documents,
            "before_rows": before_rows,
            "updated_rows": updated_rows,
            "editable_fields": list(field_defs),
            "submission": submission,
            "validated_modify_sha256": actual_hash,
            "validation_receipt": st.session_state.get(
                "modify_validation_receipt", ""
            ),
            "provenance_events": list(_events()),
        }
        if notification_email:
            kwargs["email"] = notification_email
        if encryption_password:
            kwargs["encryption_password"] = encryption_password
        job_id = enqueue_task(modify_task, kwargs)
        job_submitted_dialog(job_id, encrypted=bool(encryption_password))


if __name__ == "__main__":
    runUI()
