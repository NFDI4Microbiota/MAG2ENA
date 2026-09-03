"""Reusable metadata preparation helpers for the web and command-line workflows."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pandas as pd

from .ena_metadata import attrs_to_mag_metadata, fetch_ena_sample_metadata


_EMPTY = {"", "<na>", "none", "nan"}


def is_empty(value: object) -> bool:
    """Return whether a metadata cell has no curator-supplied value."""
    return pd.isna(value) or str(value).strip().lower() in _EMPTY


def parse_checkm_file(path_or_buffer: str | Path | object) -> pd.DataFrame:
    """Parse CheckM or CheckM2 TSV output into MAG2ENA quality columns."""
    frame = pd.read_csv(path_or_buffer, sep="\t")
    columns = set(frame.columns)
    if {"Name", "Completeness", "Contamination"}.issubset(columns):
        name_column, software = "Name", "CheckM2"
    elif {"Bin Id", "Completeness", "Contamination"}.issubset(columns):
        name_column, software = "Bin Id", "CheckM"
    else:
        raise ValueError(
            "Unrecognized quality report: expected CheckM 'Bin Id' or CheckM2 "
            "'Name' together with Completeness and Contamination columns."
        )
    return pd.DataFrame(
        {
            "sample_name": frame[name_column].astype(str),
            "completeness score": frame["Completeness"].round(2).astype(str),
            "contamination score": frame["Contamination"].round(2).astype(str),
            "completeness software": software,
        }
    )


def _gtdb_name(classification: str) -> str:
    ranks = {}
    for part in str(classification or "").split(";"):
        if "__" in part:
            rank, name = part.split("__", 1)
            ranks[rank.strip()] = name.strip().replace("_", " ")
    for rank in ("s", "g", "f", "o", "c", "p"):
        name = ranks.get(rank, "")
        if name:
            return name if rank == "s" else f"uncultured {name} bacterium"
    return "uncultured bacterium"


def parse_gtdbtk_file(path_or_buffer: str | Path | object) -> pd.DataFrame:
    """Parse a GTDB-Tk summary into MAG2ENA sample and organism columns."""
    frame = pd.read_csv(path_or_buffer, sep="\t")
    if not {"user_genome", "classification"}.issubset(frame.columns):
        raise ValueError(
            "Unrecognized GTDB-Tk report: expected user_genome and classification columns."
        )
    return pd.DataFrame(
        {
            "sample_name": frame["user_genome"].astype(str),
            "organism": frame["classification"].apply(_gtdb_name),
        }
    )


def merge_empty_cells(base: pd.DataFrame, overlay: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Merge known overlay columns by sample name without overwriting existing cells."""
    result = base.copy()
    changed = 0
    if "sample_name" not in result.columns or "sample_name" not in overlay.columns:
        raise ValueError("Both tables must contain sample_name")
    row_by_name = {
        str(value): index for index, value in result["sample_name"].items()
    }
    for _, overlay_row in overlay.iterrows():
        index = row_by_name.get(str(overlay_row["sample_name"]))
        if index is None:
            continue
        for column in overlay.columns:
            if column == "sample_name" or column not in result.columns:
                continue
            value = overlay_row[column]
            if is_empty(result.at[index, column]) and not is_empty(value):
                result.at[index, column] = str(value)
                changed += 1
    return result, changed


def autofill_parent_metadata(
    frame: pd.DataFrame,
    *,
    fetcher: Callable[[str], dict] = fetch_ena_sample_metadata,
) -> tuple[pd.DataFrame, dict, dict[str, dict]]:
    """Fill empty cells from unique parent samples and return a concise summary."""
    if "sample derived from" not in frame.columns:
        raise ValueError("Metadata table does not contain 'sample derived from'")
    result = frame.copy()
    parent_metadata = {}
    failures = {}
    fallback_accessions = []
    filled = 0
    accessions = sorted(
        {
            str(value).strip()
            for value in result["sample derived from"]
            if not is_empty(value)
        }
    )
    for accession in accessions:
        raw = fetcher(accession)
        if raw.get("_error"):
            failures[accession] = raw["_error"]
            continue
        parent_metadata[accession] = raw
        if raw.get("_metadata_source") == "ENA Portal linked-assembly fallback":
            fallback_accessions.append(accession)
        mapped = attrs_to_mag_metadata(raw)
        matching = result.index[
            result["sample derived from"].astype(str).str.strip() == accession
        ]
        for index in matching:
            for column, value in mapped.items():
                if column not in result.columns or is_empty(value):
                    continue
                if is_empty(result.at[index, column]):
                    result.at[index, column] = str(value)
                    filled += 1
    return result, {
        "parent_samples": len(accessions),
        "filled_cells": filled,
        "failed_accessions": failures,
        "fallback_accessions": fallback_accessions,
    }, parent_metadata
