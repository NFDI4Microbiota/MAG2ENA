#!/usr/bin/env python3
"""MAG2ENA command-line workflow for preparing and submitting MAGs to ENA."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import date, datetime, timedelta
import getpass
import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from zipfile import ZIP_DEFLATED, ZipFile

import pandas as pd
import polars as pl
import requests
from tqdm import tqdm

from App.utils.assembly import (
    ENA_SINGLE_CONTIG_HELP,
    chromosome_list_text,
    inspect_fasta_collection,
    require_single_contig_confirmation,
    single_contig_samples,
)
from App.utils.environment_resolver import (
    OlsEnvoClient,
    apply_environment_proposals,
    resolve_environment_batch,
)
from App.utils.metadata_assistance import (
    autofill_parent_metadata,
    is_empty,
    merge_empty_cells,
    parse_checkm_file,
    parse_gtdbtk_file,
)
from App.utils.taxonomy import resolve_species_taxon, validate_mag_taxonomy_records


VERSION = "1.0.7"
PROJECT_DIR = Path(__file__).resolve().parent
CHECKLIST_XML = PROJECT_DIR / "checklists" / "ERC000047.xml"
WEBIN_JAR = PROJECT_DIR / "App" / "webin-cli-9.0.3.jar"
WEBIN_USER_ENV = "MAG2ENA_WEBIN_USER"
WEBIN_PASSWORD_ENV = "MAG2ENA_WEBIN_PASSWORD"
ENVIRONMENT_DECISIONS = {"accept", "accepted", "approve", "approved", "yes", "true", "1"}


def load_fields_from_xml(xml_path: str | Path) -> dict:
    """Load field requirements from the ERC000047 XML checklist."""
    root = ET.parse(xml_path).getroot()
    fields = {
        "sample_name": {
            "label": "sample_name",
            "description": "Sample name (must match FASTA filename)",
            "type": "free",
            "regex": None,
            "enum": None,
            "mandatory": True,
        },
        "organism": {
            "label": "organism",
            "description": "Scientific organism name",
            "type": "free",
            "regex": None,
            "enum": None,
            "mandatory": True,
        },
        "tax_id": {
            "label": "tax_id",
            "description": "NCBI taxonomy ID",
            "type": "free",
            "regex": None,
            "enum": None,
            "mandatory": True,
        },
    }
    for field in root.findall(".//FIELD"):
        label = field.findtext("LABEL")
        if not label:
            continue
        description = field.findtext("DESCRIPTION")
        mandatory = field.findtext("MANDATORY") == "mandatory"
        regex_value = field.findtext(".//REGEX_VALUE")
        enum_values = [
            element.text
            for element in field.findall(".//TEXT_CHOICE_FIELD/TEXT_VALUE/VALUE")
            if element.text is not None
        ]
        if regex_value:
            field_type, compiled, choices = "regex", re.compile(regex_value), None
        elif enum_values:
            field_type, compiled, choices = "enum", None, enum_values
        else:
            field_type, compiled, choices = "free", None, None
        fields[label] = {
            "label": label,
            "description": description,
            "type": field_type,
            "regex": compiled,
            "enum": choices,
            "mandatory": mandatory,
        }
    fields["genome coverage"] = {
        "label": "genome coverage",
        "description": "Estimated sequencing depth",
        "type": "regex",
        "regex": re.compile(
            r"^(?:0?\.[0-9]*[1-9][0-9]*|[1-9][0-9]*(?:\.[0-9]+)?)$"
        ),
        "enum": None,
        "mandatory": True,
    }
    fields["platform"] = {
        "label": "platform",
        "description": "Sequencing platform",
        "type": "free",
        "regex": None,
        "enum": None,
        "mandatory": True,
    }
    return fields


def validate_dataframe(df: pl.DataFrame, field_defs: dict) -> None:
    """Run checklist and species-level MAG taxonomy preflight."""
    errors = []
    missing_columns = [
        label
        for label, field in field_defs.items()
        if field["mandatory"] and label not in df.columns
    ]
    if missing_columns:
        errors.append("Missing mandatory columns: " + ", ".join(missing_columns))

    indexed = df.with_row_index(offset=1)
    for column in df.columns:
        field = field_defs.get(column)
        if not field:
            continue
        values = pl.col(column).cast(pl.String)
        if field["mandatory"]:
            invalid = indexed.filter(
                values.is_null() | (values.str.strip_chars() == "")
            ).select("index").to_series().to_list()
            if invalid:
                errors.append(f"Column '{column}': empty at rows {invalid}")
        if field["type"] == "regex":
            invalid = indexed.filter(
                values.is_not_null()
                & (values.str.strip_chars() != "")
                & ~values.str.contains(field["regex"].pattern)
            ).select("index").to_series().to_list()
            if invalid:
                errors.append(f"Column '{column}': regex mismatch at rows {invalid}")
        if field["type"] == "enum":
            invalid = indexed.filter(
                values.is_not_null()
                & (values.str.strip_chars() != "")
                & ~values.is_in(field["enum"])
            ).select("index").to_series().to_list()
            if invalid:
                errors.append(f"Column '{column}': invalid choice at rows {invalid}")

    errors.extend(
        "Taxonomy at row "
        f"{error['row']} ({error['sample_name'] or 'unnamed sample'}, "
        f"tax_id {error['tax_id']}): {error['expected']}"
        for error in validate_mag_taxonomy_records(df.to_dicts())
    )
    if errors:
        raise ValueError("\n".join(errors))


def _pandas_from_polars(frame: pl.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(frame.to_dicts(), columns=frame.columns)


def _polars_from_pandas(frame: pd.DataFrame) -> pl.DataFrame:
    return pl.DataFrame(
        {
            column: [None if pd.isna(value) else str(value) for value in frame[column]]
            for column in frame.columns
        }
    )


def _resolve_missing_taxonomy(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    result = frame.copy()
    resolved = {}
    unresolved = []
    for index, row in result.iterrows():
        organism = str(row.get("organism", "") or "").strip()
        tax_id = row.get("tax_id", "")
        if not organism or not is_empty(tax_id):
            continue
        match = resolve_species_taxon(organism)
        sample_name = str(row.get("sample_name", index))
        if match:
            result.at[index, "tax_id"] = match["tax_id"]
            resolved[sample_name] = match
        else:
            unresolved.append(sample_name)
    return result, {"resolved": resolved, "unresolved": unresolved}


def _proposal_target(proposal: dict) -> str:
    parent = proposal.get("parent_accession") or "no parent"
    samples = ", ".join(map(str, proposal.get("sample_names", [])))
    return f"{parent} | {proposal['field']} | {samples}"


def write_environment_review(
    proposals: list[dict],
    output_path: str | Path,
    *,
    decisions: dict[tuple[str, str], str] | None = None,
) -> Path:
    """Write a concise four-column curator review table."""
    decision_lookup = decisions or {}
    rows = []
    for proposal in proposals:
        target = _proposal_target(proposal)
        term = proposal.get("selected_value", "")
        rows.append(
            {
                "target": target,
                "source evidence": proposal.get("evidence", ""),
                "proposed term": term,
                "decision": decision_lookup.get((target, term), ""),
            }
        )
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def read_environment_decisions(path: str | Path) -> dict[tuple[str, str], str]:
    """Read curator decisions keyed by the review target and proposed term."""
    frame = pd.read_csv(path, sep="\t", keep_default_na=False)
    required = {"target", "proposed term", "decision"}
    if not required.issubset(frame.columns):
        raise ValueError(
            "Environmental decisions TSV must contain target, proposed term, and decision columns."
        )
    return {
        (str(row["target"]), str(row["proposed term"])): str(row["decision"])
        for _, row in frame.iterrows()
    }


def _apply_environment_decisions(
    frame: pd.DataFrame,
    proposals: list[dict],
    decisions: dict[tuple[str, str], str],
    *,
    overwrite: bool,
) -> tuple[pd.DataFrame, list[dict], int]:
    accepted = []
    for proposal in proposals:
        decision = decisions.get(
            (_proposal_target(proposal), proposal.get("selected_value", "")), ""
        )
        if str(decision).strip().casefold() in ENVIRONMENT_DECISIONS:
            accepted.append(proposal["proposal_id"])
    result, audit = apply_environment_proposals(
        frame, proposals, accepted, overwrite=overwrite
    )
    return result, audit, len(accepted)


def prepare_metadata(frame: pl.DataFrame, args: argparse.Namespace) -> tuple[pl.DataFrame, dict]:
    """Apply opt-in web-equivalent metadata assistance without provenance."""
    prepared = _pandas_from_polars(frame)
    summary: dict[str, object] = {}

    if args.checkm:
        prepared, changed = merge_empty_cells(prepared, parse_checkm_file(args.checkm))
        summary["quality_cells_filled"] = changed
        print(f"CheckM/CheckM2 import filled {changed} cell(s).")
    if args.gtdbtk:
        prepared, changed = merge_empty_cells(prepared, parse_gtdbtk_file(args.gtdbtk))
        summary["organism_cells_filled"] = changed
        print(f"GTDB-Tk import filled {changed} cell(s).")

    needs_parent_records = bool(
        args.autofill_parent_metadata
        or args.environment_review_out
        or args.environment_decisions
    )
    parent_metadata = {}
    if needs_parent_records:
        autofilled, parent_summary, parent_metadata = autofill_parent_metadata(prepared)
        summary["parent_metadata"] = parent_summary
        if args.autofill_parent_metadata:
            prepared = autofilled
            print(
                "ENA parent autofill filled "
                f"{parent_summary['filled_cells']} cell(s) from "
                f"{parent_summary['parent_samples']} parent sample(s)."
            )
            if parent_summary["fallback_accessions"]:
                print(
                    "Linked-assembly fallback used for: "
                    + ", ".join(parent_summary["fallback_accessions"])
                )
            for accession, error in parent_summary["failed_accessions"].items():
                print(f"WARNING: Parent metadata unavailable for {accession}: {error}")

    if args.resolve_taxonomy:
        prepared, taxonomy_summary = _resolve_missing_taxonomy(prepared)
        summary["taxonomy"] = taxonomy_summary
        print(
            f"Taxonomy resolver filled {len(taxonomy_summary['resolved'])} tax ID(s); "
            f"{len(taxonomy_summary['unresolved'])} remained unresolved."
        )

    proposals = []
    decisions = {}
    if args.environment_decisions:
        decisions = read_environment_decisions(args.environment_decisions)
    if args.environment_review_out or args.environment_decisions:
        proposals = resolve_environment_batch(
            prepared,
            OlsEnvoClient(),
            empty_only=args.environment_scope == "empty",
            parent_metadata_by_accession=parent_metadata,
        )
        summary["environment_proposals"] = len(proposals)
        if decisions:
            prepared, audit, accepted = _apply_environment_decisions(
                prepared,
                proposals,
                decisions,
                overwrite=args.environment_scope == "all",
            )
            summary["environment_decisions_accepted"] = accepted
            summary["environment_cells_changed"] = sum(
                len(str(item.get("rows", "")).split(", ")) for item in audit
            )
            print(
                f"Applied {accepted} approved environmental proposal(s), changing "
                f"{summary['environment_cells_changed']} cell(s)."
            )
        if args.environment_review_out:
            review_path = Path(args.environment_review_out)
            decisions_path = Path(args.environment_decisions) if args.environment_decisions else None
            if decisions_path is None or review_path.resolve() != decisions_path.resolve():
                write_environment_review(proposals, review_path, decisions=decisions)
                print(f"Environmental review written to {review_path}.")

    output_path = args.prepared_metadata_out
    if args.prepare_only and not output_path:
        output_path = Path("mag2ena_prepared_metadata.tsv")
    if output_path:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        prepared.to_csv(path, sep="\t", index=False)
        summary["prepared_metadata"] = str(path.resolve())
        print(f"Prepared metadata written to {path}.")
    return _polars_from_pandas(prepared), summary


def _validate_fasta_headers(fasta_map: dict[str, Path], max_len: int = 50) -> list[dict]:
    violations = []
    for sample_name, path in fasta_map.items():
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                sequence_name = line[1:].split()[0] if line[1:].strip() else ""
                if len(sequence_name) > max_len:
                    violations.append(
                        {
                            "sample": sample_name,
                            "sequence name": sequence_name[:80]
                            + ("…" if len(sequence_name) > 80 else ""),
                            "length": len(sequence_name),
                        }
                    )
    return violations


def _deduplicate_fasta_headers(src_path: Path) -> tuple[Path, bool]:
    seen = set()
    with gzip.open(src_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            raw_name = line[1:].split()[0] if line[1:].strip() else ""
            ena_name = raw_name.split("|")[0]
            if ena_name in seen:
                break
            seen.add(ena_name)
        else:
            return src_path, False

    descriptor, temporary_name = tempfile.mkstemp(suffix=".fasta.gz")
    os.close(descriptor)
    destination = Path(temporary_name)
    counts: dict[str, int] = {}
    with gzip.open(src_path, "rt", encoding="utf-8") as source, gzip.open(
        destination, "wt", encoding="utf-8"
    ) as output:
        for line in source:
            if not line.startswith(">"):
                output.write(line)
                continue
            raw_name = line[1:].split()[0] if line[1:].strip() else ""
            rest = line[1 + len(raw_name):]
            ena_name = raw_name.split("|")[0]
            count = counts.get(ena_name, 0)
            counts[ena_name] = count + 1
            if count == 0:
                output.write(line)
            elif "|" in raw_name:
                prefix, suffix = raw_name.split("|", 1)
                output.write(f">{prefix}_{count}|{suffix}{rest}")
            else:
                output.write(f">{raw_name}_{count}{rest}")
    return destination, True


def collect_fastas(fasta_dir: Path) -> dict[str, Path]:
    return {
        path.name.removesuffix(".fasta.gz"): path.resolve()
        for path in fasta_dir.glob("*.fasta.gz")
    }


def resolve_webin_credentials(
    args: argparse.Namespace,
    *,
    environ: dict[str, str] | os._Environ = os.environ,
    input_fn=input,
    getpass_fn=getpass.getpass,
) -> tuple[str, str]:
    """Read credentials from environment variables or non-echoing prompts."""
    username = str(environ.get(WEBIN_USER_ENV, "")).strip()
    password = str(environ.get(WEBIN_PASSWORD_ENV, ""))
    if not username:
        if args.non_interactive:
            raise ValueError(f"Set {WEBIN_USER_ENV} for non-interactive submission.")
        username = input_fn("ENA Webin username: ").strip()
    if not password:
        if args.non_interactive:
            raise ValueError(f"Set {WEBIN_PASSWORD_ENV} for non-interactive submission.")
        password = getpass_fn("ENA Webin password: ")
    if not username or not password:
        raise ValueError("ENA Webin username and password are required.")
    return username, password


def _java_argfile_value(value: object) -> str:
    text = str(value)
    if "\n" in text or "\r" in text:
        raise ValueError("Webin-CLI arguments cannot contain line breaks.")
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


@contextmanager
def secure_java_argfile(arguments: list[object]):
    """Store sensitive Java arguments briefly in a mode-0600 response file."""
    descriptor, name = tempfile.mkstemp(prefix="mag2ena-webin-", suffix=".args")
    path = Path(name)
    try:
        os.fchmod(descriptor, 0o600)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            for argument in arguments:
                handle.write(_java_argfile_value(argument) + "\n")
        yield path
    finally:
        path.unlink(missing_ok=True)


def _poll_receipt(poll_url: str, auth: tuple[str, str], timeout: int = 900) -> str:
    deadline = time.monotonic() + timeout
    while True:
        response = requests.get(poll_url, auth=auth, timeout=30)
        if response.status_code == 200:
            return response.text
        if response.status_code != 202:
            response.raise_for_status()
        if time.monotonic() >= deadline:
            raise TimeoutError("Timed out waiting for the ENA metadata receipt.")
        time.sleep(5)


def _add_sample_attribute(parent: ET.Element, tag: str, value: object, units=None):
    attribute = ET.SubElement(parent, "SAMPLE_ATTRIBUTE")
    ET.SubElement(attribute, "TAG").text = str(tag)
    ET.SubElement(attribute, "VALUE").text = str(value)
    if units:
        ET.SubElement(attribute, "UNITS").text = units


def _build_metadata_xml(
    rows: list[dict], submission: dict, project_accession: str | None
) -> bytes:
    root = ET.Element("WEBIN")
    submission_set = ET.SubElement(root, "SUBMISSION_SET")
    submission_element = ET.SubElement(submission_set, "SUBMISSION")
    actions = ET.SubElement(submission_element, "ACTIONS")
    add_action = ET.SubElement(actions, "ACTION")
    ET.SubElement(add_action, "ADD")
    if submission.get("hold_until_date"):
        hold_action = ET.SubElement(actions, "ACTION")
        ET.SubElement(
            hold_action,
            "HOLD",
            {"HoldUntilDate": submission["hold_until_date"]},
        )
    if not project_accession:
        project_set = ET.SubElement(root, "PROJECT_SET")
        project = ET.SubElement(project_set, "PROJECT", {"alias": submission["study_name"]})
        ET.SubElement(project, "TITLE").text = submission["study_title"]
        ET.SubElement(project, "DESCRIPTION").text = submission["study_description"]
        submission_project = ET.SubElement(project, "SUBMISSION_PROJECT")
        ET.SubElement(submission_project, "SEQUENCING_PROJECT")

    sample_set = ET.SubElement(root, "SAMPLE_SET")
    explicit_columns = {
        "metagenomic source": None,
        "sample derived from": None,
        "project name": None,
        "completeness score": "%",
        "completeness software": None,
        "contamination score": "%",
        "binning software": None,
        "assembly quality": None,
        "binning parameters": None,
        "taxonomic identity marker": None,
        "isolation_source": None,
        "collection date": None,
        "geographic location (latitude)": "DD",
        "geographic location (longitude)": "DD",
        "broad-scale environmental context": None,
        "local environmental context": None,
        "environmental medium": None,
        "geographic location (country and/or sea)": None,
        "assembly software": None,
    }
    reserved = {"sample_name", "organism", "tax_id", "platform", "genome coverage"}
    for row in rows:
        sample = ET.SubElement(
            sample_set,
            "SAMPLE",
            {"alias": str(row["sample_name"]), "center_name": ""},
        )
        ET.SubElement(sample, "TITLE").text = (
            "This sample represents a MAG derived from the metagenomic sample "
            + str(row["sample derived from"])
        )
        sample_name = ET.SubElement(sample, "SAMPLE_NAME")
        ET.SubElement(sample_name, "TAXON_ID").text = str(row["tax_id"])
        ET.SubElement(sample_name, "SCIENTIFIC_NAME").text = str(row["organism"])
        attributes = ET.SubElement(sample, "SAMPLE_ATTRIBUTES")
        for column, units in explicit_columns.items():
            _add_sample_attribute(attributes, column, row[column], units)
        for column, value in row.items():
            if column in reserved or column in explicit_columns or column == "ENA-CHECKLIST":
                continue
            if value is None or str(value).strip() == "":
                continue
            _add_sample_attribute(attributes, column, value)
        _add_sample_attribute(attributes, "ENA-CHECKLIST", "ERC000047")
    return ET.tostring(root, encoding="UTF-8", xml_declaration=True)


def _safe_output_directory(path: Path) -> Path:
    target = path.resolve()
    if target.exists() and any(target.iterdir()):
        raise ValueError(f"Output directory is not empty: {target}")
    target.mkdir(parents=True, exist_ok=True)
    return target


def build_and_submit(
    df: pl.DataFrame,
    submission: dict,
    fasta_map: dict[str, Path],
    output_dir: Path,
) -> dict:
    """Register metadata and submit each assembly serially through Webin-CLI."""
    validate_errors = validate_mag_taxonomy_records(df.to_dicts())
    if validate_errors:
        raise ValueError(
            "Taxonomy preflight failed immediately before submission: "
            + "; ".join(error["expected"] for error in validate_errors)
        )
    output = _safe_output_directory(output_dir)
    df.write_csv(output / "metadata.submitted.tsv", separator="\t")
    inspections = inspect_fasta_collection(fasta_map)
    require_single_contig_confirmation(
        inspections, bool(submission.get("single_contig_chromosomes_confirmed"))
    )

    submitted = 0
    errors = 0
    project_accession = submission.get("study_accession")
    auth = (submission["ena_user"], submission["ena_password"])
    metadata_url = (
        "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/submit/queue"
        if submission["portal"] == "test"
        else "https://www.ebi.ac.uk/ena/submit/webin-v2/submit/queue"
    )
    success_log = output / "success.txt"
    error_log = output / "error.txt"
    batch_size = 1_000
    manifest_zip_path = output / "manifests.zip"

    with ZipFile(manifest_zip_path, "w", ZIP_DEFLATED) as manifest_zip:
        for offset in tqdm(
            range(0, len(df), batch_size),
            desc=f"Processing ENA metadata batches ({batch_size} samples per batch)",
        ):
            rows = df.slice(offset, batch_size).to_dicts()
            xml_bytes = _build_metadata_xml(rows, submission, project_accession)
            xml_name = "submit.xml" if offset == 0 else f"submit_{offset}.xml"
            (output / xml_name).write_bytes(xml_bytes)
            response = requests.post(
                metadata_url,
                headers={"Accept": "application/json"},
                files={"file": (xml_name, xml_bytes, "text/xml")},
                auth=auth,
                timeout=60,
            )
            response.raise_for_status()
            payload = response.json()
            poll_url = payload.get("_links", {}).get("poll", {}).get("href")
            if not poll_url:
                raise ValueError("ENA metadata response did not include a poll URL.")
            receipt_text = _poll_receipt(poll_url, auth)
            receipt_name = "webin_log.xml" if offset == 0 else f"webin_log_{offset}.xml"
            (output / receipt_name).write_text(receipt_text, encoding="utf-8")
            receipt = ET.fromstring(receipt_text)
            alias_to_accession = {
                sample.get("alias"): sample.get("accession")
                for sample in receipt.findall(".//SAMPLE")
                if sample.get("alias") and sample.get("accession")
            }
            existing = re.compile(r'alias: "(.+?)".*accession: "(.+?)"')
            for error in receipt.findall(".//ERROR"):
                message = error.text or "Unknown ENA metadata error"
                match = existing.search(message)
                if match:
                    alias_to_accession[match.group(1)] = match.group(2)
                else:
                    with error_log.open("a", encoding="utf-8") as handle:
                        handle.write(f"METADATA ERROR: {message}\n")
            if not project_accession:
                project = receipt.find(".//PROJECT")
                project_accession = project.get("accession") if project is not None else None
            if not project_accession:
                raise ValueError("ENA did not return a study accession.")

            for alias, sample_accession in alias_to_accession.items():
                original_fasta = fasta_map.get(alias)
                if original_fasta is None:
                    continue
                fasta_path = original_fasta
                deduplicated = None
                chromosome_path = None
                manifest_path = None
                try:
                    inspection = inspections[alias]
                    if inspection["sequence_count"] > 1:
                        fasta_path, changed = _deduplicate_fasta_headers(original_fasta)
                        if changed:
                            deduplicated = fasta_path
                    row = df.filter(pl.col("sample_name") == alias).row(0, named=True)
                    manifest = (
                        f"STUDY\t{project_accession}\n"
                        f"SAMPLE\t{sample_accession}\n"
                        f"ASSEMBLYNAME\t{alias}\n"
                        "ASSEMBLY_TYPE\tMetagenome-Assembled Genome (MAG)\n"
                        f"COVERAGE\t{row['genome coverage']}\n"
                        f"PROGRAM\t{row['assembly software']}\n"
                        f"PLATFORM\t{row['platform']}\n"
                        f"FASTA\t{fasta_path}"
                    )
                    if inspection["sequence_count"] == 1:
                        descriptor, name = tempfile.mkstemp(suffix=".gz")
                        os.close(descriptor)
                        chromosome_path = Path(name)
                        with gzip.open(chromosome_path, "wt", encoding="utf-8") as handle:
                            handle.write(chromosome_list_text(inspection))
                        manifest += f"\nCHROMOSOME_LIST\t{chromosome_path}"
                    manifest_zip.writestr(f"{alias}.manifest", manifest + "\n")
                    descriptor, name = tempfile.mkstemp(prefix="mag2ena-manifest-")
                    os.fchmod(descriptor, 0o600)
                    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
                        handle.write(manifest)
                    manifest_path = Path(name)
                    java_arguments = [
                        "-jar",
                        WEBIN_JAR,
                        "-username",
                        auth[0],
                        "-password",
                        auth[1],
                        "-context",
                        "genome",
                        "-manifest",
                        manifest_path,
                        "-submit",
                    ]
                    if submission["portal"] == "test":
                        java_arguments.append("-test")
                    with secure_java_argfile(java_arguments) as argument_file:
                        result = subprocess.run(
                            ["java", f"@{argument_file}"],
                            capture_output=True,
                            text=True,
                            check=False,
                        )
                    combined = "\n".join(
                        part for part in (result.stdout, result.stderr) if part
                    )
                    combined = combined.replace(auth[0], "[REDACTED]").replace(
                        auth[1], "[REDACTED]"
                    )
                    successful = result.returncode == 0 and "successfully" in combined.lower()
                    destination = success_log if successful else error_log
                    with destination.open("a", encoding="utf-8") as handle:
                        handle.write(f"SAMPLE : {alias}\n{combined}\n\n")
                    if successful:
                        submitted += 1
                    else:
                        errors += 1
                finally:
                    if manifest_path:
                        manifest_path.unlink(missing_ok=True)
                    if chromosome_path:
                        chromosome_path.unlink(missing_ok=True)
                    if deduplicated:
                        deduplicated.unlink(missing_ok=True)

    result_summary = {
        "submitted": submitted,
        "errors": errors,
        "project_accession": project_accession,
        "webin_cli": "9.0.3",
        "portal": submission["portal"],
    }
    (output / "result.json").write_text(
        json.dumps(result_summary, indent=2), encoding="utf-8"
    )
    print(f"Submitted successfully: {submitted}; errors: {errors}")
    print(f"Results written to: {output}")
    return result_summary


def _validate_release_date(value: str) -> str:
    try:
        parsed = date.fromisoformat(value)
    except ValueError as exc:
        raise ValueError("Release date must use YYYY-MM-DD.") from exc
    if parsed < date.today() or parsed > date.today() + timedelta(days=730):
        raise ValueError("Release date must be between today and 730 days from today.")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Curate, validate, and submit metagenome-assembled genomes to the "
            "European Nucleotide Archive."
        )
    )
    parser.add_argument("--version", action="version", version=f"MAG2ENA {VERSION}")
    parser.add_argument("--metadata", required=True, type=Path, help="Input metadata TSV.")
    parser.add_argument(
        "--fasta-dir",
        type=Path,
        help="Directory containing one <sample_name>.fasta.gz file per MAG.",
    )
    parser.add_argument(
        "--checkm", type=Path, help="Optional CheckM storage.tsv or CheckM2 quality_report.tsv."
    )
    parser.add_argument(
        "--gtdbtk", type=Path, help="Optional GTDB-Tk bacterial or archaeal summary TSV."
    )
    parser.add_argument(
        "--autofill-parent-metadata",
        action="store_true",
        help="Fill empty parent-derived cells from ENA, with strict linked-assembly fallback.",
    )
    parser.add_argument(
        "--resolve-taxonomy",
        action="store_true",
        help="Resolve missing tax IDs from organism names, accepting species-rank taxa only.",
    )
    parser.add_argument(
        "--environment-review-out",
        type=Path,
        help="Write four-column environmental proposals for curator decisions.",
    )
    parser.add_argument(
        "--environment-decisions",
        type=Path,
        help="Apply proposals marked accept/yes/true in a prior review TSV.",
    )
    parser.add_argument(
        "--environment-scope",
        choices=["empty", "all"],
        default="empty",
        help="Resolve empty fields only (default) or normalize all three environmental fields.",
    )
    parser.add_argument(
        "--prepared-metadata-out",
        type=Path,
        help="Write metadata after assistance and curator decisions.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare metadata and review files without requesting credentials or submitting.",
    )
    parser.add_argument(
        "--non-interactive",
        action="store_true",
        help=f"Require credentials in {WEBIN_USER_ENV}/{WEBIN_PASSWORD_ENV}; never prompt.",
    )
    parser.add_argument(
        "--portal",
        choices=["test", "prod"],
        default="test",
        help="ENA portal (default: test).",
    )
    parser.add_argument("--study-accession", help="Existing ENA study accession.")
    parser.add_argument("--study-name", help="Alias for a new study.")
    parser.add_argument("--study-title", help="New study title (20-250 characters).")
    parser.add_argument(
        "--study-description", help="New study description (20-4,000 characters)."
    )
    default_release = (date.today() + timedelta(days=729)).isoformat()
    parser.add_argument(
        "--release-date",
        default=default_release,
        help=f"Hold-until date in YYYY-MM-DD (default: {default_release}).",
    )
    parser.add_argument(
        "--single-contig-as-chromosome",
        action="store_true",
        help=(
            "Confirm every one-sequence MAG is a complete chromosome. Incomplete "
            "single-contig MAGs require an ENA exception and remain blocked."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="New or empty results directory (default: timestamped directory).",
    )
    return parser


def _validate_study_arguments(args: argparse.Namespace) -> None:
    if args.study_accession:
        return
    if not args.study_name:
        raise ValueError("Study name is required when creating a new study.")
    if not 20 <= len(args.study_title or "") <= 250:
        raise ValueError("Study title must be between 20 and 250 characters.")
    if not 20 <= len(args.study_description or "") <= 4_000:
        raise ValueError("Study description must be between 20 and 4,000 characters.")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if not args.metadata.is_file():
            raise ValueError(f"Metadata file not found: {args.metadata}")
        frame = pl.read_csv(args.metadata, separator="\t", infer_schema_length=0)
        frame, _ = prepare_metadata(frame, args)
        if args.prepare_only:
            print("Metadata preparation complete; no ENA submission was attempted.")
            return 0

        if args.fasta_dir is None or not args.fasta_dir.is_dir():
            raise ValueError("--fasta-dir is required for submission and must exist.")
        _validate_study_arguments(args)
        release_date = _validate_release_date(args.release_date)
        validate_dataframe(frame, load_fields_from_xml(CHECKLIST_XML))

        fasta_map = collect_fastas(args.fasta_dir)
        expected = {str(value) for value in frame["sample_name"]}
        missing = expected - set(fasta_map)
        if missing:
            raise ValueError("Missing FASTA files for: " + ", ".join(sorted(missing)))
        fasta_map = {name: path for name, path in fasta_map.items() if name in expected}
        violations = _validate_fasta_headers(fasta_map)
        if violations:
            details = "; ".join(
                f"{item['sample']}: {item['sequence name']} ({item['length']} chars)"
                for item in violations
            )
            raise ValueError(f"FASTA sequence names exceed the ENA 50-character limit: {details}")
        inspections = inspect_fasta_collection(fasta_map)
        detected_single_contigs = single_contig_samples(inspections)
        require_single_contig_confirmation(
            inspections, args.single_contig_as_chromosome
        )
        if detected_single_contigs:
            print(
                "Confirmed complete chromosome-level one-sequence MAGs: "
                + ", ".join(detected_single_contigs)
            )
            print(ENA_SINGLE_CONTIG_HELP)

        username, password = resolve_webin_credentials(args)
        output_dir = args.output_dir or Path(
            "mag2ena_cli_results_" + datetime.now().strftime("%Y%m%d_%H%M%S")
        )
        result = build_and_submit(
            frame,
            {
                "study_accession": args.study_accession,
                "study_name": args.study_name,
                "study_title": args.study_title,
                "study_description": args.study_description,
                "ena_user": username,
                "ena_password": password,
                "portal": args.portal,
                "hold_until_date": release_date,
                "single_contig_chromosomes_confirmed": args.single_contig_as_chromosome,
            },
            fasta_map,
            output_dir,
        )
        return 0 if result["errors"] == 0 else 1
    except (OSError, ValueError, requests.RequestException, ET.ParseError) as exc:
        print(f"ERROR: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
