import streamlit as st
import pandas as pd
import polars as pl
import re
import xml.etree.ElementTree as ET
import gzip
import tempfile
import zipfile
from pathlib import Path
import subprocess
import os
import requests
import time
import json
from datetime import datetime, date, timedelta
from urllib.parse import quote
from rq import get_current_job
from utils.tasks import enqueue_task, manager, scrub_job_credentials
from utils.db import TaskResultManager, TaskStatus
from utils.encryption import encrypt_job_folder
from utils.ena_metadata import fetch_ena_sample_metadata
from utils.geography import render_geographic_overview
from utils.assembly import (
    ENA_SINGLE_CONTIG_HELP,
    chromosome_list_text,
    inspect_fasta_collection,
    require_single_contig_confirmation,
    single_contig_samples,
)
from utils.environment_resolver import (
    OlsEnvoClient,
    apply_environment_proposals,
    build_resolution_groups,
    extract_sample_accession as _extract_sample_accession,
    proposals_to_review_dataframe,
    resolve_environment_batch,
)
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
from utils.taxonomy import (
    fetch_ena_taxon,
    mag_taxon_validation_error,
    validate_mag_taxonomy_records,
)

# =========================================================
# CONFIG
# =========================================================

_HERE = Path(__file__).parent
_APP_DIR = _HERE.parent
_PROJECT_DIR = _APP_DIR.parent

CHECKLIST_XML = str(_PROJECT_DIR / "checklists" / "ERC000047.xml")
WEBIN_JAR = str(_APP_DIR / "webin-cli-9.0.3.jar")
SUBMISSIONS_DIR = _APP_DIR / "jobs"
EXAMPLES_DIR = _PROJECT_DIR / "examples"
MAG2ENA_VERSION = "1.0.7"

# =========================================================
# METADATA HELPERS
# =========================================================

_DATE_REGEX = re.compile(
    r'^[12][0-9]{3}(-(0[1-9]|1[0-2])(-(0[1-9]|[12][0-9]|3[01])'
    r'(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?'
    r'(/[0-9]{4}(-[0-9]{2}(-[0-9]{2}'
    r'(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?)?$'
)


@st.cache_data(ttl=300, show_spinner=False)
def _fetch_ena_sample_metadata(sample_acc: str) -> dict:
    """Fetch a sample record, with strict linked-assembly recovery if needed."""
    return fetch_ena_sample_metadata(sample_acc)


def _fix_coordinate_string(coord_str: str):
    """Parse lat_lon attribute (e.g. '1.5 N 40.2 E') into (lat_str, lon_str)."""
    if not coord_str:
        return None, None
    coord_str = str(coord_str).strip()
    match = re.search(
        r'(\d*\.?\d+)\s*([NS])\s+(\d*\.?\d+)\s*([EW])',
        coord_str,
        re.IGNORECASE,
    )
    if match:
        lat = float(match.group(1)) * (-1 if match.group(2).upper() == "S" else 1)
        lon = float(match.group(3)) * (-1 if match.group(4).upper() == "W" else 1)
        return str(round(lat, 6)), str(round(lon, 6))
    # Bare "lat lon" without direction letters
    parts = coord_str.split()
    if len(parts) == 2:
        try:
            return str(float(parts[0])), str(float(parts[1]))
        except ValueError:
            pass
    return None, None


def _normalize_date(date_str: str) -> str:
    """Normalise a date string to the ISO 8601 format ENA expects."""
    if not date_str or date_str.strip() in ("", "missing"):
        return ""
    d = re.sub(r'^(\d{4})$', r'\1-01-01', date_str.strip())
    d = re.sub(r'^(\d{4})-(\d{2})$', r'\1-\2-01', d)
    d = re.sub(r'^(\d{4})(\d{2})(\d{2})$', r'\1-\2-\3', d)
    return d if _DATE_REGEX.fullmatch(d) else ""


def _attrs_to_mag_metadata(attrs: dict) -> dict:
    """Map raw ENA sample attributes to MAG metadata column names."""
    env_broad = attrs.get("env_broad_scale") or attrs.get("env_biome") or ""
    env_local = attrs.get("env_local_scale") or attrs.get("env_feature") or ""
    env_medium = attrs.get("env_medium") or attrs.get("env_material") or ""

    if not env_broad and not env_local and not env_medium:
        gold = attrs.get("GOLD Ecosystem Classification", "")
        if gold:
            parts = gold.split(" | ")
            env_broad = parts[1] if len(parts) > 1 else ""
            env_local = parts[2] if len(parts) > 2 else ""
            env_medium = parts[3] if len(parts) > 3 else ""

    geo_loc_name = attrs.get("geo_loc_name", "")
    country = geo_loc_name.split(":")[0].strip() if geo_loc_name else ""

    lat, lon = _fix_coordinate_string(attrs.get("lat_lon", ""))

    return {
        "metagenomic source": attrs.get("organism", ""),
        "isolation_source": attrs.get("isolation_source", ""),
        "collection date": _normalize_date(attrs.get("collection_date", "")),
        "geographic location (latitude)": lat or "",
        "geographic location (longitude)": lon or "",
        "broad-scale environmental context": env_broad,
        "local environmental context": env_local,
        "environmental medium": env_medium,
        "geographic location (country and/or sea)": country,
    }


@st.cache_data(ttl=300, show_spinner=False)
def _ena_taxonomy_search(query: str) -> list[dict]:
    if len(query.strip()) < 3:
        return []
    try:
        r = requests.get(
            f"https://www.ebi.ac.uk/ena/taxonomy/rest/suggest-for-search/{quote(query.strip())}",
            timeout=6,
        )
        r.raise_for_status()
        return r.json()
    except Exception:
        return []


@st.cache_data(ttl=300, show_spinner=False)
def _fetch_ena_taxon(tax_id: str) -> dict:
    """Return taxon info from ENA taxonomy API for a given tax ID, or {} on failure."""
    return fetch_ena_taxon(tax_id)


def _parse_checkm_file(file) -> pd.DataFrame | None:
    """Parse CheckM or CheckM2 output TSV into sample_name + quality columns."""
    try:
        df = pd.read_csv(file, sep="\t")
        cols = df.columns.tolist()
        if "Name" in cols and "Completeness" in cols and "Contamination" in cols:
            name_col, software = "Name", "CheckM2"
        elif "Bin Id" in cols and "Completeness" in cols and "Contamination" in cols:
            name_col, software = "Bin Id", "CheckM"
        else:
            return None
        return pd.DataFrame({
            "sample_name": df[name_col].astype(str),
            "completeness score": df["Completeness"].round(2).astype(str),
            "contamination score": df["Contamination"].round(2).astype(str),
            "completeness software": software,
        })
    except Exception:
        return None


def _parse_gtdbtk_file(file) -> pd.DataFrame | None:
    """Parse GTDB-Tk summary TSV into sample_name + organism columns."""
    try:
        df = pd.read_csv(file, sep="\t")
        if "user_genome" not in df.columns or "classification" not in df.columns:
            return None

        def _gtdb_to_name(classification: str) -> str:
            parts = {}
            for p in classification.split(";"):
                if "__" in p:
                    rank, name = p.split("__", 1)
                    parts[rank.strip()] = name.strip().replace("_", " ")
            for rank in ("s", "g", "f", "o", "c", "p"):
                name = parts.get(rank, "")
                if name:
                    return name if rank == "s" else f"uncultured {name} bacterium"
            return "uncultured bacterium"

        return pd.DataFrame({
            "sample_name": df["user_genome"].astype(str),
            "organism": df["classification"].apply(_gtdb_to_name),
        })
    except Exception:
        return None


def _merge_into_metadata(base: pd.DataFrame, overlay: pd.DataFrame) -> pd.DataFrame:
    """Merge overlay columns into base, keyed on sample_name. Only fills empty cells."""
    result = base.copy()
    for col in overlay.columns:
        if col == "sample_name" or col not in result.columns:
            continue
        for _, orow in overlay.iterrows():
            mask = result["sample_name"] == orow["sample_name"]
            if not mask.any():
                continue
            idx = result.index[mask][0]
            current = result.at[idx, col]
            if pd.isna(current) or str(current).strip() in ("", "<NA>", "None", "nan"):
                result.at[idx, col] = str(orow[col])
    return result


# =========================================================
# XML CHECKLIST
# =========================================================

def _field(label, description, mandatory, ftype, *, regex=None, enum=None):
    return {"label": label, "description": description, "type": ftype,
            "regex": regex, "enum": enum, "mandatory": mandatory}


def load_fields_from_xml(xml_path: str, *, include_assembly_fields: bool = True) -> dict:

    tree = ET.parse(xml_path)
    root = tree.getroot()

    fields = {
        "sample_name": _field("sample_name", "Sample name (must match fasta filename)", True, "free"),
        "organism":    _field("organism",    "Scientific organism name",                True, "free"),
        "tax_id":      _field("tax_id",      "NCBI taxonomy ID",                        True, "free"),
    }

    for field in root.findall(".//FIELD"):

        label = field.findtext("LABEL")
        description = field.findtext("DESCRIPTION")
        mandatory = field.findtext("MANDATORY") == "mandatory"

        regex_elem = field.findtext(".//REGEX_VALUE")
        enum_elems = field.findall(".//TEXT_CHOICE_FIELD/TEXT_VALUE/VALUE")

        if regex_elem:
            fields[label] = _field(label, description, mandatory, "regex",
                                   regex=re.compile(regex_elem))
        elif enum_elems:
            fields[label] = _field(label, description, mandatory, "enum",
                                   enum=[e.text for e in enum_elems])
        else:
            fields[label] = _field(label, description, mandatory, "free")

    if include_assembly_fields:
        fields["genome coverage"] = _field(
            "genome coverage", "Estimated sequencing depth", True, "regex",
            regex=re.compile(r"^(?:0?\.[0-9]*[1-9][0-9]*|[1-9][0-9]*(?:\.[0-9]+)?)$"),
        )
        fields["platform"] = _field("platform", "Sequencing platform", True, "free")

    return fields

# =========================================================
# VALIDATION
# =========================================================

def validate_dataframe(df: pd.DataFrame, field_defs: dict) -> pd.DataFrame:

    errors = []

    for col in df.columns:

        if col not in field_defs:
            continue

        field = field_defs[col]

        for idx, value in df[col].items():

            value_str = "" if pd.isna(value) else str(value).strip()

            if field["mandatory"] and value_str == "":
                errors.append({
                    "row": idx + 1,
                    "field": col,
                    "value": "",
                    "expected": "Mandatory field — cannot be empty"
                })
                continue

            if value_str == "":
                continue

            if field["type"] == "regex":

                if not field["regex"].fullmatch(value_str):
                    errors.append({
                        "row": idx + 1,
                        "field": col,
                        "value": value_str,
                        "expected": f"Pattern: {field['regex'].pattern}"
                    })

            elif field["type"] == "enum":

                if value_str not in field["enum"]:
                    errors.append({
                        "row": idx + 1,
                        "field": col,
                        "value": value_str,
                        "expected": ", ".join(field["enum"])
                    })

            if col == "tax_id" and value_str:
                taxon = _fetch_ena_taxon(value_str)
                expected = mag_taxon_validation_error(value_str, taxon)
                if expected:
                    errors.append({
                        "row": idx + 1,
                        "field": col,
                        "value": value_str,
                        "expected": expected,
                    })

    return pd.DataFrame(errors)

# =========================================================
# STREAMLIT TABLE CONFIG
# =========================================================

def build_column_config(field_defs: dict):

    column_config = {}

    for name, field in field_defs.items():

        label = (
            f"{field['label']} *"
            if field["mandatory"]
            else field["label"]
        )

        if field["type"] == "enum":

            column_config[name] = st.column_config.SelectboxColumn(
                label=label,
                options=field["enum"],
                help=field["description"]
            )

        else:

            column_config[name] = st.column_config.TextColumn(
                label=label,
                help=field["description"]
            )

    return column_config

# =========================================================
# DATAFRAME
# =========================================================

def initialize_empty_dataframe(field_defs):

    return pd.DataFrame(
        [{col: None for col in field_defs.keys()}],
        dtype="string"
    )

def load_tsv_into_schema(tsv_file, field_defs):

    schema_cols = list(field_defs.keys())

    df = pd.read_csv(tsv_file, sep="\t", dtype="string")

    return df.reindex(columns=schema_cols).astype("string")

# =========================================================
# FASTA
# =========================================================

def _validate_fasta_headers(fasta_map: dict, max_len: int = 50) -> list[dict]:
    """Return one entry per header whose sequence name exceeds max_len characters."""
    violations = []
    for sample_name, path in fasta_map.items():
        with gzip.open(path, "rt") as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                seq_name = line[1:].split()[0] if line[1:].strip() else ""
                if len(seq_name) > max_len:
                    violations.append({
                        "sample": sample_name,
                        "sequence name": seq_name[:80] + ("…" if len(seq_name) > 80 else ""),
                        "length": len(seq_name),
                    })
    return violations


def _deduplicate_fasta_headers(src_path: Path) -> tuple[Path, bool, list[dict]]:
    """
    Detect and fix duplicate sequence names in a gzipped FASTA.
    ENA's webin-cli considers only the part before the first '|' as the sequence
    name for uniqueness checks, so SAMPLE123|XYUO and SAMPLE123|AAUO collide.
    Appends _N to the prefix segment (before '|') of later duplicates.
    Returns the output path, whether it changed, and the sequence-name rename map.
    """
    seen: set[str] = set()
    needs_dedup = False
    with gzip.open(src_path, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                continue
            raw_name = line[1:].split()[0] if line[1:].strip() else ""
            ena_name = raw_name.split("|")[0]
            if ena_name in seen:
                needs_dedup = True
                break
            seen.add(ena_name)

    if not needs_dedup:
        return src_path, False, []

    tmp = tempfile.NamedTemporaryFile(suffix=".fasta.gz", delete=False)
    tmp.close()

    seen_counts: dict[str, int] = {}
    rename_map = []
    with gzip.open(src_path, "rt") as fin, gzip.open(tmp.name, "wt") as fout:
        for line in fin:
            if line.startswith(">"):
                raw_name = line[1:].split()[0] if line[1:].strip() else ""
                rest = line[1 + len(raw_name):]
                ena_name = raw_name.split("|")[0]
                count = seen_counts.get(ena_name, 0)
                seen_counts[ena_name] = count + 1
                if count > 0:
                    if "|" in raw_name:
                        prefix, suffix = raw_name.split("|", 1)
                        new_name = f"{prefix}_{count}|{suffix}"
                    else:
                        new_name = f"{raw_name}_{count}"
                    fout.write(f">{new_name}{rest}")
                    rename_map.append({"before": raw_name, "after": new_name})
                else:
                    fout.write(line)
            else:
                fout.write(line)

    return Path(tmp.name), True, rename_map


def persist_fastas_temp(uploaded_files):

    fasta_map = {}

    for uf in uploaded_files:

        if not uf.name.endswith(".fasta.gz"):
            continue

        sample_name = uf.name.removesuffix(".fasta.gz")

        if sample_name in fasta_map:
            raise ValueError(f"Duplicate FASTA file: {sample_name}")

        tmp = tempfile.NamedTemporaryFile(
            suffix=".fasta.gz",
            delete=False
        )

        tmp.write(uf.getbuffer())
        tmp.flush()
        tmp.close()

        fasta_map[sample_name] = Path(tmp.name).resolve()

    return fasta_map

# =========================================================
# SUBMISSION WORKER (runs in RQ worker process)
# =========================================================

def _redact_credentials(text, *credentials):
    """Remove credential values from content before it is written to job logs."""
    redacted = "" if text is None else str(text)
    for credential in credentials:
        if credential:
            redacted = redacted.replace(str(credential), "[REDACTED]")
    return redacted


def _submission_task_impl(
    df_records,
    submission,
    fasta_map_str,
    email=None,
    provenance_events=None,
):
    """
    RQ worker function. All arguments must be JSON-serializable.
    df_records: list of dicts from pl.DataFrame.to_dicts()
    fasta_map_str: {sample_name: absolute_path_str}
    """
    job = get_current_job()
    job_id = job.id if job else datetime.now().strftime("%Y%m%d_%H%M%S")

    # RQ must carry credentials while the job is waiting, but the worker already
    # has its local arguments now. Remove credentials from the persisted job data
    # immediately so finished/failed jobs cannot expose them later.
    scrub_job_credentials(job)

    manager.store_start(job_id, TaskStatus.RUNNING)

    log_dir = SUBMISSIONS_DIR / job_id
    log_dir.mkdir(parents=True, exist_ok=True)

    provenance_ledger = list(provenance_events or [])
    _job_provenance_event(
        provenance_ledger,
        log_dir,
        stage="submission",
        action="submission_started",
        entity_type="submission",
        source=submission.get("portal", ""),
        after="running",
        method="MAG2ENA RQ worker",
        method_version=MAG2ENA_VERSION,
        details={"job_id": job_id},
    )

    df = pl.DataFrame(df_records)
    taxonomy_errors = validate_mag_taxonomy_records(
        df_records,
        fetcher=_fetch_ena_taxon,
    )
    if taxonomy_errors:
        _job_provenance_event(
            provenance_ledger,
            log_dir,
            stage="validation",
            action="taxonomy_preflight_failed",
            entity_type="metadata",
            source="ENA taxonomy API",
            after="failed",
            method="MAG taxonomy preflight",
            method_version=MAG2ENA_VERSION,
            decision="submission blocked",
            details={"errors": taxonomy_errors},
        )
        details = "; ".join(
            f"row {error['row']} ({error['sample_name'] or 'unnamed sample'}, "
            f"tax_id {error['tax_id']}): {error['expected']}"
            for error in taxonomy_errors
        )
        raise ValueError(
            "Taxonomy preflight failed before ENA submission: " + details
        )
    _job_provenance_event(
        provenance_ledger,
        log_dir,
        stage="validation",
        action="taxonomy_preflight_passed",
        entity_type="metadata",
        source="ENA taxonomy API",
        after="passed",
        method="MAG taxonomy preflight",
        method_version=MAG2ENA_VERSION,
        decision="continue",
        details={
            "sample_count": len(df_records),
            "unique_tax_id_count": len(
                {str(row.get("tax_id", "")).strip() for row in df_records}
            ),
            "metagenome_taxa_allowed": False,
            "required_rank": "species",
        },
    )

    fasta_map = {k: Path(v) for k, v in fasta_map_str.items()}
    fasta_inspections = inspect_fasta_collection(fasta_map)
    confirmed_single_contigs = require_single_contig_confirmation(
        fasta_inspections,
        bool(submission.get("single_contig_chromosomes_confirmed")),
    )

    samples_submitted = 0
    samples_error = 0

    RESERVED_COLUMNS = {
        "sample_name",
        "organism",
        "tax_id",
        "ENA-CHECKLIST",
    }

    df = df.with_columns(
        pl.lit("ERC000047").alias("ENA-CHECKLIST")
    )

    project_accession = submission["study_accession"]

    root = ET.Element("WEBIN")

    submission_set = ET.fromstring("""
    <SUBMISSION_SET>
        <SUBMISSION>
            <ACTIONS>
                <ACTION>
                    <ADD/>
                </ACTION>
            </ACTIONS>
        </SUBMISSION>
    </SUBMISSION_SET>
    """)

    hold_until_date = submission.get("hold_until_date")
    if hold_until_date:
        actions = submission_set.find(".//ACTIONS")
        hold_action = ET.SubElement(actions, "ACTION")
        ET.SubElement(hold_action, "HOLD", {"HoldUntilDate": hold_until_date})

    root.append(submission_set)

    if not project_accession:

        project_set = ET.fromstring(f"""
        <PROJECT_SET>
            <PROJECT alias="{submission["study_name"]}">
            <TITLE>{submission["study_title"]}</TITLE>
            <DESCRIPTION>{submission["study_description"]}</DESCRIPTION>
            <SUBMISSION_PROJECT>
                <SEQUENCING_PROJECT/>
            </SUBMISSION_PROJECT>
            </PROJECT>
        </PROJECT_SET>
        """)

        root.append(project_set)

    sample_set = ET.SubElement(root, "SAMPLE_SET")

    for row in df.iter_rows(named=True):

        sample = ET.SubElement(sample_set, "SAMPLE", {
            "alias": str(row["sample_name"]),
            "center_name": ""
        })

        ET.SubElement(sample, "TITLE").text = (
            "This sample represents a MAG derived from the metagenomic sample "
            + str(row["sample derived from"])
        )

        sname = ET.SubElement(sample, "SAMPLE_NAME")

        ET.SubElement(sname, "TAXON_ID").text = str(row["tax_id"])
        ET.SubElement(sname, "SCIENTIFIC_NAME").text = str(row["organism"])

        attrs = ET.SubElement(sample, "SAMPLE_ATTRIBUTES")

        def add_attr(tag, value, units=None):

            sa = ET.SubElement(attrs, "SAMPLE_ATTRIBUTE")

            ET.SubElement(sa, "TAG").text = tag
            ET.SubElement(sa, "VALUE").text = str(value)

            if units:
                ET.SubElement(sa, "UNITS").text = units

        explicit_cols = {
            "metagenomic source",
            "sample derived from",
            "project name",
            "completeness score",
            "completeness software",
            "contamination score",
            "binning software",
            "assembly quality",
            "binning parameters",
            "taxonomic identity marker",
            "isolation_source",
            "collection date",
            "geographic location (latitude)",
            "geographic location (longitude)",
            "broad-scale environmental context",
            "local environmental context",
            "environmental medium",
            "geographic location (country and/or sea)",
            "assembly software",
            "platform",
            "genome coverage",
        }

        add_attr("metagenomic source", row["metagenomic source"])
        add_attr("sample derived from", row["sample derived from"])
        add_attr("project name", row["project name"])
        add_attr("completeness score", row["completeness score"], "%")
        add_attr("completeness software", row["completeness software"])
        add_attr("contamination score", row["contamination score"], "%")
        add_attr("binning software", row["binning software"])
        add_attr("assembly quality", row["assembly quality"])
        add_attr("binning parameters", row["binning parameters"])
        add_attr("taxonomic identity marker", row["taxonomic identity marker"])
        add_attr("isolation_source", row["isolation_source"])
        add_attr("collection date", row["collection date"])

        add_attr(
            "geographic location (latitude)",
            row["geographic location (latitude)"],
            "DD"
        )

        add_attr(
            "geographic location (longitude)",
            row["geographic location (longitude)"],
            "DD"
        )

        add_attr(
            "broad-scale environmental context",
            row["broad-scale environmental context"]
        )

        add_attr(
            "local environmental context",
            row["local environmental context"]
        )

        add_attr(
            "environmental medium",
            row["environmental medium"]
        )

        add_attr(
            "geographic location (country and/or sea)",
            row["geographic location (country and/or sea)"]
        )

        add_attr(
            "assembly software",
            row["assembly software"]
        )

        for col, value in row.items():

            if col in RESERVED_COLUMNS:
                continue

            if col in explicit_cols:
                continue

            if value is None or str(value).strip() == "":
                continue

            add_attr(col, value)

        add_attr("ENA-CHECKLIST", row["ENA-CHECKLIST"])

    xml_bytes = ET.tostring(
        root,
        encoding="UTF-8",
        xml_declaration=True
    )

    with open(log_dir / "submit.xml", "wb") as f:
        f.write(xml_bytes)

    _job_provenance_event(
        provenance_ledger,
        log_dir,
        stage="submission",
        action="submission_xml_generated",
        entity_type="file",
        source="validated metadata",
        after=sha256_bytes(xml_bytes),
        method="WEBIN XML generation",
        method_version=MAG2ENA_VERSION,
        details={
            "filename": "submit.xml",
            "sha256": sha256_bytes(xml_bytes),
            "size_bytes": len(xml_bytes),
            "sample_count": len(df),
            "checklist": "ERC000047",
        },
    )

    if submission["portal"] == "Testing":
        url = "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/submit/queue"
    else:
        url = "https://www.ebi.ac.uk/ena/submit/webin-v2/submit/queue"

    auth = (
        submission["ena_user"],
        submission["ena_password"]
    )

    headers = {
        "Accept": "application/json"
    }

    files = {
        "file": ("submit.xml", xml_bytes, "text/xml")
    }

    response = requests.post(
        url,
        headers=headers,
        files=files,
        auth=auth
    )

    response.raise_for_status()

    response_json = response.json()

    poll_url = response_json["_links"]["poll"]["href"]

    response = requests.get(poll_url, auth=auth)

    while response.status_code != 200:
        time.sleep(5)
        response = requests.get(poll_url, auth=auth)

    xml_text = _redact_credentials(response.text, *auth)

    xml_log = log_dir / "webin_log.xml"

    with open(xml_log, "w", encoding="utf-8") as f:
        f.write(xml_text)

    root = ET.fromstring(xml_text)

    alias_to_accession = {}

    for sample in root.findall(".//SAMPLE"):

        accession = sample.get("accession")

        if accession:
            alias_to_accession[
                sample.get("alias")
            ] = accession

    error_regex = re.compile(
        r'alias: "(.+?)".*accession: "(.+?)"'
    )

    for err in root.findall(".//ERROR"):

        error_msg = err.text

        match = error_regex.search(error_msg)

        if match:

            alias, accession = match.groups()

            alias_to_accession[alias] = accession

    if not project_accession:

        for project in root.findall(".//PROJECT"):
            project_accession = project.get("accession")

    receipt_success = str(root.get("success", "")).lower() == "true"
    submission_accessions = [
        item.get("accession", "")
        for item in root.findall(".//SUBMISSION")
        if item.get("accession")
    ]
    _job_provenance_event(
        provenance_ledger,
        log_dir,
        stage="submission",
        action="ena_metadata_receipt_received",
        entity_type="receipt",
        source=url,
        after="success" if receipt_success else "failure",
        method="ENA Webin REST asynchronous submission",
        decision="accepted" if receipt_success else "rejected",
        details={
            "job_id": job_id,
            "portal": submission.get("portal", ""),
            "poll_url": poll_url,
            "project_accession": project_accession or "",
            "submission_accessions": submission_accessions,
            "sample_accessions": alias_to_accession,
            "receipt_sha256": sha256_bytes(xml_text.encode("utf-8")),
        },
    )

    manifests = {}

    for alias, sample_accession in alias_to_accession.items():

        if sample_accession is None:
            samples_error += 1
            _job_provenance_event(
                provenance_ledger,
                log_dir,
                stage="submission",
                action="assembly_submission_skipped",
                entity_type="assembly",
                sample_name=alias,
                source="ENA metadata receipt",
                after="no sample accession",
                method="MAG2ENA submission worker",
                method_version=MAG2ENA_VERSION,
                decision="skipped",
            )
            with open(log_dir / "error.txt", "a") as f:
                f.write(
                    f"SAMPLE : {alias}\n"
                    "No sample accession returned — FASTA not submitted.\n"
                )
            continue

        orig_fasta_path = fasta_map[alias]
        dedup_fasta_path = None

        inspection = fasta_inspections[alias]
        seq_count = inspection["sequence_count"]

        fasta_path = orig_fasta_path
        rename_map = []
        if seq_count > 1:
            fasta_path, modified, rename_map = _deduplicate_fasta_headers(orig_fasta_path)
            if modified:
                dedup_fasta_path = fasta_path

        original_fasta_sha256 = sha256_path(orig_fasta_path)
        prepared_fasta_sha256 = sha256_path(fasta_path)
        _job_provenance_event(
            provenance_ledger,
            log_dir,
            stage="fasta",
            action="fasta_prepared",
            entity_type="file",
            sample_name=alias,
            field="FASTA",
            source=f"{alias}.fasta.gz",
            before=original_fasta_sha256,
            after=prepared_fasta_sha256,
            method="ENA FASTA preprocessor",
            method_version=MAG2ENA_VERSION,
            details={
                "sequence_count": seq_count,
                "duplicate_headers_renamed": len(rename_map),
                "rename_map": rename_map,
                "output_changed": bool(dedup_fasta_path),
            },
        )

        sample_row = df.filter(pl.col("sample_name") == alias)

        base_manifest = f"""STUDY   {project_accession}
SAMPLE   {sample_accession}
ASSEMBLYNAME   {alias}
ASSEMBLY_TYPE   Metagenome-Assembled Genome (MAG)
COVERAGE   {sample_row["genome coverage"].item()}
PROGRAM   {sample_row["assembly software"].item()}
PLATFORM   {sample_row["platform"].item()}
FASTA   {fasta_path}"""

        chromosome_gz_path = None

        if seq_count == 1:
            chromosome_list = chromosome_list_text(inspection)

            with tempfile.NamedTemporaryFile(
                mode="wb+",
                suffix=".gz",
                delete=False
            ) as tmp_chromosome_gz:

                with gzip.GzipFile(
                    fileobj=tmp_chromosome_gz,
                    mode="wb"
                ) as gz_file:

                    gz_file.write(
                        chromosome_list.encode("utf-8")
                    )

                chromosome_gz_path = tmp_chromosome_gz.name

            base_manifest += (
                f"\nCHROMOSOME_LIST   {chromosome_gz_path}"
            )
            _job_provenance_event(
                provenance_ledger,
                log_dir,
                stage="fasta",
                action="chromosome_list_generated",
                entity_type="file",
                sample_name=alias,
                field="CHROMOSOME_LIST",
                source=f"{alias}.fasta.gz",
                after=sha256_path(chromosome_gz_path),
                method="single-contig chromosome-list generation",
                method_version=MAG2ENA_VERSION,
                details={
                    "sequence_name": inspection["first_sequence_name"],
                    "chromosome_name": "1",
                    "role": "chromosome",
                    "user_confirmed_complete_chromosome": (
                        alias in confirmed_single_contigs
                    ),
                },
            )

        manifests[alias] = base_manifest
        manifest_sha256 = sha256_bytes(base_manifest.encode("utf-8"))
        _job_provenance_event(
            provenance_ledger,
            log_dir,
            stage="submission",
            action="webin_manifest_generated",
            entity_type="file",
            sample_name=alias,
            field="manifest",
            source="validated metadata and prepared FASTA",
            after=manifest_sha256,
            method="Webin-CLI manifest generation",
            method_version=MAG2ENA_VERSION,
            details={
                "study_accession": project_accession or "",
                "sample_accession": sample_accession,
                "sha256": manifest_sha256,
                "chromosome_list_included": bool(chromosome_gz_path),
            },
        )

        with tempfile.NamedTemporaryFile(
            mode="w+",
            delete=False
        ) as tmp_manifest:

            tmp_manifest.write(base_manifest)

            manifest_path = tmp_manifest.name

        try:

            password_env_name = "MAG2ENA_WEBIN_PASSWORD"
            webin_env = os.environ.copy()
            webin_env[password_env_name] = auth[1]

            cmd = [
                "java",
                "-jar",
                WEBIN_JAR,
                "-username",
                auth[0],
                f"-passwordEnv={password_env_name}",
                "-context",
                "genome",
                "-manifest",
                manifest_path,
                "-submit",
            ]

            if submission["portal"] == "Testing":
                cmd.append("-test")

            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    env=webin_env,
                )
            finally:
                webin_env.pop(password_env_name, None)

            safe_stdout = _redact_credentials(result.stdout, *auth)

            success = "successfully" in safe_stdout.lower()

            logfile = (
                log_dir / "success.txt"
                if success
                else log_dir / "error.txt"
            )

            with open(logfile, "a") as f:
                f.write(
                    f"SAMPLE : {alias}\n"
                    f"{safe_stdout}\n"
                )

            if success:
                samples_submitted += 1
            else:
                samples_error += 1

            _job_provenance_event(
                provenance_ledger,
                log_dir,
                stage="submission",
                action="assembly_submitted",
                entity_type="assembly",
                sample_name=alias,
                source="ENA Webin-CLI",
                after="success" if success else "failure",
                method="Webin-CLI genome submission",
                method_version=Path(WEBIN_JAR).stem,
                decision="accepted" if success else "rejected",
                details={
                    "portal": submission.get("portal", ""),
                    "project_accession": project_accession or "",
                    "sample_accession": sample_accession,
                    "manifest_sha256": manifest_sha256,
                    "prepared_fasta_sha256": prepared_fasta_sha256,
                    "return_code": result.returncode,
                },
            )

        finally:

            Path(manifest_path).unlink(missing_ok=True)

            if chromosome_gz_path:
                Path(chromosome_gz_path).unlink(missing_ok=True)

            if dedup_fasta_path:
                dedup_fasta_path.unlink(missing_ok=True)

    with zipfile.ZipFile(log_dir / "manifests.zip", "w", zipfile.ZIP_DEFLATED) as zf:
        for alias, content in manifests.items():
            zf.writestr(f"{alias}.manifest", content)

    manifests_path = log_dir / "manifests.zip"
    _job_provenance_event(
        provenance_ledger,
        log_dir,
        stage="submission",
        action="manifest_archive_generated",
        entity_type="file",
        source="per-MAG Webin-CLI manifests",
        after=sha256_path(manifests_path),
        method="ZIP archive generation",
        method_version=MAG2ENA_VERSION,
        details={
            "filename": manifests_path.name,
            "sha256": sha256_path(manifests_path),
            "manifest_count": len(manifests),
        },
    )

    # Clean up temp FASTA files
    for path in fasta_map.values():
        try:
            path.unlink(missing_ok=True)
        except Exception:
            pass

    task_result = {
        "submitted": samples_submitted,
        "errors": samples_error,
    }

    with open(log_dir / "result.json", "w") as f:
        json.dump(task_result, f)

    result_path = log_dir / "result.json"
    _job_provenance_event(
        provenance_ledger,
        log_dir,
        stage="submission",
        action="submission_completed",
        entity_type="submission",
        source=submission.get("portal", ""),
        after="completed",
        method="MAG2ENA RQ worker",
        method_version=MAG2ENA_VERSION,
        decision="completed with errors" if samples_error else "completed",
        details={
            "job_id": job_id,
            "submitted": samples_submitted,
            "errors": samples_error,
            "result_sha256": sha256_path(result_path),
        },
    )

    return task_result


def submission_task(
    df_records,
    submission,
    fasta_map_str,
    email=None,
    encryption_password=None,
    provenance_events=None,
):
    """Run an ENA job and always encrypt any resulting artifacts when requested."""
    job = get_current_job()
    job_id = job.id if job else datetime.now().strftime("%Y%m%d_%H%M%S")
    scrub_job_credentials(job)

    try:
        return _submission_task_impl(
            df_records,
            submission,
            fasta_map_str,
            email=email,
            provenance_events=provenance_events,
        )
    except Exception as exc:
        # Keep the latest worker ledger and add a terminal failure event without
        # recording credentials or submission payloads.
        log_dir = SUBMISSIONS_DIR / job_id
        log_dir.mkdir(parents=True, exist_ok=True)
        ledger = list(provenance_events or [])
        json_path = log_dir / "provenance.json"
        try:
            if json_path.exists():
                loaded = json.loads(json_path.read_text(encoding="utf-8"))
                if isinstance(loaded, list):
                    ledger = loaded
        except Exception:
            pass

        safe_error = _redact_credentials(
            str(exc),
            submission.get("ena_user"),
            submission.get("ena_password"),
        )
        try:
            _job_provenance_event(
                ledger,
                log_dir,
                stage="submission",
                action="submission_failed",
                entity_type="submission",
                source=submission.get("portal", ""),
                after="failed",
                method="MAG2ENA RQ worker",
                method_version=MAG2ENA_VERSION,
                decision="failed",
                details={
                    "job_id": job_id,
                    "error_type": type(exc).__name__,
                    "error": safe_error,
                },
            )
        except Exception:
            # Provenance must never replace the original submission exception.
            pass
        raise
    finally:
        if encryption_password:
            log_dir = SUBMISSIONS_DIR / job_id
            if log_dir.is_dir() and any(log_dir.iterdir()):
                encrypt_job_folder(
                    log_dir,
                    encryption_password,
                    allowed_root=SUBMISSIONS_DIR,
                )

def _step_progress_html(current_step: int) -> str:
    steps = [(1, "Metadata"), (2, "FASTA Files"), (3, "Submit to ENA")]
    parts = []
    for i, label in steps:
        if i < current_step:
            cls, icon = "sp-done", "✓"
        elif i == current_step:
            cls, icon = "sp-active", str(i)
        else:
            cls, icon = "sp-pending", str(i)
        parts.append(
            f'<div class="sp-item {cls}">'
            f'<span class="sp-num">{icon}</span>'
            f'<span class="sp-label">{label}</span>'
            f'</div>'
        )
        if i < 3:
            parts.append('<div class="sp-line"></div>')
    return f'<div class="step-progress">{"".join(parts)}</div>'


def _metadata_fingerprint(df: pd.DataFrame) -> str:
    """Return a stable fingerprint used to reject stale resolver proposals."""
    normalised = df.fillna("").astype(str)
    hashed = pd.util.hash_pandas_object(normalised, index=True)
    return str(int(hashed.sum()))


def _ensure_provenance_state() -> list[dict]:
    if "provenance_events" not in st.session_state:
        st.session_state.provenance_events = []
    return st.session_state.provenance_events


def _append_provenance(events) -> None:
    events = list(events)
    if not events:
        return
    ledger = _ensure_provenance_state()
    ledger.extend(events)
    st.session_state.provenance_events = ledger


def _record_metadata_changes(
    before: pd.DataFrame,
    after: pd.DataFrame,
    *,
    action: str,
    source: str,
    method: str,
    method_version: str = MAG2ENA_VERSION,
    decision: str = "",
    details: dict | None = None,
) -> list[dict]:
    events = dataframe_change_events(
        before,
        after,
        stage="metadata",
        action=action,
        source=source,
        method=method,
        method_version=method_version,
        decision=decision,
        details=details,
    )
    _append_provenance(events)
    return events


def _reset_provenance(*, action: str = "workflow_started", source: str = "web app"):
    st.session_state.provenance_events = [
        make_event(
            stage="workflow",
            action=action,
            entity_type="workflow",
            source=source,
            method="MAG2ENA web workflow",
            method_version=MAG2ENA_VERSION,
        )
    ]
    st.session_state.pop("_provenance_fasta_hashes", None)
    st.session_state.pop("_provenance_validated_fasta_key", None)


def _record_validation(df: pd.DataFrame, status: str, errors=None) -> None:
    if isinstance(errors, pd.DataFrame):
        error_records = errors.fillna("").to_dict(orient="records")
    else:
        error_records = list(errors or [])
    _append_provenance(
        [
            make_event(
                stage="validation",
                action="metadata_validated",
                entity_type="metadata_table",
                source="ERC000047",
                after=status,
                method="checklist-derived pre-submission validation",
                method_version=MAG2ENA_VERSION,
                decision=status,
                details={
                    "checklist": "ERC000047",
                    "checklist_sha256": sha256_path(CHECKLIST_XML),
                    "metadata_sha256": dataframe_sha256(df),
                    "rows": len(df),
                    "error_count": len(error_records),
                    "errors": error_records,
                },
            )
        ]
    )


def _record_fasta_inputs(fasta_map: dict, *, source: str) -> None:
    known_hashes = st.session_state.get("_provenance_fasta_hashes", {})
    events = []
    for sample_name, path in fasta_map.items():
        digest = sha256_path(path)
        if known_hashes.get(sample_name) == digest:
            continue
        known_hashes[sample_name] = digest
        events.append(
            make_event(
                stage="fasta",
                action="fasta_uploaded",
                entity_type="file",
                sample_name=sample_name,
                field="FASTA",
                source=source,
                after=digest,
                method="gzip-compressed FASTA upload",
                method_version=MAG2ENA_VERSION,
                details={
                    "filename": f"{sample_name}.fasta.gz",
                    "size_bytes": Path(path).stat().st_size,
                    "sha256": digest,
                },
            )
        )
    st.session_state._provenance_fasta_hashes = known_hashes
    _append_provenance(events)


def _job_provenance_event(
    ledger: list[dict],
    log_dir: Path,
    **event_kwargs,
) -> dict:
    event = make_event(**event_kwargs)
    ledger.append(event)
    write_provenance_files(ledger, log_dir)
    return event


def _environment_proposal_event(
    proposal: dict,
    *,
    action: str,
    decision: str,
    before=None,
) -> dict:
    candidates = proposal.get("candidates", [])
    return make_event(
        stage="metadata",
        action=action,
        entity_type="ontology_proposal",
        field=proposal.get("field", ""),
        source="OLS4 / ENVO",
        before=(
            proposal.get("original_value", "") if before is None else before
        ),
        after=proposal.get("selected_value", ""),
        method="field-aware environmental resolver",
        method_version=proposal.get("ontology_version", MAG2ENA_VERSION),
        decision=decision,
        details={
            "proposal_id": proposal.get("proposal_id", ""),
            "sample_names": proposal.get("sample_names", []),
            "row_indices": proposal.get("row_indices", []),
            "parent_accession": proposal.get("parent_accession", ""),
            "source_evidence": proposal.get("evidence", ""),
            "selected_label": proposal.get("selected_label", ""),
            "selected_id": proposal.get("selected_id", ""),
            "confidence": proposal.get("confidence", ""),
            "score": proposal.get("score", 0.0),
            "reason": proposal.get("reason", ""),
            "alternatives": [
                candidate.get("value", "") for candidate in candidates[1:6]
            ],
        },
    )


def _clear_environment_resolver_state(*, clear_history: bool = False):
    for key in (
        "envo_proposals",
        "envo_review_df",
        "envo_source_fingerprint",
        "envo_resolver_scope",
    ):
        st.session_state.pop(key, None)
    st.session_state.envo_review_key = st.session_state.get("envo_review_key", 0) + 1
    if clear_history:
        st.session_state.pop("envo_undo_df", None)


@st.dialog("Job submitted")
def job_submitted_dialog(job_id, encrypted=False):
    st.session_state["last_job_id"] = job_id
    st.success("Job submitted to the queue successfully.")
    st.markdown("Copy your **Job ID** and paste it into the **Jobs** page to track progress:")
    st.code(job_id, language=None)
    st.caption("This ID has been saved — it will be pre-filled when you visit the Jobs page.")
    if encrypted:
        st.warning(
            "This job's results will be encrypted. Save the results password securely: "
            "MAG2ENA does not store it and cannot recover or reset it."
        )

# =========================================================
# UI
# =========================================================

def runUI():

    if not Path(WEBIN_JAR).exists():
        st.error(f"Missing Webin CLI jar: {WEBIN_JAR}")
        st.stop()

    field_defs = load_fields_from_xml(CHECKLIST_XML)

    if "metadata_df" not in st.session_state:
        st.session_state.metadata_df = initialize_empty_dataframe(field_defs)

    if "editor_key" not in st.session_state:
        st.session_state.editor_key = 0

    if not _ensure_provenance_state():
        _reset_provenance()

    _step_placeholder = st.empty()
    _initial_step = 2 if "validated_df" in st.session_state else 1
    _step_placeholder.markdown(_step_progress_html(_initial_step), unsafe_allow_html=True)

    st.info(
        "This page walks you through submitting Metagenome-Assembled Genomes (MAGs) to ENA. "
        "Fill in your metadata, upload one `.fasta.gz` file per MAG, enter your ENA Webin "
        "credentials, and submit — all XML and Webin-CLI manifests are generated automatically. "
        "**Limit: 1,000 MAGs per web job.** For larger collections, use the included MAG2ENA "
        "CLI, which batches metadata registration automatically and submits assemblies "
        "serially.",
        icon="ℹ️",
    )

    # ── Step 1: Metadata ──────────────────────────────────────────────────────
    st.subheader("Step 1 — Metadata")

    use_example = st.toggle(
        "Load example data",
        key="use_example_toggle",
        help=(
            "Populate the metadata table with two example MAGs and automatically load their "
            "FASTA files — lets you test the full submission workflow without your own data."
        ),
    )

    if use_example and not st.session_state.get("_example_active"):
        example_df = load_tsv_into_schema(
            EXAMPLES_DIR / "metadata.tsv", field_defs
        )
        _reset_provenance(action="example_workflow_loaded", source="bundled example")
        _record_metadata_changes(
            initialize_empty_dataframe(field_defs),
            example_df,
            action="metadata_imported",
            source="examples/metadata.tsv",
            method="TSV schema import",
            details={
                "sha256": sha256_path(EXAMPLES_DIR / "metadata.tsv"),
                "rows": len(example_df),
            },
        )
        st.session_state.metadata_df = example_df
        st.session_state.pop("validated_df", None)
        st.session_state.pop("tsv_filename", None)
        st.session_state.pop("_example_fasta_map", None)
        _clear_environment_resolver_state(clear_history=True)
        st.session_state.editor_key += 1
        st.session_state._example_active = True
        st.rerun()
    elif not use_example and st.session_state.get("_example_active"):
        _reset_provenance(action="workflow_reset", source="example toggle")
        st.session_state.metadata_df = initialize_empty_dataframe(field_defs)
        st.session_state.pop("validated_df", None)
        st.session_state.pop("_example_fasta_map", None)
        _clear_environment_resolver_state(clear_history=True)
        st.session_state.editor_key += 1
        st.session_state._example_active = False
        st.rerun()

    col_upload, col_template = st.columns([3, 1])

    with col_upload:
        uploaded_tsv = st.file_uploader(
            "Upload metadata TSV (optional)",
            type=["tsv"],
            help="Upload an existing TSV to pre-fill the table. You can continue editing it below.",
            disabled=use_example,
        )

    with col_template:
        template_df = initialize_empty_dataframe(field_defs)
        tsv_bytes = template_df.to_csv(sep="\t", index=False).encode()
        st.download_button(
            "Download template",
            data=tsv_bytes,
            file_name="metadata_template.tsv",
            mime="text/tab-separated-values",
            width="stretch",
            help="Download an empty TSV template with the required columns."
        )

    # Only reload from TSV when a new file is uploaded (different filename),
    # so that subsequent user edits to the table are not overwritten on re-render.
    if not use_example and uploaded_tsv:
        if st.session_state.get("tsv_filename") != uploaded_tsv.name:
            before_import = st.session_state.metadata_df.copy()
            imported_df = load_tsv_into_schema(uploaded_tsv, field_defs)
            file_bytes = uploaded_tsv.getvalue()
            _append_provenance(
                [
                    make_event(
                        stage="metadata",
                        action="metadata_file_uploaded",
                        entity_type="file",
                        source=uploaded_tsv.name,
                        after=sha256_bytes(file_bytes),
                        method="TSV upload",
                        method_version=MAG2ENA_VERSION,
                        details={
                            "filename": uploaded_tsv.name,
                            "size_bytes": len(file_bytes),
                            "sha256": sha256_bytes(file_bytes),
                            "rows": len(imported_df),
                        },
                    )
                ]
            )
            _record_metadata_changes(
                before_import,
                imported_df,
                action="metadata_imported",
                source=uploaded_tsv.name,
                method="TSV schema import",
                details={"sha256": sha256_bytes(file_bytes)},
            )
            st.session_state.metadata_df = imported_df
            st.session_state.tsv_filename = uploaded_tsv.name
            st.session_state.pop("validated_df", None)
            _clear_environment_resolver_state(clear_history=True)
            st.session_state.editor_key += 1

    row_count = len(st.session_state.metadata_df)
    st.caption(
        f"{row_count} row{'s' if row_count != 1 else ''} · "
        "Columns marked with * are mandatory · "
        "Hover column headers for field descriptions"
    )

    previous_metadata_df = st.session_state.metadata_df.copy()
    edited_df = st.data_editor(
        st.session_state.metadata_df,
        num_rows="dynamic",
        column_config=build_column_config(field_defs),
        hide_index=True,
        width="stretch",
        key=f"metadata_editor_{st.session_state.editor_key}",
    )

    # Persist edits back to session state so re-renders don't discard them.
    if not edited_df.equals(previous_metadata_df):
        _record_metadata_changes(
            previous_metadata_df,
            edited_df,
            action="manual_edit",
            source="metadata table",
            method="Streamlit data editor",
            decision="user entered",
        )
        st.session_state.pop("validated_df", None)
        _clear_environment_resolver_state(clear_history=True)
    st.session_state.metadata_df = edited_df

    render_geographic_overview(
        st.session_state.metadata_df,
        key="home_geographic_overview",
    )

    # ── Step 1b: Metadata assistance ─────────────────────────────────────────
    with st.expander("Metadata assistance"):

        tab_ena, tab_tax, tab_envo, tab_fill, tab_import = st.tabs([
            "ENA autofill",
            "Taxonomy resolver",
            "Environmental resolver",
            "Fill column",
            "Import quality files",
        ])

        # ── ENA autofill ──────────────────────────────────────────────────────
        with tab_ena:
            st.caption(
                "Look up the original metagenome in ENA using the accession(s) in "
                "`sample derived from` and autofill environmental metadata. "
                "Fills: `metagenomic source`, `isolation_source`, `collection date`, "
                "`geographic location (latitude/longitude)`, "
                "`broad-scale environmental context`, `local environmental context`, "
                "`environmental medium`, `geographic location (country and/or sea)`. "
                "Only empty cells are overwritten — existing values are preserved."
            )

            derived_col = st.session_state.metadata_df.get(
                "sample derived from", pd.Series(dtype="string")
            )
            unique_accs = sorted({
                str(a).strip()
                for a in derived_col.dropna()
                if str(a).strip() not in ("", "<NA>", "None", "nan")
            })

            if unique_accs:
                st.info(
                    f"Found {len(unique_accs)} unique accession(s) in `sample derived from`: "
                    + ", ".join(f"`{a}`" for a in unique_accs[:6])
                    + ("…" if len(unique_accs) > 6 else "")
                )
            else:
                st.warning(
                    "No accessions found in `sample derived from`. "
                    "Fill in that column first, then come back here."
                )

            col_prev_input, col_prev_btn = st.columns([3, 1])
            with col_prev_input:
                preview_acc = st.text_input(
                    "Preview a single accession",
                    placeholder="e.g. ERS1234567 or SRS1981904",
                    key="ena_preview_acc",
                    help=(
                        "Enter an ENA sample accession to inspect which environmental "
                        "metadata MAG2ENA can reuse before applying it to the table."
                    ),
                )
            with col_prev_btn:
                st.write("")
                preview_clicked = st.button(
                    "Preview",
                    key="btn_ena_preview",
                    width="stretch",
                )

            if preview_clicked:
                acc = preview_acc.strip()
                if not acc:
                    st.warning("Enter an accession to preview.")
                else:
                    with st.spinner(f"Querying ENA for {acc}…"):
                        raw = _fetch_ena_sample_metadata(acc)
                    if "_error" in raw:
                        st.error(f"Could not fetch metadata: {raw['_error']}")
                        st.session_state.pop("ena_preview_result", None)
                        st.session_state.pop("ena_preview_source", None)
                    else:
                        st.session_state.ena_preview_result = _attrs_to_mag_metadata(raw)
                        st.session_state.ena_preview_source = {
                            "source": raw.get("_metadata_source", ""),
                            "source_url": raw.get("_metadata_source_url", ""),
                            "assemblies": raw.get("_fallback_assembly_accessions", []),
                            "conflicts": raw.get("_fallback_conflicts", {}),
                        }

            preview_result = st.session_state.get("ena_preview_result")
            if preview_result:
                preview_source = st.session_state.get("ena_preview_source", {})
                if preview_source.get("source"):
                    st.caption(f"Metadata source: {preview_source['source']}")
                if preview_source.get("assemblies"):
                    st.caption(
                        "Linked assembly record(s): "
                        + ", ".join(preview_source["assemblies"])
                    )
                if preview_source.get("conflicts"):
                    st.warning(
                        "Conflicting values across linked assemblies were not used: "
                        + ", ".join(sorted(preview_source["conflicts"]))
                    )
                preview_rows = [
                    {"field": k, "value": v}
                    for k, v in preview_result.items()
                    if v
                ]
                if preview_rows:
                    st.dataframe(
                        pd.DataFrame(preview_rows),
                        hide_index=True,
                        width="stretch",
                    )
                else:
                    st.info("No fillable attributes found for this accession.")

            st.divider()

            if st.button(
                "Autofill from all accessions",
                key="btn_ena_autofill",
                type="primary",
                width="stretch",
                disabled=not unique_accs,
            ):
                df_mod = st.session_state.metadata_df.copy()
                cells_filled = 0
                failed_accs = []
                metadata_by_accession = {}
                fallback_accs = []

                with st.spinner(
                    f"Fetching metadata for {len(unique_accs)} accession(s)…"
                ):
                    for acc in unique_accs:
                        raw = _fetch_ena_sample_metadata(acc)
                        if "_error" in raw:
                            failed_accs.append(acc)
                            continue
                        metadata_by_accession[acc] = raw
                        if raw.get("_fallback_assembly_accessions"):
                            fallback_accs.append(acc)
                        mapped = _attrs_to_mag_metadata(raw)
                        mask = df_mod["sample derived from"] == acc
                        for col, val in mapped.items():
                            if col not in df_mod.columns or not val:
                                continue
                            for idx in df_mod.index[mask]:
                                current = df_mod.at[idx, col]
                                if pd.isna(current) or str(current).strip() in (
                                    "", "<NA>", "None", "nan"
                                ):
                                    df_mod.at[idx, col] = val
                                    cells_filled += 1

                autofill_events = _record_metadata_changes(
                    st.session_state.metadata_df,
                    df_mod,
                    action="ena_parent_autofill",
                    source="ENA Browser and Portal APIs",
                    method="parent-record mapping and normalization",
                    details={
                        "only_empty_cells": True,
                        "assembly_fallback_accessions": fallback_accs,
                    },
                )
                for event in autofill_events:
                    row_index = event.get("details", {}).get("row_index")
                    try:
                        parent_accession = str(
                            df_mod.at[int(row_index), "sample derived from"]
                        ).strip()
                    except (KeyError, TypeError, ValueError):
                        parent_accession = ""
                    raw_source = metadata_by_accession.get(parent_accession, {})
                    if raw_source.get("_metadata_source"):
                        event["source"] = raw_source["_metadata_source"]
                    event["details"].update(
                        {
                            "parent_accession": parent_accession,
                            "source_url": raw_source.get("_metadata_source_url", ""),
                            "fallback_assembly_accessions": raw_source.get(
                                "_fallback_assembly_accessions", []
                            ),
                            "omitted_fallback_conflicts": raw_source.get(
                                "_fallback_conflicts", {}
                            ),
                        }
                    )
                st.session_state.metadata_df = df_mod
                st.session_state.pop("validated_df", None)
                st.session_state.pop("ena_preview_result", None)
                st.session_state.pop("ena_preview_source", None)
                _clear_environment_resolver_state(clear_history=True)
                st.session_state.editor_key += 1

                if cells_filled:
                    st.toast(f"Filled {cells_filled} cell(s) from ENA.", icon="✅")
                if failed_accs:
                    st.warning(
                        "Could not fetch metadata for: "
                        + ", ".join(f"`{a}`" for a in failed_accs)
                    )
                if fallback_accs:
                    st.toast(
                        "Recovered parent metadata through exact linked assemblies for "
                        f"{len(fallback_accs)} accession(s).",
                        icon="🔗",
                    )
                if not cells_filled and not failed_accs:
                    st.info("No empty cells to fill — all matching fields already have values.")
                st.rerun()

        # ── Taxonomy resolver ─────────────────────────────────────────────────
        with tab_tax:
            st.caption(
                "Search the ENA taxonomy database by organism name to find the correct "
                "scientific name and NCBI taxon ID."
            )
            st.info(
                "**MAG taxonomy policy:** each `tax_id` must resolve to exactly species "
                "rank. Genus or higher-rank assignments, subspecies/strain records, and "
                "metagenome taxonomies are blocked before submission. If your assignment "
                "is below species rank, use its species-level parent taxon."
            )
            tax_query = st.text_input(
                "Organism name",
                key="tax_query",
                placeholder="e.g. uncultured Firmicutes bacterium",
                help="Search by a scientific organism name to find its ENA/NCBI taxon ID.",
            )
            col_search, col_batch = st.columns(2)

            with col_search:
                if st.button("Search taxonomy", key="btn_tax_search", width="stretch"):
                    if tax_query:
                        with st.spinner("Querying ENA taxonomy…"):
                            st.session_state.tax_results = _ena_taxonomy_search(tax_query)
                    else:
                        st.warning("Enter an organism name first.")

            with col_batch:
                if st.button(
                    "Batch resolve tax IDs",
                    key="btn_batch_tax",
                    width="stretch",
                    help=(
                        "For every row with `organism` filled but `tax_id` empty, "
                        "query the ENA taxonomy API and fill the best match."
                    ),
                ):
                    df_mod = st.session_state.metadata_df.copy()
                    resolved, unresolved = 0, []
                    tax_resolution_by_index = {}
                    with st.spinner("Resolving tax IDs…"):
                        for idx, row in df_mod.iterrows():
                            org = str(row.get("organism", "")).strip()
                            tid = str(row.get("tax_id", "")).strip()
                            if org and tid in ("", "<NA>", "None", "nan"):
                                hits = _ena_taxonomy_search(org)
                                if hits:
                                    exact = next(
                                        (
                                            h for h in hits
                                            if h.get("scientificName", "").lower() == org.lower()
                                        ),
                                        hits[0],
                                    )
                                    df_mod.at[idx, "tax_id"] = str(exact["taxId"])
                                    tax_resolution_by_index[str(idx)] = {
                                        "query": org,
                                        "selected_scientific_name": exact.get(
                                            "scientificName", ""
                                        ),
                                        "selected_tax_id": exact.get("taxId", ""),
                                        "rank": exact.get("rank", ""),
                                        "match": (
                                            "exact scientific name"
                                            if exact.get("scientificName", "").lower()
                                            == org.lower()
                                            else "first suggested candidate"
                                        ),
                                    }
                                    resolved += 1
                                else:
                                    unresolved.append(org)
                    tax_events = _record_metadata_changes(
                        st.session_state.metadata_df,
                        df_mod,
                        action="taxonomy_resolved",
                        source="ENA taxonomy API",
                        method="batch taxonomy suggestion",
                    )
                    for event in tax_events:
                        row_index = event.get("details", {}).get("row_index", "")
                        event["details"].update(
                            tax_resolution_by_index.get(str(row_index), {})
                        )
                    st.session_state.metadata_df = df_mod
                    st.session_state.pop("validated_df", None)
                    _clear_environment_resolver_state(clear_history=True)
                    st.session_state.editor_key += 1
                    if resolved:
                        st.toast(f"Resolved {resolved} tax ID(s).", icon="✅")
                    if unresolved:
                        st.warning(
                            f"Could not resolve {len(unresolved)} organism(s): "
                            + ", ".join(unresolved[:5])
                            + ("…" if len(unresolved) > 5 else "")
                        )
                    st.rerun()

            results = st.session_state.get("tax_results")
            if results is not None:
                if results:
                    res_df = pd.DataFrame([
                        {
                            "tax_id": r.get("taxId", ""),
                            "Scientific name": r.get("scientificName", ""),
                            "Rank": r.get("rank", ""),
                        }
                        for r in results[:10]
                    ])
                    st.dataframe(res_df, hide_index=True, width="stretch")

                    col_sel, col_fill = st.columns([3, 2])
                    with col_sel:
                        selected_idx = st.selectbox(
                            "Apply result",
                            options=range(len(results[:10])),
                            format_func=lambda i: (
                                f"{results[i].get('scientificName')} "
                                f"(taxId: {results[i].get('taxId')})"
                            ),
                            key="tax_select",
                        )
                    with col_fill:
                        st.write("")
                        if st.button("Fill all rows", key="btn_tax_fill", width="stretch"):
                            chosen = results[selected_idx]
                            df_mod = st.session_state.metadata_df.copy()
                            df_mod["organism"] = str(chosen["scientificName"])
                            df_mod["tax_id"] = str(chosen["taxId"])
                            _record_metadata_changes(
                                st.session_state.metadata_df,
                                df_mod,
                                action="taxonomy_selected",
                                source="ENA taxonomy API",
                                method="user-selected taxonomy candidate",
                                decision="user approved",
                                details={
                                    "query": tax_query,
                                    "selected_scientific_name": chosen.get(
                                        "scientificName", ""
                                    ),
                                    "selected_tax_id": chosen.get("taxId", ""),
                                    "rank": chosen.get("rank", ""),
                                },
                            )
                            st.session_state.metadata_df = df_mod
                            st.session_state.pop("validated_df", None)
                            st.session_state.pop("tax_results", None)
                            _clear_environment_resolver_state(clear_history=True)
                            st.session_state.editor_key += 1
                            st.toast("Filled organism and tax_id for all rows.", icon="✅")
                            st.rerun()
                else:
                    st.info("No results found. Try a broader search term.")

        # ── Environmental context resolver ───────────────────────────────────
        with tab_envo:
            resolver_accessions = sorted(
                {
                    accession
                    for value in st.session_state.metadata_df.get(
                        "sample derived from", pd.Series(dtype=str)
                    )
                    if (accession := _extract_sample_accession(value))
                }
            )
            st.caption(
                "Resolve environmental metadata in batches. MAG2ENA constrains biome and "
                "environmental-material terms to their ENVO branches, and checks local-context "
                "candidates for category conflicts before you review and apply them. Explicit "
                "ENVO identifiers in parent-sample titles and descriptions are reused only when "
                "their surrounding language identifies the corresponding MIxS field."
            )

            scope = st.radio(
                "Resolution scope",
                ["Fill empty fields", "Review and normalize all fields"],
                horizontal=True,
                key="envo_scope_control",
                help=(
                    "The default never overwrites populated cells. Normalization also checks "
                    "existing labels and identifiers and may replace accepted values with the "
                    "canonical `preferred label [ontology:identifier]` representation."
                ),
            )
            if (
                st.session_state.get("envo_proposals")
                and st.session_state.get("envo_resolver_scope") != scope
            ):
                _clear_environment_resolver_state()
            empty_only = scope == "Fill empty fields"
            resolution_groups = build_resolution_groups(
                st.session_state.metadata_df,
                empty_only=empty_only,
            )
            st.caption(
                f"{len(resolution_groups)} distinct context group(s) across "
                f"{len(st.session_state.metadata_df)} metadata row(s). Identical evidence from "
                "the same parent sample is resolved once and applied as a group."
            )

            col_resolve, col_undo = st.columns([3, 1])
            with col_resolve:
                resolve_clicked = st.button(
                    "Generate batch proposals",
                    key="btn_envo_batch_resolve",
                    type="primary",
                    width="stretch",
                    disabled=not resolution_groups,
                )
            with col_undo:
                if st.button(
                    "Undo last apply",
                    key="btn_envo_undo",
                    width="stretch",
                    disabled="envo_undo_df" not in st.session_state,
                ):
                    restored_df = st.session_state.pop("envo_undo_df")
                    _record_metadata_changes(
                        st.session_state.metadata_df,
                        restored_df,
                        action="environment_resolution_undone",
                        source="environmental resolver",
                        method="undo last curator-approved apply",
                        decision="user reversed",
                    )
                    st.session_state.metadata_df = restored_df
                    st.session_state.pop("validated_df", None)
                    _clear_environment_resolver_state(clear_history=True)
                    st.session_state.editor_key += 1
                    st.toast("Restored metadata from before the last ENVO apply.", icon="↩️")
                    st.rerun()

            if resolve_clicked:
                with st.spinner(
                    "Reading parent evidence, searching field-specific ontology scopes, "
                    "and ranking candidates…"
                ):
                    parent_metadata_by_accession = {}
                    failed_parent_evidence = []
                    fallback_parent_evidence = []
                    for accession in resolver_accessions:
                        raw = _fetch_ena_sample_metadata(accession)
                        if "_error" in raw:
                            failed_parent_evidence.append(accession)
                        else:
                            parent_metadata_by_accession[accession] = raw
                            if raw.get("_fallback_assembly_accessions"):
                                fallback_parent_evidence.append((accession, raw))
                    if fallback_parent_evidence:
                        _append_provenance(
                            [
                                make_event(
                                    stage="metadata",
                                    action="parent_metadata_recovered_from_assembly",
                                    entity_type="parent_sample",
                                    field="sample derived from",
                                    source=raw.get("_metadata_source", "ENA Portal API"),
                                    before=raw.get("_direct_sample_error", ""),
                                    after=accession,
                                    method="exact linked-assembly metadata fallback",
                                    method_version=MAG2ENA_VERSION,
                                    decision="used agreed indexed values",
                                    details={
                                        "parent_accession": accession,
                                        "source_url": raw.get("_metadata_source_url", ""),
                                        "assembly_accessions": raw.get(
                                            "_fallback_assembly_accessions", []
                                        ),
                                        "omitted_conflicts": raw.get(
                                            "_fallback_conflicts", {}
                                        ),
                                        "purpose": "environmental resolver evidence",
                                    },
                                )
                                for accession, raw in fallback_parent_evidence
                            ]
                        )
                    client = OlsEnvoClient()
                    proposals = resolve_environment_batch(
                        st.session_state.metadata_df,
                        client,
                        empty_only=empty_only,
                        parent_metadata_by_accession=parent_metadata_by_accession,
                    )
                if failed_parent_evidence:
                    st.warning(
                        "Parent title/description evidence could not be retrieved for: "
                        + ", ".join(f"`{value}`" for value in failed_parent_evidence)
                        + ". Proposals for those parents use the existing table evidence only."
                    )
                if fallback_parent_evidence:
                    st.info(
                        "Parent evidence was recovered through exact linked assemblies for: "
                        + ", ".join(
                            f"`{accession}`"
                            for accession, _ in fallback_parent_evidence
                        )
                        + ". Conflicting indexed values, if any, were omitted."
                    )
                st.session_state.envo_proposals = proposals
                st.session_state.envo_review_df = proposals_to_review_dataframe(proposals)
                _append_provenance(
                    [
                        _environment_proposal_event(
                            proposal,
                            action="environment_term_proposed",
                            decision="pending curator review",
                        )
                        for proposal in proposals
                    ]
                )
                st.session_state.envo_source_fingerprint = _metadata_fingerprint(
                    st.session_state.metadata_df
                )
                st.session_state.envo_resolver_scope = scope
                st.session_state.envo_review_key = (
                    st.session_state.get("envo_review_key", 0) + 1
                )

            proposals = st.session_state.get("envo_proposals")
            if proposals:
                counts = pd.Series([p["confidence"] for p in proposals]).value_counts()
                summary = " · ".join(f"{name}: {count}" for name, count in counts.items())
                st.info(f"Generated {len(proposals)} proposal(s) · {summary}")
                st.caption(
                    "No proposal is applied automatically. Check `accept` only after reviewing "
                    "the source evidence, field constraint, confidence, and selected value."
                )

                review_df = st.data_editor(
                    st.session_state.envo_review_df,
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
                    column_config={
                        "accept": st.column_config.CheckboxColumn(
                            "Accept", default=False, width="small"
                        ),
                        "proposal_id": st.column_config.TextColumn("ID", width="small"),
                        "rows": st.column_config.NumberColumn("Rows", width="small"),
                        "score": st.column_config.NumberColumn(
                            "Score", min_value=0.0, max_value=1.0, format="%.3f"
                        ),
                    },
                    key=f"envo_review_editor_{st.session_state.get('envo_review_key', 0)}",
                )
                st.session_state.envo_review_df = review_df

                # Keep decisions with the proposal records so alternative selection and
                # audit export survive Streamlit reruns.
                decisions = dict(zip(review_df["proposal_id"], review_df["accept"]))
                for proposal in proposals:
                    proposal["accept"] = bool(decisions.get(proposal["proposal_id"], False))
                st.session_state.envo_proposals = proposals

                col_accept_safe, col_clear = st.columns(2)
                with col_accept_safe:
                    if st.button(
                        "Select verified and high-confidence",
                        key="btn_envo_accept_safe",
                        width="stretch",
                    ):
                        safe = review_df["confidence"].isin(["Verified", "High"])
                        has_value = review_df["proposed value"].astype(str).str.len() > 0
                        review_df.loc[safe & has_value, "accept"] = True
                        st.session_state.envo_review_df = review_df
                        st.session_state.envo_review_key += 1
                        st.rerun()
                with col_clear:
                    if st.button(
                        "Clear selections",
                        key="btn_envo_clear_accept",
                        width="stretch",
                    ):
                        review_df["accept"] = False
                        st.session_state.envo_review_df = review_df
                        st.session_state.envo_review_key += 1
                        st.rerun()

                candidate_proposals = [p for p in proposals if p["candidates"]]
                if candidate_proposals:
                    with st.expander("Inspect alternatives or change a proposal"):
                        proposal_lookup = {
                            proposal["proposal_id"]: proposal
                            for proposal in candidate_proposals
                        }
                        selected_proposal_id = st.selectbox(
                            "Proposal",
                            options=list(proposal_lookup),
                            format_func=lambda proposal_id: (
                                f"{proposal_id} · {proposal_lookup[proposal_id]['field']} · "
                                f"{proposal_lookup[proposal_id]['row_count']} row(s)"
                            ),
                            key="envo_alternative_proposal",
                        )
                        selected_proposal = proposal_lookup[selected_proposal_id]
                        candidate_options = selected_proposal["candidates"]
                        candidate_index = st.selectbox(
                            "Candidate",
                            options=range(len(candidate_options)),
                            format_func=lambda index: (
                                f"{candidate_options[index]['value']} · "
                                f"{candidate_options[index]['score']:.3f} · "
                                f"{candidate_options[index]['reason']}"
                            ),
                            key=f"envo_alternative_candidate_{selected_proposal_id}",
                        )
                        chosen = candidate_options[candidate_index]
                        if chosen.get("description"):
                            st.caption(chosen["description"][:500])
                        if st.button(
                            "Use this candidate",
                            key="btn_envo_use_alternative",
                            width="stretch",
                        ):
                            previous_value = selected_proposal.get("selected_value", "")
                            selected_proposal["selected_value"] = chosen["value"]
                            selected_proposal["selected_label"] = chosen["label"]
                            selected_proposal["selected_id"] = chosen["obo_id"]
                            selected_proposal["score"] = chosen["score"]
                            selected_proposal["confidence"] = "User selected"
                            selected_proposal["reason"] = chosen["reason"]
                            _append_provenance(
                                [
                                    _environment_proposal_event(
                                        selected_proposal,
                                        action="environment_alternative_selected",
                                        decision="curator selected candidate",
                                        before=previous_value,
                                    )
                                ]
                            )
                            st.session_state.envo_proposals = proposals
                            replacement = proposals_to_review_dataframe(proposals)
                            previous_decisions = dict(
                                zip(review_df["proposal_id"], review_df["accept"])
                            )
                            replacement["accept"] = replacement["proposal_id"].map(
                                previous_decisions
                            ).fillna(False)
                            st.session_state.envo_review_df = replacement
                            st.session_state.envo_review_key += 1
                            st.rerun()

                accepted_ids = review_df.loc[review_df["accept"], "proposal_id"].tolist()
                if st.button(
                    f"Apply {len(accepted_ids)} accepted proposal(s)",
                    key="btn_envo_apply",
                    type="primary",
                    width="stretch",
                    disabled=not accepted_ids,
                ):
                    current_fingerprint = _metadata_fingerprint(st.session_state.metadata_df)
                    if current_fingerprint != st.session_state.get("envo_source_fingerprint"):
                        st.error(
                            "The metadata table changed after these proposals were generated. "
                            "Generate a fresh batch before applying them."
                        )
                    else:
                        _append_provenance(
                            [
                                _environment_proposal_event(
                                    proposal,
                                    action="environment_proposal_reviewed",
                                    decision=(
                                        "curator approved"
                                        if proposal.get("proposal_id") in accepted_ids
                                        and proposal.get("selected_value")
                                        else "not selected"
                                    ),
                                )
                                for proposal in proposals
                            ]
                        )
                        st.session_state.envo_undo_df = st.session_state.metadata_df.copy()
                        df_mod, audit = apply_environment_proposals(
                            st.session_state.metadata_df,
                            proposals,
                            accepted_ids,
                            overwrite=(
                                st.session_state.get("envo_resolver_scope")
                                == "Review and normalize all fields"
                            ),
                        )
                        resolver_events = _record_metadata_changes(
                            st.session_state.metadata_df,
                            df_mod,
                            action="environment_term_applied",
                            source="OLS4 / ENVO",
                            method="field-aware environmental resolver",
                            decision="curator approved",
                        )
                        audit_lookup = {}
                        for item in audit:
                            for sample_name in str(item.get("sample_names", "")).split(", "):
                                audit_lookup[(sample_name, item.get("field", ""))] = item
                        for event in resolver_events:
                            audit_item = audit_lookup.get(
                                (event.get("sample_name", ""), event.get("field", "")),
                                {},
                            )
                            event["details"].update(
                                {
                                    "proposal_id": audit_item.get("proposal_id", ""),
                                    "parent_accession": audit_item.get(
                                        "parent_accession", ""
                                    ),
                                    "source_evidence": audit_item.get(
                                        "source_evidence", ""
                                    ),
                                    "confidence": audit_item.get("confidence", ""),
                                    "score": audit_item.get("score", ""),
                                    "reason": audit_item.get("reason", ""),
                                    "ontology_version": audit_item.get(
                                        "ontology_version", ""
                                    ),
                                    "alternatives": audit_item.get("alternatives", ""),
                                }
                            )
                        st.session_state.metadata_df = df_mod
                        st.session_state.pop("validated_df", None)
                        _clear_environment_resolver_state()
                        st.session_state.editor_key += 1
                        changed_cells = sum(
                            len(str(item["rows"]).split(", ")) for item in audit
                        )
                        st.toast(
                            f"Applied {len(audit)} resolution group(s) to "
                            f"{changed_cells} cell(s).",
                            icon="✅",
                        )
                        st.rerun()

        # ── Fill column ───────────────────────────────────────────────────────
        with tab_fill:
            st.caption(
                "Apply a single value across an entire column — useful for study-wide "
                "constants such as `platform`, `assembly software`, or `collection date`."
            )
            fillable_cols = [c for c in st.session_state.metadata_df.columns if c != "sample_name"]
            fill_col = st.selectbox("Column", options=fillable_cols, key="fill_col")
            fill_empty_only = st.checkbox("Fill empty cells only", value=True, key="fill_empty_only")

            field_def = field_defs.get(fill_col, {})
            if field_def.get("type") == "enum":
                fill_val = st.selectbox("Value", options=field_def["enum"], key="fill_val_enum")
            else:
                fill_val = st.text_input("Value", key="fill_val_text")

            if st.button("Apply to column", key="btn_fill_apply", type="primary"):
                df_mod = st.session_state.metadata_df.copy()
                if fill_empty_only:
                    mask = df_mod[fill_col].isna() | (
                        df_mod[fill_col].astype(str).str.strip().isin(["", "<NA>", "None", "nan"])
                    )
                    df_mod.loc[mask, fill_col] = fill_val
                else:
                    df_mod[fill_col] = fill_val
                _record_metadata_changes(
                    st.session_state.metadata_df,
                    df_mod,
                    action="column_filled",
                    source="Fill column tool",
                    method="bulk metadata update",
                    decision="user applied",
                    details={"empty_cells_only": fill_empty_only},
                )
                st.session_state.metadata_df = df_mod
                st.session_state.pop("validated_df", None)
                _clear_environment_resolver_state(clear_history=True)
                st.session_state.editor_key += 1
                st.toast(f"Applied value to `{fill_col}`.", icon="✅")
                st.rerun()

        # ── Import quality files ──────────────────────────────────────────────
        with tab_import:
            st.caption(
                "Auto-populate metadata columns from bioinformatics tool outputs. "
                "Only empty cells are overwritten — existing values are preserved."
            )
            col_checkm, col_gtdbtk = st.columns(2)

            with col_checkm:
                col_title, col_ex1, col_ex2 = st.columns([3, 2, 2])
                with col_title:
                    st.markdown("**CheckM / CheckM2**")
                with col_ex1:
                    checkm1_path = EXAMPLES_DIR / "checkm_example.tsv"
                    if checkm1_path.exists():
                        st.download_button(
                            "CheckM example",
                            data=checkm1_path.read_bytes(),
                            file_name="checkm_example.tsv",
                            mime="text/tab-separated-values",
                            width="stretch",
                            help="Download an example CheckM v1 output (storage.tsv format).",
                        )
                with col_ex2:
                    checkm2_path = EXAMPLES_DIR / "checkm2_example.tsv"
                    if checkm2_path.exists():
                        st.download_button(
                            "CheckM2 example",
                            data=checkm2_path.read_bytes(),
                            file_name="checkm2_example.tsv",
                            mime="text/tab-separated-values",
                            width="stretch",
                            help="Download an example CheckM2 output (quality_report.tsv format).",
                        )
                st.caption("Fills: `completeness score`, `contamination score`, `completeness software`")
                checkm_file = st.file_uploader(
                    "quality_report.tsv or storage.tsv",
                    type=["tsv", "txt"],
                    key="checkm_upload",
                    help=(
                        "Upload CheckM2 quality_report.tsv or CheckM storage.tsv. "
                        "Bin names must match values in the sample_name column."
                    ),
                )
                if checkm_file:
                    parsed_qc = _parse_checkm_file(checkm_file)
                    if parsed_qc is None:
                        st.error("Unrecognised format. Expected CheckM or CheckM2 output.")
                    elif st.button("Apply CheckM data", key="btn_checkm"):
                        checkm_bytes = checkm_file.getvalue()
                        df_mod = _merge_into_metadata(st.session_state.metadata_df, parsed_qc)
                        matched = int(parsed_qc["sample_name"].isin(df_mod["sample_name"]).sum())
                        detected_software = " | ".join(
                            sorted(parsed_qc["completeness software"].unique())
                        )
                        _append_provenance(
                            [
                                make_event(
                                    stage="metadata",
                                    action="quality_file_imported",
                                    entity_type="file",
                                    source=checkm_file.name,
                                    after=sha256_bytes(checkm_bytes),
                                    method="CheckM output parser",
                                    method_version=MAG2ENA_VERSION,
                                    details={
                                        "filename": checkm_file.name,
                                        "size_bytes": len(checkm_bytes),
                                        "sha256": sha256_bytes(checkm_bytes),
                                        "detected_software": detected_software,
                                        "matched_samples": matched,
                                    },
                                )
                            ]
                        )
                        _record_metadata_changes(
                            st.session_state.metadata_df,
                            df_mod,
                            action="quality_metrics_imported",
                            source=checkm_file.name,
                            method="CheckM/CheckM2 parser",
                            details={
                                "source_file_sha256": sha256_bytes(checkm_bytes),
                                "detected_software": detected_software,
                            },
                        )
                        st.session_state.metadata_df = df_mod
                        st.session_state.pop("validated_df", None)
                        _clear_environment_resolver_state(clear_history=True)
                        st.session_state.editor_key += 1
                        st.toast(f"Applied quality scores for {matched} sample(s).", icon="✅")
                        st.rerun()

            with col_gtdbtk:
                col_title, col_ex = st.columns([3, 2])
                with col_title:
                    st.markdown("**GTDB-Tk**")
                with col_ex:
                    gtdbtk_path = EXAMPLES_DIR / "gtdbtk_example.tsv"
                    if gtdbtk_path.exists():
                        st.download_button(
                            "GTDB-Tk example",
                            data=gtdbtk_path.read_bytes(),
                            file_name="gtdbtk_example.tsv",
                            mime="text/tab-separated-values",
                            width="stretch",
                            help="Download an example GTDB-Tk summary TSV (bac120 format).",
                        )
                st.caption("Fills: `organism` — use *Taxonomy resolver* afterwards to auto-fill `tax_id`")
                gtdbtk_file = st.file_uploader(
                    "gtdbtk.bac120.summary.tsv or gtdbtk.ar53.summary.tsv",
                    type=["tsv"],
                    key="gtdbtk_upload",
                    help=(
                        "Upload a GTDB-Tk bacterial or archaeal summary. The user_genome "
                        "values must match values in the sample_name column."
                    ),
                )
                if gtdbtk_file:
                    parsed_gtdb = _parse_gtdbtk_file(gtdbtk_file)
                    if parsed_gtdb is None:
                        st.error("Unrecognised format. Expected GTDB-Tk summary TSV.")
                    elif st.button("Apply GTDB-Tk data", key="btn_gtdbtk"):
                        gtdbtk_bytes = gtdbtk_file.getvalue()
                        df_mod = _merge_into_metadata(st.session_state.metadata_df, parsed_gtdb)
                        matched = int(parsed_gtdb["sample_name"].isin(df_mod["sample_name"]).sum())
                        _append_provenance(
                            [
                                make_event(
                                    stage="metadata",
                                    action="taxonomy_file_imported",
                                    entity_type="file",
                                    source=gtdbtk_file.name,
                                    after=sha256_bytes(gtdbtk_bytes),
                                    method="GTDB-Tk output parser",
                                    method_version=MAG2ENA_VERSION,
                                    details={
                                        "filename": gtdbtk_file.name,
                                        "size_bytes": len(gtdbtk_bytes),
                                        "sha256": sha256_bytes(gtdbtk_bytes),
                                        "matched_samples": matched,
                                    },
                                )
                            ]
                        )
                        _record_metadata_changes(
                            st.session_state.metadata_df,
                            df_mod,
                            action="organism_names_imported",
                            source=gtdbtk_file.name,
                            method="GTDB-Tk rank-aware parser",
                            details={
                                "source_file_sha256": sha256_bytes(gtdbtk_bytes)
                            },
                        )
                        st.session_state.metadata_df = df_mod
                        st.session_state.pop("validated_df", None)
                        _clear_environment_resolver_state(clear_history=True)
                        st.session_state.editor_key += 1
                        st.toast(
                            f"Applied organism names for {matched} sample(s). "
                            "Run Taxonomy resolver to fill tax IDs.",
                            icon="✅",
                        )
                        st.rerun()

    with st.expander("Provenance and history"):
        provenance_events = _ensure_provenance_state()
        history_df = events_to_dataframe(provenance_events)
        st.caption(
            "A unified, append-only history of metadata changes, curator decisions, "
            "validation, file processing, and submission actions. Credentials are never recorded."
        )
        if history_df.empty:
            st.info("No provenance events have been recorded yet.")
        else:
            col_action, col_sample = st.columns(2)
            with col_action:
                selected_actions = st.multiselect(
                    "Filter actions",
                    options=sorted(history_df["action"].dropna().unique()),
                    key="provenance_action_filter",
                )
            with col_sample:
                sample_options = sorted(
                    value
                    for value in history_df["sample_name"].dropna().unique()
                    if value
                )
                selected_samples = st.multiselect(
                    "Filter MAGs",
                    options=sample_options,
                    key="provenance_sample_filter",
                )
            filtered_history = history_df
            if selected_actions:
                filtered_history = filtered_history[
                    filtered_history["action"].isin(selected_actions)
                ]
            if selected_samples:
                filtered_history = filtered_history[
                    filtered_history["sample_name"].isin(selected_samples)
                ]
            display_columns = [
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
                filtered_history[display_columns],
                hide_index=True,
                width="stretch",
                height=min(480, 80 + 35 * max(len(filtered_history), 1)),
            )
            col_tsv, col_json = st.columns(2)
            with col_tsv:
                st.download_button(
                    "Download provenance.tsv",
                    data=events_to_tsv(provenance_events).encode("utf-8"),
                    file_name="provenance.tsv",
                    mime="text/tab-separated-values",
                    width="stretch",
                    key="download_provenance_tsv",
                )
            with col_json:
                st.download_button(
                    "Download provenance.json",
                    data=events_to_json(provenance_events).encode("utf-8"),
                    file_name="provenance.json",
                    mime="application/json",
                    width="stretch",
                    key="download_provenance_json",
                )

    col_validate, col_reset = st.columns([4, 1])

    with col_reset:
        if st.button("Reset table", width="stretch", help="Clear all rows and start over."):
            _reset_provenance(action="workflow_reset", source="Reset table button")
            st.session_state.metadata_df = initialize_empty_dataframe(field_defs)
            st.session_state.pop("tsv_filename", None)
            st.session_state.pop("validated_df", None)
            st.session_state.pop("_fasta_ready", None)
            _clear_environment_resolver_state(clear_history=True)
            st.session_state.editor_key += 1
            st.rerun()

    with col_validate:
        validate_clicked = st.button(
            "Validate metadata",
            type="primary",
            width="stretch",
            help=(
                "Check mandatory fields, formats, controlled vocabularies, and taxon IDs "
                "before continuing to FASTA upload. MAG tax IDs assigned to metagenome "
                "taxonomies are blocked."
            ),
        )

    if validate_clicked:

        non_empty = edited_df.dropna(how="all")

        if edited_df.empty or non_empty.empty:
            _record_validation(
                edited_df,
                "failed",
                [{"expected": "Metadata table cannot be empty"}],
            )
            st.error("Metadata table cannot be empty.")
            st.stop()

        if len(edited_df) > 1000:
            _record_validation(
                edited_df,
                "failed",
                [{"expected": "At most 1000 metadata rows per web submission"}],
            )
            st.error(
                f"Submission is limited to 1000 rows. "
                f"Your table has {len(edited_df)} rows."
            )
            st.stop()

        error_df = validate_dataframe(edited_df, field_defs)

        if not error_df.empty:
            _record_validation(edited_df, "failed", error_df)
            st.error(f"{len(error_df)} validation error(s) found:")
            st.dataframe(
                error_df,
                hide_index=True,
                width="stretch",
            )
            st.stop()

        st.success(
            f"All {len(edited_df)} sample(s) validated successfully."
        )

        _record_validation(edited_df, "passed")
        st.session_state.validated_df = edited_df

    if "validated_df" not in st.session_state:
        return

    # ── Step 2: FASTA files ───────────────────────────────────────────────────
    st.divider()

    st.subheader("Step 2 — FASTA files")

    n_samples = len(st.session_state.validated_df)

    fasta_map = {}

    if use_example:
        st.caption(
            f"Example FASTA files are loaded automatically — "
            f"no upload needed when using example data."
        )

        # Build temp copies of example FASTAs so the originals are never deleted by the worker.
        if "_example_fasta_map" not in st.session_state or not all(
            Path(p).exists() for p in st.session_state._example_fasta_map.values()
        ):
            fasta_dir = EXAMPLES_DIR / "fasta"
            tmp_map = {}
            for fp in sorted(fasta_dir.glob("*.fasta.gz")):
                sample_name = fp.name.removesuffix(".fasta.gz")
                tmp = tempfile.NamedTemporaryFile(suffix=".fasta.gz", delete=False)
                tmp.write(fp.read_bytes())
                tmp.flush()
                tmp.close()
                tmp_map[sample_name] = str(Path(tmp.name).resolve())
            st.session_state._example_fasta_map = tmp_map

        fasta_map = {k: Path(v) for k, v in st.session_state._example_fasta_map.items()}
        names = ", ".join(f"`{k}.fasta.gz`" for k in sorted(fasta_map))
        st.success(f"Example FASTA files loaded: {names}")

    else:
        st.caption(
            f"Upload one `.fasta.gz` file per sample. "
            f"Expected {n_samples} file(s) — filenames must match the `sample_name` column."
        )

        uploaded_fastas = st.file_uploader(
            "Upload FASTA.GZ files",
            type=["gz"],
            accept_multiple_files=True,
            help=(
                "Upload one gzip-compressed .fasta.gz file per metadata row. Each filename "
                "without .fasta.gz must exactly match its sample_name value."
            ),
        )

        if uploaded_fastas:

            try:
                fasta_map = persist_fastas_temp(uploaded_fastas)
            except Exception as e:
                st.error(str(e))
                st.stop()

            expected = set(
                st.session_state.validated_df["sample_name"].dropna()
            )
            uploaded = set(fasta_map)
            missing = expected - uploaded
            extra = uploaded - expected

            if missing:
                st.error(
                    f"Missing FASTA files for: {', '.join(sorted(missing))}"
                )
                st.stop()

            if extra:
                st.warning(
                    f"Extra FASTA files not in the metadata (will be ignored): "
                    f"{', '.join(sorted(extra))}"
                )
                for sample_name in extra:
                    extra_path = fasta_map.pop(sample_name, None)
                    if extra_path is not None:
                        Path(extra_path).unlink(missing_ok=True)

            st.success(
                f"All {len(fasta_map)} FASTA file(s) uploaded and matched successfully."
            )

    if not fasta_map:
        st.session_state.pop("_fasta_ready", None)
        return

    _record_fasta_inputs(
        fasta_map,
        source="bundled example" if use_example else "browser upload",
    )

    with st.spinner("Validating FASTA headers…"):
        header_violations = _validate_fasta_headers(fasta_map)

    if header_violations:
        _append_provenance(
            [
                make_event(
                    stage="validation",
                    action="fasta_headers_validated",
                    entity_type="file_collection",
                    source="uploaded FASTA files",
                    after="failed",
                    method="ENA sequence-name preflight",
                    method_version=MAG2ENA_VERSION,
                    decision="failed",
                    details={
                        "maximum_sequence_name_length": 50,
                        "violations": header_violations,
                    },
                )
            ]
        )
        st.error(
            f"{len(header_violations)} sequence name(s) exceed 50 characters "
            f"(ENA limit). Rename them in your FASTA files and re-upload."
        )
        st.dataframe(
            pd.DataFrame(header_violations),
            hide_index=True,
            width="stretch",
        )
        return

    try:
        fasta_inspections = inspect_fasta_collection(fasta_map)
    except ValueError as exc:
        st.error(str(exc))
        return

    single_contig_names = single_contig_samples(fasta_inspections)
    single_contig_chromosomes_confirmed = not single_contig_names
    if single_contig_names:
        formatted_names = ", ".join(
            f"`{sample_name}`" for sample_name in single_contig_names
        )
        st.warning(
            f"Single-contig MAG(s) detected: {formatted_names}. "
            f"{ENA_SINGLE_CONTIG_HELP}",
            icon="⚠️",
        )
        single_contig_chromosomes_confirmed = st.checkbox(
            "I confirm that every single-contig MAG listed above represents a "
            "complete chromosome.",
            help=(
                "MAG2ENA will submit each confirmed sequence as a chromosome-level "
                "assembly and generate its ENA chromosome list. Do not confirm for "
                "an incomplete single-contig assembly; contact the ENA helpdesk instead."
            ),
        )
        if not single_contig_chromosomes_confirmed:
            st.info(
                "Submission is paused until the complete-chromosome assertion is "
                "confirmed. For an incomplete one-sequence MAG, contact "
                "datasubs@ebi.ac.uk to request a single-contig exception."
            )
            return

    fasta_validation_key = tuple(
        sorted(st.session_state.get("_provenance_fasta_hashes", {}).items())
    )
    if st.session_state.get("_provenance_validated_fasta_key") != fasta_validation_key:
        _append_provenance(
            [
                make_event(
                    stage="validation",
                    action="fasta_headers_validated",
                    entity_type="file_collection",
                    source="uploaded FASTA files",
                    after="passed",
                    method="ENA sequence-name preflight",
                    method_version=MAG2ENA_VERSION,
                    decision="passed",
                    details={
                        "maximum_sequence_name_length": 50,
                        "file_count": len(fasta_map),
                    },
                )
            ]
        )
        st.session_state._provenance_validated_fasta_key = fasta_validation_key

    st.session_state["_fasta_ready"] = True
    _step_placeholder.markdown(_step_progress_html(3), unsafe_allow_html=True)

    # ── Step 3: ENA submission ────────────────────────────────────────────────
    st.divider()
    st.subheader("Step 3 — ENA submission")

    st.warning(
        "**New ENA Webin account?** Accounts typically require at least **24 hours** after "
        "registration before they can submit data. If you receive an authentication or "
        "authorisation error, wait 24 h and try again.",
        icon="⚠️",
    )

    study_mode = st.radio(
        "Submission mode",
        ["Create new study", "Existing study accession"],
        horizontal=True,
        help=(
            "Create and register a new ENA study with this submission, or attach the MAGs "
            "to a study that already has an ENA project accession."
        ),
    )

    notification_container = st.container()

    with st.form("submission_form", clear_on_submit=True):

        if study_mode == "Create new study":

            study_name = st.text_input(
                "Study name",
                help=(
                    "Short alias used to identify the new study in ENA. Choose a name "
                    "that is unique within your Webin account."
                ),
            )
            study_title = st.text_input(
                "Study title",
                max_chars=250,
                help="Public-facing title for the ENA study; enter 20 to 250 characters.",
            )
            study_description = st.text_area(
                "Study description",
                max_chars=4000,
                help=(
                    "Public summary of the study's purpose, samples, and context; "
                    "enter 20 to 4,000 characters."
                ),
            )
            study_accession = None

        else:

            study_accession = st.text_input(
                "Study accession",
                placeholder="e.g. PRJEB12345",
                help=(
                    "ENA project accession for an existing study to which your Webin "
                    "account has submission access."
                ),
            )
            study_name = None
            study_title = None
            study_description = None

        st.info(
            "Your Webin credentials are used only to authenticate this submission. "
            "They are not written to MAG2ENA's database, job directory, or downloadable logs. "
            "The temporary queue copy is removed as soon as the worker starts processing the job.",
            icon="🔒",
        )

        col_user, col_pass = st.columns(2)

        with col_user:
            ena_user = st.text_input(
                "ENA Webin username",
                placeholder="e.g. Webin-12345",
                help="Username of the ENA Webin account that will own the submission.",
            )

        with col_pass:
            ena_password = st.text_input(
                "ENA Webin password",
                type="password",
                help=(
                    "Password used to authenticate this submission. It is not written "
                    "to MAG2ENA's permanent database or job output files."
                ),
            )

        portal = st.radio(
            "Submission portal",
            ["Testing", "Production"],
            horizontal=True,
            help=(
                "Testing sends the job to ENA's test service. Use it first, then choose "
                "Production only when you are ready to create the real submission."
            ),
        )

        _max_release = date.today() + timedelta(days=730)
        release_date = st.date_input(
            "Public release date",
            value=_max_release,
            min_value=date.today() + timedelta(days=1),
            max_value=_max_release,
            format="YYYY-MM-DD",
            help=(
                "Date when submitted data become publicly available. "
                "ENA allows at most 2 years from today."
            ),
        )

        col_email, col_encrypt = st.columns(2)

        with col_email:
            notification_email = st.text_input(
                "Notification email (optional)",
                help="Enter your email to receive a notification when the submission completes."
            )

        with col_encrypt:
            encryption_password = st.text_input(
                "Password to encrypt submission (optional)",
                type="password",
                help=(
                    "Only with this password can the job results be accessed. "
                    "Not even the administrators can view encrypted submissions."
                ),
            )

        submitted = st.form_submit_button(
            "Submit to ENA",
            type="primary",
            width="stretch",
            help=(
                "Queue the validated metadata and matched FASTA files for submission "
                "to the selected ENA portal."
            ),
        )

    if submitted:

        if not encryption_password:
            encryption_password = None

        if study_mode == "Create new study":

            if not study_name:
                with notification_container:
                    st.error("Study name is required.")
                st.stop()

            title_len = len(study_title or "")
            description_len = len(study_description or "")

            if title_len < 20 or title_len > 250:
                with notification_container:
                    st.error(
                        "Study title must be between 20 and 250 characters."
                    )
                st.stop()

            if description_len < 20 or description_len > 4000:
                with notification_container:
                    st.error(
                        "Study description must be between 20 and 4,000 characters."
                    )
                st.stop()

        if not ena_user or not ena_password:
            with notification_container:
                st.error(
                    "Please enter a valid ENA Webin username and password."
                )
            st.stop()

        submission = {
            "study_name": study_name,
            "study_title": study_title,
            "study_description": study_description,
            "study_accession": study_accession,
            "ena_user": ena_user,
            "ena_password": ena_password,
            "portal": portal,
            "hold_until_date": release_date.strftime("%Y-%m-%d"),
            "single_contig_chromosomes_confirmed": (
                single_contig_chromosomes_confirmed
            ),
        }

        fasta_map_str = {k: str(v) for k, v in fasta_map.items()}

        df_records = pl.from_pandas(
            st.session_state.validated_df
        ).to_dicts()

        _append_provenance(
            [
                make_event(
                    stage="submission",
                    action="submission_queued",
                    entity_type="submission",
                    source=portal,
                    after="queued",
                    method="MAG2ENA asynchronous submission",
                    method_version=MAG2ENA_VERSION,
                    decision="user submitted",
                    details={
                        "portal": portal,
                        "study_mode": study_mode,
                        "study_accession": study_accession or "",
                        "study_name": study_name or "",
                        "hold_until_date": release_date.strftime("%Y-%m-%d"),
                        "sample_count": len(st.session_state.validated_df),
                        "single_contig_samples": single_contig_names,
                        "single_contig_chromosomes_confirmed": (
                            single_contig_chromosomes_confirmed
                        ),
                        "metadata_sha256": dataframe_sha256(
                            st.session_state.validated_df
                        ),
                        "fasta_sha256": st.session_state.get(
                            "_provenance_fasta_hashes", {}
                        ),
                    },
                )
            ]
        )

        fn_kwargs = {
            "df_records": df_records,
            "submission": submission,
            "fasta_map_str": fasta_map_str,
            "provenance_events": list(_ensure_provenance_state()),
        }

        if notification_email:
            fn_kwargs["email"] = notification_email

        if encryption_password:
            fn_kwargs["encryption_password"] = encryption_password

        job_id = enqueue_task(submission_task, fn_kwargs)

        _step_placeholder.markdown(_step_progress_html(4), unsafe_allow_html=True)

        with notification_container:
            st.success("Submission queued successfully!")
            job_submitted_dialog(job_id, encrypted=bool(encryption_password))

# =========================================================
# MAIN
# =========================================================

if __name__ == "__main__":
    runUI()
