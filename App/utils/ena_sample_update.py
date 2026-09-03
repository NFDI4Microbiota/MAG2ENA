"""Helpers for finding and updating existing ENA MAG metadata."""

from __future__ import annotations

import hashlib
import re
import time
import xml.etree.ElementTree as ET
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import quote

import requests

BATCH_SIZE = 1000
MAG_ASSEMBLY_TYPE = "Metagenome-Assembled Genome (MAG)"
PORTAL_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
_REPORT_XML_CHUNK_SIZE = 300
_REPORT_FETCH_WORKERS = 4

_STUDY_ACCESSION_RE = re.compile(
    r"(?:PRJ(?:EB|NA|DB)\d+|[ESD]RP\d+)", re.IGNORECASE
)
_SAMPLE_ACCESSION_RE = re.compile(
    r"(?:[ESD]RS\d{6,}|SAM[END][AG]?\d+|EGAN\d{11})", re.IGNORECASE
)
_EMPTY_VALUES = {"", "<na>", "nan", "none"}
_UNITS_FOR_NEW_ATTRIBUTES = {
    "completeness score": "%",
    "contamination score": "%",
    "geographic location (latitude)": "DD",
    "geographic location (longitude)": "DD",
}


class EnaSampleUpdateError(RuntimeError):
    """Raised when an ENA discovery, retrieval, or update request fails."""


# Accessions and study discovery
def normalize_study_accession(value: str) -> str:
    match = _STUDY_ACCESSION_RE.fullmatch(str(value or "").strip())
    return match.group(0).upper() if match else ""


def parse_sample_accessions(value: str | Iterable[str]) -> list[str]:
    """Return unique, normalized sample accessions in input order."""
    if isinstance(value, str):
        candidates = re.split(r"[\s,;]+", value)
    else:
        candidates = list(value)
    accessions = []
    seen = set()
    for candidate in candidates:
        text = str(candidate or "").strip().upper()
        if not text:
            continue
        if not _SAMPLE_ACCESSION_RE.fullmatch(text):
            raise ValueError(f"Unsupported ENA sample accession: {candidate}")
        if text not in seen:
            seen.add(text)
            accessions.append(text)
    return accessions


def _canonical_study_accession(
    study_accession: str,
    *,
    search_url: str,
    auth: tuple[str, str] | None,
    client,
    timeout: int,
) -> str:
    """Resolve an ERP/SRP/DRP accession to the BioProject used by assemblies."""
    if study_accession.startswith("PRJ"):
        return study_accession
    response = client.get(
        search_url,
        params={
            "result": "study",
            "query": f'secondary_study_accession="{study_accession}"',
            "fields": "study_accession",
            "format": "json",
            "limit": 10,
        },
        auth=auth,
        timeout=timeout,
    )
    response.raise_for_status()
    rows = response.json()
    projects = {
        str(row.get("study_accession", "") or "").strip().upper()
        for row in rows
        if isinstance(row, dict) and row.get("study_accession")
    }
    if len(projects) != 1:
        raise EnaSampleUpdateError(
            f"Could not resolve {study_accession} to one ENA BioProject accession"
        )
    return projects.pop()


def _portal_query_for_study(study_accession: str) -> str:
    escaped = study_accession.replace('"', "")
    return (
        f'study_accession="{escaped}" AND '
        f'assembly_type="{MAG_ASSEMBLY_TYPE}"'
    )


def _discover_mag_samples_from_reports(
    study_accession: str,
    *,
    portal: str,
    username: str,
    password: str,
    client,
    timeout: int,
) -> list[dict]:
    """Use the authenticated Reports service when Portal search is unavailable."""
    base = _reports_base(portal)
    auth = (username, password)
    response = client.get(
        f"{base}/analyses/{quote(study_accession)}",
        params={"format": "json", "max-results": 100_000},
        auth=auth,
        timeout=max(timeout, 120),
    )
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, list):
        raise EnaSampleUpdateError(
            "ENA Reports study discovery returned an unexpected response"
        )

    reports = []
    for item in payload:
        report = item.get("report", {}) if isinstance(item, dict) else {}
        if str(report.get("analysisType", "")) != "SEQUENCE_ASSEMBLY":
            continue
        analysis_accession = str(report.get("id", "") or "").strip().upper()
        if analysis_accession:
            reports.append((analysis_accession, report))

    reports_by_accession = dict(reports)
    batches = []
    for offset in range(0, len(reports), _REPORT_XML_CHUNK_SIZE):
        accessions = [
            accession
            for accession, _ in reports[offset : offset + _REPORT_XML_CHUNK_SIZE]
        ]
        batches.append((offset, accessions))

    def fetch_analysis_xml(accessions: list[str]) -> ET.Element:
        xml_response = client.get(
            f"{base}/analyses/xml/{quote(','.join(accessions), safe=',')}",
            auth=auth,
            timeout=max(timeout, 120),
        )
        xml_response.raise_for_status()
        try:
            root = ET.fromstring(xml_response.content)
        except ET.ParseError as exc:
            raise EnaSampleUpdateError(
                "ENA Reports returned invalid analysis XML"
            ) from exc
        return root

    # Four workers keep large private studies responsive while remaining below
    # ENA's recommendation of no more than ten simultaneous API searches.
    roots = {}
    with ThreadPoolExecutor(
        max_workers=min(_REPORT_FETCH_WORKERS, max(1, len(batches)))
    ) as executor:
        future_to_offset = {
            executor.submit(fetch_analysis_xml, accessions): offset
            for offset, accessions in batches
        }
        for future in as_completed(future_to_offset):
            roots[future_to_offset[future]] = future.result()

    discovered: dict[str, dict] = {}
    for offset in sorted(roots):
        root = roots[offset]
        for analysis in root.iter():
            if _local_name(analysis) != "ANALYSIS":
                continue
            analysis_accession = str(analysis.get("accession", "") or "").upper()
            report = reports_by_accession.get(analysis_accession, {})
            assembly = _find_first(analysis, "SEQUENCE_ASSEMBLY")
            assembly_type = (
                _find_direct(assembly, "TYPE") if assembly is not None else None
            )
            if (
                assembly_type is None
                or str(assembly_type.text or "").strip() != MAG_ASSEMBLY_TYPE
            ):
                continue
            sample_ref = _find_first(analysis, "SAMPLE_REF")
            sample_accession = str(
                sample_ref.get("accession", "") if sample_ref is not None else ""
            ).strip().upper()
            if not re.fullmatch(r"[ESD]RS\d+", sample_accession):
                continue
            name = _find_direct(assembly, "NAME")
            assembly_name = str(name.text or "").strip() if name is not None else ""
            item = discovered.setdefault(
                sample_accession,
                {
                    "sample_accession": sample_accession,
                    "biosample_accession": "",
                    "assembly_accessions": [],
                    "assembly_names": [],
                    "scientific_name": "",
                    "tax_id": "",
                    "study_accession": str(
                        report.get("projectId")
                        or report.get("studyId")
                        or study_accession
                    ),
                },
            )
            if analysis_accession not in item["assembly_accessions"]:
                item["assembly_accessions"].append(analysis_accession)
            if assembly_name and assembly_name not in item["assembly_names"]:
                item["assembly_names"].append(assembly_name)

    return sorted(discovered.values(), key=lambda item: item["sample_accession"])


def discover_mag_samples_by_study(
    study_accession: str,
    *,
    portal: str = "Production",
    username: str = "",
    password: str = "",
    client=requests,
    timeout: int = 60,
) -> list[dict]:
    """Return the distinct MAG samples linked to an ENA study.

    Discovery uses the assembly index and an exact MAG assembly-type filter,
    preventing environmental parent samples, read samples, isolates, and SAGs
    in the same study from entering the update editor.
    """
    accession = normalize_study_accession(study_accession)
    if not accession:
        raise ValueError(f"Unsupported ENA study accession: {study_accession}")

    # The study and assembly type are already fixed by the query. Scientific
    # name and taxonomy are loaded from authoritative sample XML with the batch,
    # so excluding them roughly halves the large discovery response.
    fields = (
        "assembly_accession,assembly_name,sample_accession,"
        "secondary_sample_accession"
    )
    auth = (username, password) if username and password else None
    search_url = (
        PORTAL_SEARCH_URL.replace("www.ebi.ac.uk", "wwwdev.ebi.ac.uk")
        if str(portal).casefold() == "testing"
        else PORTAL_SEARCH_URL
    )
    try:
        canonical_accession = _canonical_study_accession(
            accession,
            search_url=search_url,
            auth=auth,
            client=client,
            timeout=timeout,
        )
        response = client.get(
            search_url,
            params={
                "result": "assembly",
                "query": _portal_query_for_study(canonical_accession),
                "fields": fields,
                "format": "json",
                # ENA uses limit=0 for all results. The UI still loads batches of 1,000.
                "limit": 0,
            },
            auth=auth,
            timeout=timeout,
        )
        response.raise_for_status()
        payload = response.json()
    except (requests.RequestException, ValueError) as exc:
        if username and password:
            try:
                return _discover_mag_samples_from_reports(
                    accession,
                    portal=portal,
                    username=username,
                    password=password,
                    client=client,
                    timeout=timeout,
                )
            except (requests.RequestException, ValueError, EnaSampleUpdateError):
                pass
        raise EnaSampleUpdateError("ENA study discovery failed.") from exc

    if not isinstance(payload, list):
        raise EnaSampleUpdateError("ENA study discovery returned an unexpected response")

    discovered: dict[str, dict] = {}
    for row in payload:
        if not isinstance(row, dict):
            continue
        assembly_type = str(row.get("assembly_type", "")).strip()
        if assembly_type and assembly_type != MAG_ASSEMBLY_TYPE:
            continue
        primary = str(row.get("sample_accession", "") or "").strip().upper()
        secondary = str(
            row.get("secondary_sample_accession", "") or ""
        ).strip().upper()
        # ERS/SRS/DRS is the native ENA object identifier used for MODIFY.
        update_accession = secondary if re.fullmatch(r"[ESD]RS\d+", secondary) else primary
        if not _SAMPLE_ACCESSION_RE.fullmatch(update_accession):
            continue
        item = discovered.setdefault(
            update_accession,
            {
                "sample_accession": update_accession,
                "biosample_accession": primary if primary != update_accession else "",
                "assembly_accessions": [],
                "assembly_names": [],
                "scientific_name": str(row.get("scientific_name", "") or ""),
                "tax_id": str(row.get("tax_id", "") or ""),
                "study_accession": str(
                    row.get("study_accession", "") or canonical_accession
                ),
            },
        )
        assembly_accession = str(row.get("assembly_accession", "") or "").strip()
        assembly_name = str(row.get("assembly_name", "") or "").strip()
        if assembly_accession and assembly_accession not in item["assembly_accessions"]:
            item["assembly_accessions"].append(assembly_accession)
        if assembly_name and assembly_name not in item["assembly_names"]:
            item["assembly_names"].append(assembly_name)

    return sorted(discovered.values(), key=lambda item: item["sample_accession"])


def sample_candidates_from_accessions(accessions: str | Iterable[str]) -> list[dict]:
    return [
        {
            "sample_accession": accession,
            "biosample_accession": "",
            "assembly_accessions": [],
            "assembly_names": [],
            "scientific_name": "",
            "tax_id": "",
            "study_accession": "",
        }
        for accession in parse_sample_accessions(accessions)
    ]


# Authenticated sample retrieval
def _dropbox_base(portal: str) -> str:
    host = "wwwdev.ebi.ac.uk" if str(portal).casefold() == "testing" else "www.ebi.ac.uk"
    return f"https://{host}/ena/submit/drop-box"


def _webin_v2_base(portal: str) -> str:
    host = "wwwdev.ebi.ac.uk" if str(portal).casefold() == "testing" else "www.ebi.ac.uk"
    return f"https://{host}/ena/submit/webin-v2"


def _reports_base(portal: str) -> str:
    host = "wwwdev.ebi.ac.uk" if str(portal).casefold() == "testing" else "www.ebi.ac.uk"
    return f"https://{host}/ena/submit/report"


def verify_sample_account_access(
    sample_accession: str,
    *,
    portal: str,
    username: str,
    password: str,
    client=requests,
    timeout: int = 30,
) -> None:
    """Verify that a sample is listed for the authenticated Webin account.

    The drop-box XML endpoint also returns public records owned by other
    submitters, so successful retrieval alone does not imply MODIFY access.
    ENA's authenticated Reports Service is account-scoped and is therefore
    checked before a record is offered for editing.
    """
    accessions = parse_sample_accessions([sample_accession])
    if not accessions:
        raise ValueError("A sample accession is required")
    accession = accessions[0]
    url = f"{_reports_base(portal)}/samples/{quote(accession)}"
    try:
        response = client.get(
            url,
            params={"format": "json", "max-results": 1},
            auth=(username, password),
            timeout=timeout,
        )
        status_code = getattr(response, "status_code", 200)
        if status_code in {401, 403}:
            raise EnaSampleUpdateError(
                "ENA rejected the Webin credentials for the selected portal"
            )
        if status_code == 404:
            raise EnaSampleUpdateError(
                f"{accession} is not available in this Webin account's Samples Report "
                "on the selected portal. Check whether the sample belongs to Production "
                "rather than Testing (or vice versa), then verify that the account owns it."
            )
        response.raise_for_status()
    except EnaSampleUpdateError:
        raise
    except requests.RequestException as exc:
        raise EnaSampleUpdateError(
            f"Could not verify update access for {accession} through ENA's "
            f"authenticated Reports Service: {exc}"
        ) from exc

    body = str(getattr(response, "text", "") or "").strip()
    if not body and getattr(response, "content", b""):
        body = bytes(response.content).decode("utf-8", errors="replace").strip()
    if body.casefold() in {"", "[]", "{}", "null"}:
        raise EnaSampleUpdateError(
            f"{accession} is not available in this Webin account's Samples Report on "
            "the selected portal. Check whether the sample belongs to Production rather "
            "than Testing (or vice versa), then verify that the account owns it."
        )


def _local_name(element: ET.Element) -> str:
    return element.tag.rsplit("}", 1)[-1]


def _find_first(element: ET.Element, name: str) -> ET.Element | None:
    for child in element.iter():
        if _local_name(child) == name:
            return child
    return None


def _find_direct(element: ET.Element, name: str) -> ET.Element | None:
    for child in list(element):
        if _local_name(child) == name:
            return child
    return None


def _find_sample(root: ET.Element) -> ET.Element | None:
    return root if _local_name(root) == "SAMPLE" else _find_first(root, "SAMPLE")


def canonical_sample_hash(sample_xml: str | bytes) -> str:
    """Return a whitespace-insensitive hash for optimistic concurrency checks."""
    content = sample_xml.decode("utf-8") if isinstance(sample_xml, bytes) else sample_xml
    try:
        canonical = ET.canonicalize(content, strip_text=True)
    except (ET.ParseError, TypeError, ValueError):
        canonical = content.strip()
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


# Sample XML parsing and editing
def parse_sample_document(content: str | bytes, *, requested_accession: str = "") -> dict:
    raw = content.decode("utf-8") if isinstance(content, bytes) else str(content)
    try:
        root = ET.fromstring(raw)
    except ET.ParseError as exc:
        raise EnaSampleUpdateError(f"ENA returned invalid sample XML: {exc}") from exc
    sample = _find_sample(root)
    if sample is None:
        raise EnaSampleUpdateError("ENA response does not contain a SAMPLE object")

    sample_xml = ET.tostring(sample, encoding="unicode")
    accession = str(sample.get("accession") or requested_accession).strip().upper()
    alias = str(sample.get("alias") or "").strip()
    center_name = str(sample.get("center_name") or "").strip()
    title_element = _find_direct(sample, "TITLE")
    sample_name_element = _find_direct(sample, "SAMPLE_NAME")
    tax_id_element = (
        _find_direct(sample_name_element, "TAXON_ID") if sample_name_element is not None else None
    )
    scientific_element = (
        _find_direct(sample_name_element, "SCIENTIFIC_NAME")
        if sample_name_element is not None
        else None
    )

    attributes: dict[str, list[dict]] = {}
    attribute_parent = _find_direct(sample, "SAMPLE_ATTRIBUTES")
    if attribute_parent is not None:
        for attribute in list(attribute_parent):
            if _local_name(attribute) != "SAMPLE_ATTRIBUTE":
                continue
            tag_element = _find_direct(attribute, "TAG")
            value_element = _find_direct(attribute, "VALUE")
            units_element = _find_direct(attribute, "UNITS")
            tag = str(tag_element.text or "").strip() if tag_element is not None else ""
            if not tag:
                continue
            attributes.setdefault(tag, []).append(
                {
                    "value": str(value_element.text or "").strip()
                    if value_element is not None
                    else "",
                    "units": str(units_element.text or "").strip()
                    if units_element is not None
                    else "",
                }
            )

    conflicts = {
        tag: sorted({item["value"] for item in values})
        for tag, values in attributes.items()
        if len({item["value"] for item in values}) > 1
    }
    return {
        "requested_accession": str(requested_accession or accession).strip().upper(),
        "accession": accession,
        "alias": alias,
        "center_name": center_name,
        "title": str(title_element.text or "").strip() if title_element is not None else "",
        "tax_id": str(tax_id_element.text or "").strip()
        if tax_id_element is not None
        else "",
        "scientific_name": str(scientific_element.text or "").strip()
        if scientific_element is not None
        else "",
        "attributes": attributes,
        "conflicting_duplicate_tags": conflicts,
        "sample_xml": sample_xml,
        "canonical_hash": canonical_sample_hash(sample_xml),
    }


def fetch_sample_document(
    sample_accession: str,
    *,
    portal: str,
    username: str,
    password: str,
    client=requests,
    timeout: int = 30,
) -> dict:
    accessions = parse_sample_accessions([sample_accession])
    if not accessions:
        raise ValueError("A sample accession is required")
    accession = accessions[0]
    verify_sample_account_access(
        accession,
        portal=portal,
        username=username,
        password=password,
        client=client,
        timeout=timeout,
    )
    url = f"{_dropbox_base(portal)}/samples/{quote(accession)}"
    try:
        response = client.get(url, auth=(username, password), timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise EnaSampleUpdateError(
            f"Could not retrieve {accession} from the authenticated ENA service: {exc}"
        ) from exc
    document = parse_sample_document(response.content, requested_accession=accession)
    document["source_url"] = getattr(response, "url", url)
    return document


def document_to_row(document: dict, field_names: Iterable[str]) -> dict:
    row = {field: "" for field in field_names}
    row["sample_accession"] = document.get("accession", "")
    row["sample_name"] = document.get("alias", "")
    row["organism"] = document.get("scientific_name", "")
    row["tax_id"] = document.get("tax_id", "")
    for tag, values in document.get("attributes", {}).items():
        if tag in row and values:
            row[tag] = values[0].get("value", "")
    return row


def _is_empty(value) -> bool:
    return str(value if value is not None else "").strip().casefold() in _EMPTY_VALUES


def _set_direct_text(parent: ET.Element, name: str, value: str) -> None:
    child = _find_direct(parent, name)
    if child is None:
        child = ET.SubElement(parent, name)
    child.text = str(value)


def apply_row_to_sample(document: dict, row: dict, editable_fields: Iterable[str]) -> ET.Element:
    """Apply one edited row while preserving immutable and unknown XML content."""
    sample = ET.fromstring(document["sample_xml"])
    if str(sample.get("accession") or document.get("accession", "")).upper() != str(
        document.get("accession", "")
    ).upper():
        raise EnaSampleUpdateError("Sample accession changed while preparing the update")

    sample_name_element = _find_direct(sample, "SAMPLE_NAME")
    proposed_tax_id = str(row.get("tax_id", "") or "").strip()
    proposed_scientific_name = str(row.get("organism", "") or "").strip()
    if (
        proposed_tax_id != str(document.get("tax_id", "") or "").strip()
        or proposed_scientific_name
        != str(document.get("scientific_name", "") or "").strip()
    ):
        if sample_name_element is None:
            sample_name_element = ET.SubElement(sample, "SAMPLE_NAME")
        if proposed_tax_id != str(document.get("tax_id", "") or "").strip():
            _set_direct_text(sample_name_element, "TAXON_ID", proposed_tax_id)
        if (
            proposed_scientific_name
            != str(document.get("scientific_name", "") or "").strip()
        ):
            scientific_element = _find_direct(sample_name_element, "SCIENTIFIC_NAME")
            if proposed_scientific_name:
                _set_direct_text(
                    sample_name_element, "SCIENTIFIC_NAME", proposed_scientific_name
                )
            elif scientific_element is not None:
                sample_name_element.remove(scientific_element)

    attributes_parent = _find_direct(sample, "SAMPLE_ATTRIBUTES")
    if attributes_parent is None:
        attributes_parent = ET.SubElement(sample, "SAMPLE_ATTRIBUTES")

    for field in editable_fields:
        if field in {"sample_accession", "sample_name", "organism", "tax_id"}:
            continue
        new_value = str(row.get(field, "") or "").strip()
        original_values = document.get("attributes", {}).get(field, [])
        original_value = (
            str(original_values[0].get("value", "") or "").strip()
            if original_values
            else ""
        )
        # Only touch edited values. Everything else should remain as ENA returned it.
        if new_value == original_value:
            continue
        matches = []
        for attribute in list(attributes_parent):
            if _local_name(attribute) != "SAMPLE_ATTRIBUTE":
                continue
            tag_element = _find_direct(attribute, "TAG")
            if tag_element is not None and str(tag_element.text or "").strip() == field:
                matches.append(attribute)

        if _is_empty(new_value):
            for attribute in matches:
                attributes_parent.remove(attribute)
            continue

        if matches:
            # Preserve repeated identical attributes and their units. Conflicting
            # duplicates are blocked by the UI and worker before this point.
            for attribute in matches:
                _set_direct_text(attribute, "VALUE", new_value)
        else:
            attribute = ET.SubElement(attributes_parent, "SAMPLE_ATTRIBUTE")
            _set_direct_text(attribute, "TAG", field)
            _set_direct_text(attribute, "VALUE", new_value)
            units = _UNITS_FOR_NEW_ATTRIBUTES.get(field)
            if units:
                _set_direct_text(attribute, "UNITS", units)
    return sample


# ENA update documents and receipts
def _submission_actions(validate_only: bool) -> ET.Element:
    submission_set = ET.Element("SUBMISSION_SET")
    submission = ET.SubElement(submission_set, "SUBMISSION")
    actions = ET.SubElement(submission, "ACTIONS")
    action = ET.SubElement(actions, "ACTION")
    ET.SubElement(action, "MODIFY")
    if validate_only:
        validation_action = ET.SubElement(actions, "ACTION")
        ET.SubElement(validation_action, "VALIDATE")
    return submission_set


def build_modify_xml(
    documents: Iterable[dict],
    rows: Iterable[dict],
    editable_fields: Iterable[str],
    *,
    validate_only: bool = False,
) -> bytes:
    documents = list(documents)
    rows = list(rows)
    editable_fields = list(editable_fields)
    if not documents or not rows:
        raise EnaSampleUpdateError(
            "The metadata update is empty; load at least one sample and change a value"
        )
    if not editable_fields:
        raise EnaSampleUpdateError("No editable metadata fields were provided")
    documents_by_accession = {
        str(document.get("accession", "")).upper(): document for document in documents
    }
    root = ET.Element("WEBIN")
    root.append(_submission_actions(validate_only))
    sample_set = ET.SubElement(root, "SAMPLE_SET")
    used = set()
    changed = False
    for row in rows:
        accession = str(row.get("sample_accession", "")).strip().upper()
        document = documents_by_accession.get(accession)
        if not document:
            raise EnaSampleUpdateError(f"No authoritative XML loaded for {accession}")
        if accession in used:
            raise EnaSampleUpdateError(f"Duplicate sample accession in update: {accession}")
        used.add(accession)
        updated_sample = apply_row_to_sample(document, row, editable_fields)
        updated_xml = ET.tostring(updated_sample, encoding="unicode")
        original_hash = document.get("canonical_hash") or canonical_sample_hash(
            document["sample_xml"]
        )
        if canonical_sample_hash(updated_xml) != original_hash:
            changed = True
        sample_set.append(updated_sample)
    if not changed:
        raise EnaSampleUpdateError("The metadata update contains no changed values")
    return ET.tostring(root, encoding="UTF-8", xml_declaration=True)


def build_rollback_xml(
    documents: Iterable[dict], *, validate_only: bool = False
) -> bytes:
    """Build a MODIFY document that restores the authoritative sample XML."""
    documents = list(documents)
    if not documents:
        raise EnaSampleUpdateError("The metadata rollback is empty")

    root = ET.Element("WEBIN")
    root.append(_submission_actions(validate_only))
    sample_set = ET.SubElement(root, "SAMPLE_SET")
    used = set()
    for document in documents:
        accession = str(document.get("accession", "") or "").strip().upper()
        if not accession:
            raise EnaSampleUpdateError("A rollback sample has no accession")
        if accession in used:
            raise EnaSampleUpdateError(
                f"Duplicate sample accession in rollback: {accession}"
            )
        used.add(accession)
        try:
            sample = ET.fromstring(document["sample_xml"])
        except (KeyError, ET.ParseError) as exc:
            raise EnaSampleUpdateError(
                f"Could not restore the authoritative XML for {accession}"
            ) from exc
        sample_set.append(sample)

    return ET.tostring(root, encoding="UTF-8", xml_declaration=True)


def build_sample_set_xml(documents: Iterable[dict]) -> bytes:
    root = ET.Element("SAMPLE_SET")
    for document in documents:
        root.append(ET.fromstring(document["sample_xml"]))
    return ET.tostring(root, encoding="UTF-8", xml_declaration=True)


def rows_diff(before_rows: Iterable[dict], after_rows: Iterable[dict]) -> list[dict]:
    before = {
        str(row.get("sample_accession", "")).upper(): row for row in before_rows
    }
    changes = []
    for row in after_rows:
        accession = str(row.get("sample_accession", "")).upper()
        old = before.get(accession, {})
        columns = list(dict.fromkeys([*old.keys(), *row.keys()]))
        for field in columns:
            if field in {"sample_accession", "sample_name"}:
                continue
            previous = str(old.get(field, "") or "").strip()
            proposed = str(row.get(field, "") or "").strip()
            if previous == proposed:
                continue
            changes.append(
                {
                    "sample_accession": accession,
                    "sample_name": str(row.get("sample_name", "") or ""),
                    "field": field,
                    "current_ena_value": previous,
                    "proposed_value": proposed,
                }
            )
    return changes


def receipt_success(receipt_xml: str | bytes) -> tuple[bool, list[str]]:
    content = receipt_xml.decode("utf-8") if isinstance(receipt_xml, bytes) else receipt_xml
    try:
        root = ET.fromstring(content)
    except ET.ParseError as exc:
        raise EnaSampleUpdateError(f"ENA returned an invalid receipt: {exc}") from exc
    messages = [
        str(element.text or "").strip()
        for element in root.iter()
        if _local_name(element) in {"ERROR", "WARNING", "INFO"}
        and str(element.text or "").strip()
    ]
    return str(root.get("success", "")).casefold() == "true", messages


def submit_modify_xml(
    xml_bytes: bytes,
    *,
    portal: str,
    username: str,
    password: str,
    asynchronous: bool,
    client=requests,
    timeout: int = 120,
    poll_timeout: int = 1800,
) -> str:
    """Submit a Webin MODIFY document and return its receipt XML."""
    if not bytes(xml_bytes or b"").strip():
        raise EnaSampleUpdateError("The ENA metadata-update document is empty")
    try:
        root = ET.fromstring(xml_bytes)
    except ET.ParseError as exc:
        raise EnaSampleUpdateError(
            f"The ENA metadata-update document is invalid XML: {exc}"
        ) from exc
    element_names = {_local_name(element) for element in root.iter()}
    if "MODIFY" not in element_names or "SAMPLE" not in element_names:
        raise EnaSampleUpdateError(
            "The ENA metadata-update document must contain MODIFY and at least one SAMPLE"
        )
    if not str(username or "").strip() or not str(password or "").strip():
        raise EnaSampleUpdateError(
            "Both the Webin username and password are required for an ENA metadata update"
        )
    path = "/submit/queue" if asynchronous else "/submit"
    url = f"{_webin_v2_base(portal)}{path}"
    # The synchronous route returns the requested receipt XML. The queue route
    # first returns a JSON document containing its poll link, even though the
    # submitted document itself remains XML.
    headers = {
        "Accept": "application/json" if asynchronous else "application/xml",
        "Content-Type": "application/xml",
    }
    auth = (username, password)
    try:
        response = client.post(
            url,
            data=xml_bytes,
            headers=headers,
            auth=auth,
            timeout=timeout,
        )
        if not response.ok:
            detail = " ".join(str(response.text or "").split())[:1000]
            message = f"ENA returned HTTP {response.status_code} for {path}"
            if detail:
                message += f": {detail}"
            raise EnaSampleUpdateError(message)
        if not asynchronous:
            return response.text

        try:
            payload = response.json()
            poll_url = payload["_links"]["poll"]["href"]
        except (ValueError, KeyError, TypeError) as exc:
            raise EnaSampleUpdateError(
                "ENA did not return a poll URL for the queued update"
            ) from exc

        deadline = time.monotonic() + poll_timeout
        while time.monotonic() < deadline:
            polled = client.get(
                poll_url,
                headers={"Accept": "application/xml"},
                auth=auth,
                timeout=timeout,
            )
            if polled.status_code == 200:
                return polled.text
            if polled.status_code not in {202, 404}:
                polled.raise_for_status()
            time.sleep(3)
        raise EnaSampleUpdateError("Timed out waiting for the ENA metadata-update receipt")
    except requests.RequestException as exc:
        raise EnaSampleUpdateError(f"ENA metadata update request failed: {exc}") from exc
