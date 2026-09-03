"""Retrieve ENA sample metadata with a strict assembly-index fallback."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET

import requests


ENA_BROWSER_XML_URL = "https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
ENA_PORTAL_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"

_SAMPLE_ACCESSION_RE = re.compile(
    r"(?:[ESD]RS\d{6,}|SAM[END][AG]?\d+|EGAN\d{11})",
    re.IGNORECASE,
)

_DATE_REGEX = re.compile(
    r"^[12][0-9]{3}(-(0[1-9]|1[0-2])(-(0[1-9]|[12][0-9]|3[01])"
    r"(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?"
    r"(/[0-9]{4}(-[0-9]{2}(-[0-9]{2}"
    r"(T[0-9]{2}:[0-9]{2}(:[0-9]{2})?Z?([+-][0-9]{1,2})?)?)?)?)?$"
)

_ASSEMBLY_FIELDS = (
    "assembly_accession",
    "assembly_name",
    "sample_accession",
    "secondary_sample_accession",
    "sample_title",
    "description",
    "scientific_name",
    "tax_id",
    "broad_scale_environmental_context",
    "local_environmental_context",
    "environmental_medium",
    "isolation_source",
    "collection_date",
    "country",
    "location",
)

_PORTAL_TO_SAMPLE_ATTRIBUTE = {
    "scientific_name": "organism",
    "tax_id": "tax_id",
    "broad_scale_environmental_context": "env_broad_scale",
    "local_environmental_context": "env_local_scale",
    "environmental_medium": "env_medium",
    "isolation_source": "isolation_source",
    "collection_date": "collection_date",
    "country": "geo_loc_name",
    "location": "lat_lon",
}


def normalize_sample_accession(value: str) -> str:
    """Return a normalized supported sample accession or an empty string."""
    match = _SAMPLE_ACCESSION_RE.fullmatch(str(value or "").strip())
    return match.group(0).upper() if match else ""


def normalize_collection_date(value: str) -> str:
    """Normalize partial ENA collection dates to an accepted ISO 8601 value."""
    text = str(value or "").strip()
    if not text or text.lower() == "missing":
        return ""
    normalized = re.sub(r"^(\d{4})$", r"\1-01-01", text)
    normalized = re.sub(r"^(\d{4})-(\d{2})$", r"\1-\2-01", normalized)
    normalized = re.sub(r"^(\d{4})(\d{2})(\d{2})$", r"\1-\2-\3", normalized)
    return normalized if _DATE_REGEX.fullmatch(normalized) else ""


def split_coordinates(value: str) -> tuple[str, str]:
    """Convert ENA ``lat_lon`` text to signed decimal latitude/longitude."""
    text = str(value or "").strip()
    if not text:
        return "", ""
    match = re.search(
        r"(\d*\.?\d+)\s*([NS])\s+(\d*\.?\d+)\s*([EW])",
        text,
        re.IGNORECASE,
    )
    if match:
        latitude = float(match.group(1)) * (-1 if match.group(2).upper() == "S" else 1)
        longitude = float(match.group(3)) * (-1 if match.group(4).upper() == "W" else 1)
        return str(round(latitude, 6)), str(round(longitude, 6))
    parts = text.split()
    if len(parts) == 2:
        try:
            return str(float(parts[0])), str(float(parts[1]))
        except ValueError:
            pass
    return "", ""


def attrs_to_mag_metadata(attrs: dict) -> dict:
    """Map an ENA parent-sample record to ERC000047-shaped MAG columns."""
    env_broad = attrs.get("env_broad_scale") or attrs.get("env_biome") or ""
    env_local = attrs.get("env_local_scale") or attrs.get("env_feature") or ""
    env_medium = attrs.get("env_medium") or attrs.get("env_material") or ""

    if not env_broad and not env_local and not env_medium:
        gold = str(attrs.get("GOLD Ecosystem Classification", "") or "")
        if gold:
            parts = gold.split(" | ")
            env_broad = parts[1] if len(parts) > 1 else ""
            env_local = parts[2] if len(parts) > 2 else ""
            env_medium = parts[3] if len(parts) > 3 else ""

    location_name = str(attrs.get("geo_loc_name", "") or "")
    country = location_name.split(":")[0].strip() if location_name else ""
    latitude, longitude = split_coordinates(attrs.get("lat_lon", ""))

    return {
        "metagenomic source": attrs.get("organism", ""),
        "isolation_source": attrs.get("isolation_source", ""),
        "collection date": normalize_collection_date(attrs.get("collection_date", "")),
        "geographic location (latitude)": latitude,
        "geographic location (longitude)": longitude,
        "broad-scale environmental context": env_broad,
        "local environmental context": env_local,
        "environmental medium": env_medium,
        "geographic location (country and/or sea)": country,
    }


def parse_sample_xml(content: bytes) -> dict:
    """Parse a single ENA sample record into the attribute mapping used by MAG2ENA."""
    root = ET.fromstring(content)
    sample = root if root.tag.rsplit("}", 1)[-1] == "SAMPLE" else root.find(".//SAMPLE")
    if sample is None:
        raise ValueError("the ENA XML response does not contain a SAMPLE record")

    attrs = {}
    for attr in sample.findall(".//SAMPLE_ATTRIBUTE"):
        tag = (attr.findtext("TAG") or "").strip()
        value = (attr.findtext("VALUE") or "").strip()
        if tag:
            attrs[tag] = value

    attrs["_sample_title"] = (sample.findtext("TITLE") or "").strip()
    attrs["_sample_description"] = (
        (sample.findtext("DESCRIPTION") or "").strip().strip('"')
    )
    scientific_name = (
        sample.findtext(".//SAMPLE_NAME/SCIENTIFIC_NAME") or ""
    ).strip()
    if not attrs.get("organism") and scientific_name:
        attrs["organism"] = scientific_name
    attrs["_metadata_source"] = "ENA Browser sample XML"
    return attrs


def _agreed_value(rows: list[dict], field: str) -> tuple[str, list[str]]:
    """Return a value only when all non-empty indexed values agree."""
    values = []
    for row in rows:
        value = str(row.get(field, "") or "").strip()
        if value and value not in values:
            values.append(value)
    if len(values) == 1:
        return values[0], []
    if len(values) > 1:
        return "", values
    return "", []


def fetch_sample_from_linked_assemblies(
    sample_accession: str,
    *,
    client=requests,
    timeout: int = 12,
) -> dict:
    """Recover agreed indexed metadata from assemblies exactly linked to a sample."""
    accession = normalize_sample_accession(sample_accession)
    if not accession:
        raise ValueError(f"Not a supported INSDC sample accession: {sample_accession}")

    query = (
        f'(sample_accession="{accession}" OR '
        f'secondary_sample_accession="{accession}")'
    )
    response = client.get(
        ENA_PORTAL_SEARCH_URL,
        params={
            "result": "assembly",
            "query": query,
            "fields": ",".join(_ASSEMBLY_FIELDS),
            "format": "json",
            "limit": 100,
        },
        timeout=timeout,
    )
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, list):
        raise ValueError("ENA Portal assembly search returned an unexpected response")

    exact_rows = []
    seen_assemblies = set()
    for row in payload:
        if not isinstance(row, dict):
            continue
        linked_accessions = {
            str(row.get("sample_accession", "") or "").strip().upper(),
            str(row.get("secondary_sample_accession", "") or "").strip().upper(),
        }
        if accession not in linked_accessions:
            continue
        assembly_key = str(row.get("assembly_accession", "") or "").strip()
        dedup_key = assembly_key or repr(sorted(row.items()))
        if dedup_key in seen_assemblies:
            continue
        seen_assemblies.add(dedup_key)
        exact_rows.append(row)

    if not exact_rows:
        raise ValueError(
            f"no indexed assembly is exactly linked to sample {accession}"
        )

    attrs = {}
    conflicts = {}
    for portal_field, sample_attribute in _PORTAL_TO_SAMPLE_ATTRIBUTE.items():
        value, conflicting_values = _agreed_value(exact_rows, portal_field)
        if value:
            attrs[sample_attribute] = value
        elif conflicting_values:
            conflicts[sample_attribute] = conflicting_values

    title, title_conflicts = _agreed_value(exact_rows, "sample_title")
    description, description_conflicts = _agreed_value(exact_rows, "description")
    attrs["_sample_title"] = title
    attrs["_sample_description"] = description
    if title_conflicts:
        conflicts["_sample_title"] = title_conflicts
    if description_conflicts:
        conflicts["_sample_description"] = description_conflicts

    attrs["_metadata_source"] = "ENA Portal linked-assembly fallback"
    attrs["_metadata_source_url"] = getattr(response, "url", ENA_PORTAL_SEARCH_URL)
    attrs["_fallback_assembly_accessions"] = sorted(
        {
            str(row.get("assembly_accession", "") or "").strip()
            for row in exact_rows
            if str(row.get("assembly_accession", "") or "").strip()
        }
    )
    attrs["_fallback_conflicts"] = conflicts
    return attrs


def fetch_ena_sample_metadata(
    sample_accession: str,
    *,
    client=requests,
    timeout: int = 12,
) -> dict:
    """Fetch sample XML, falling back to metadata agreed across linked assemblies."""
    accession = normalize_sample_accession(sample_accession)
    if not accession:
        return {"_error": f"Not a supported INSDC sample accession: {sample_accession}"}

    direct_url = ENA_BROWSER_XML_URL.format(accession=accession)
    direct_error = None
    try:
        response = client.get(direct_url, timeout=timeout)
        response.raise_for_status()
        attrs = parse_sample_xml(response.content)
        attrs["_metadata_source_url"] = getattr(response, "url", direct_url)
        return attrs
    except (requests.RequestException, ET.ParseError, ValueError, OSError) as exc:
        direct_error = str(exc)

    try:
        attrs = fetch_sample_from_linked_assemblies(
            accession,
            client=client,
            timeout=timeout,
        )
        attrs["_direct_sample_error"] = direct_error
        return attrs
    except (requests.RequestException, ValueError, OSError) as exc:
        return {
            "_error": (
                f"Direct ENA sample retrieval failed ({direct_error}); "
                f"linked-assembly fallback failed ({exc})"
            )
        }
