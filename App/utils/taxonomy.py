"""ENA taxonomy lookup and MAG-specific preflight rules."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from urllib.parse import quote

import requests


def search_ena_taxonomy(
    query: str,
    *,
    client=requests,
    timeout: int = 8,
) -> list[dict]:
    """Return ENA taxonomy suggestions for a scientific-name query."""
    value = str(query or "").strip()
    if len(value) < 3:
        return []
    try:
        response = client.get(
            "https://www.ebi.ac.uk/ena/taxonomy/rest/suggest-for-search/"
            f"{quote(value)}",
            timeout=timeout,
        )
        response.raise_for_status()
        payload = response.json()
        return payload if isinstance(payload, list) else []
    except Exception:
        return []


def resolve_species_taxon(
    organism: str,
    *,
    searcher: Callable[[str], list[Mapping[str, object]]] = search_ena_taxonomy,
    fetcher: Callable[[str], Mapping[str, object]] | None = None,
) -> dict | None:
    """Select the first exact/preferred suggestion that passes MAG taxonomy policy."""
    if fetcher is None:
        fetcher = fetch_ena_taxon
    query = str(organism or "").strip()
    suggestions = list(searcher(query))
    suggestions.sort(
        key=lambda item: (
            str(item.get("scientificName", "")).casefold() != query.casefold(),
        )
    )
    seen = set()
    for suggestion in suggestions:
        tax_id = str(suggestion.get("taxId", "") or "").strip()
        if not tax_id or tax_id in seen:
            continue
        seen.add(tax_id)
        record = dict(fetcher(tax_id))
        if mag_taxon_validation_error(tax_id, record) is None:
            return {
                "tax_id": tax_id,
                "scientific_name": str(
                    record.get("scientificName")
                    or suggestion.get("scientificName")
                    or query
                ).strip(),
            }
    return None


def fetch_ena_taxon(
    tax_id: str,
    *,
    client=requests,
    timeout: int = 8,
) -> dict:
    """Return an ENA taxonomy record or a structured lookup failure."""
    value = str(tax_id).strip()
    try:
        response = client.get(
            f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{value}",
            timeout=timeout,
        )
        if response.status_code == 404:
            return {"_not_found": True}
        response.raise_for_status()
        payload = response.json()
        return payload if isinstance(payload, dict) else {
            "_error": "ENA taxonomy returned an unexpected response",
        }
    except Exception as exc:
        return {"_error": str(exc)}


def _true(value: object) -> bool:
    return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def is_metagenome_taxon(taxon: Mapping[str, object]) -> bool:
    """Identify taxa Webin-CLI must not accept as the taxonomy of a MAG."""
    if "metagenome" in taxon:
        return _true(taxon.get("metagenome"))

    # Defensive fallbacks for older or incomplete taxonomy responses.
    scientific_name = str(taxon.get("scientificName", "")).strip().lower()
    lineage = str(taxon.get("lineage", "")).lower()
    return (
        scientific_name == "metagenome"
        or scientific_name.endswith(" metagenome")
        or "metagenomes;" in lineage
    )


def mag_taxon_validation_error(
    tax_id: str,
    taxon: Mapping[str, object],
    *,
    require_species: bool = True,
) -> str | None:
    """Return the reason a taxonomy record cannot be used for a MAG."""
    if taxon.get("_not_found"):
        return "Valid NCBI tax ID (not found in ENA taxonomy)"
    if taxon.get("_error"):
        return f"Could not verify tax ID: {taxon['_error']}"

    scientific_name = str(taxon.get("scientificName", "")).strip()
    if is_metagenome_taxon(taxon):
        label = f" ({scientific_name})" if scientific_name else ""
        return (
            "MAG tax ID cannot be a metagenome taxonomy"
            f"{label}; select the organism taxon assigned to the MAG"
        )

    if str(taxon.get("submittable", "true")).strip().lower() == "false":
        return "Tax ID is not submittable in ENA taxonomy"

    rank = str(taxon.get("rank", "unknown")).strip().lower()
    if require_species and rank != "species":
        return f"Tax ID must correspond to a species rank (got: {rank or 'unknown'})"

    return None


def validate_mag_taxonomy_records(
    records: Iterable[Mapping[str, object]],
    *,
    fetcher: Callable[[str], Mapping[str, object]] = fetch_ena_taxon,
    require_species: bool = True,
) -> list[dict]:
    """Validate every MAG taxon and return row-oriented preflight failures."""
    errors = []
    cache: dict[str, Mapping[str, object]] = {}
    for row_number, record in enumerate(records, start=1):
        tax_id = str(record.get("tax_id", "")).strip()
        if not tax_id:
            continue
        if tax_id not in cache:
            cache[tax_id] = fetcher(tax_id)
        taxon = cache[tax_id]
        expected = mag_taxon_validation_error(
            tax_id,
            taxon,
            require_species=require_species,
        )
        if expected:
            errors.append(
                {
                    "row": row_number,
                    "sample_name": str(record.get("sample_name", "")).strip(),
                    "tax_id": tax_id,
                    "expected": expected,
                }
            )
    return errors
