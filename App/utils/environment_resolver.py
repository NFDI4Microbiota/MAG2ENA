"""Resolve MIxS environmental fields against ENVO/OLS for curator review."""

from __future__ import annotations

from datetime import datetime, timezone
from difflib import SequenceMatcher
import re
import time
from typing import Iterable

import pandas as pd
import requests


ENVIRONMENT_FIELDS = (
    "broad-scale environmental context",
    "local environmental context",
    "environmental medium",
)

ENVO_ROOTS = {
    "broad-scale environmental context": {
        "obo_id": "ENVO:00000428",
        "iri": "http://purl.obolibrary.org/obo/ENVO_00000428",
        "name": "biome",
    },
    "local environmental context": {
        "obo_id": "ENVO:00002297",
        # Current ENVO does not place all valid smaller-grain entities (for
        # example, lake) below its historical environmental-feature class.
        # Local-context searches are therefore ENVO-wide, with biome/material
        # candidates detected and penalised below rather than silently omitted.
        "iri": None,
        "name": "local-context ENVO",
    },
    "environmental medium": {
        "obo_id": "ENVO:00010483",
        "iri": "http://purl.obolibrary.org/obo/ENVO_00010483",
        "name": "environmental material",
    },
}

_EMPTY_VALUES = {"", "<na>", "none", "nan"}
_SAMPLE_ACCESSION_RE = re.compile(
    r"(?:[ESD]RS\d{6,}|SAM[END][AG]?\d+|EGAN\d{11})",
    re.IGNORECASE,
)
_CURIE_RE = re.compile(r"^(?P<prefix>[A-Za-z][A-Za-z0-9_]+):(?P<local>[A-Za-z0-9_]+)$")
_STRUCTURED_TERM_RE = re.compile(
    r"^(?P<label>.+?)\s*\[(?P<curie>[A-Za-z][A-Za-z0-9_]+:[A-Za-z0-9_]+)\]\s*$"
)
_ENVO_MENTION_RE = re.compile(
    r"(?P<label>[A-Za-z][A-Za-z0-9 /_\[\]-]{1,120}?)\s*"
    r"(?:\((?P<paren>ENVO:\d+)\)|\[(?P<bracket>ENVO:\d+)\])",
    re.IGNORECASE,
)

_LOCAL_ROLE_WORDS = {
    "aquifer",
    "borehole",
    "cave",
    "channel",
    "estuary",
    "floor",
    "habitat",
    "layer",
    "plume",
    "reef",
    "shelf",
    "shore",
    "site",
    "vent",
    "zone",
}
_MATERIAL_ROLE_WORDS = {
    "air",
    "ash",
    "biofilm",
    "brine",
    "feces",
    "faeces",
    "fluid",
    "ice",
    "mat",
    "mud",
    "sediment",
    "sludge",
    "soil",
    "stool",
    "water",
}
_PROCESS_WORDS = {
    "formation process",
    "material transformation process",
    "rejection",
    "process",
}


class ResolverServiceError(RuntimeError):
    """Raised when OLS cannot answer a resolver request."""


# Parse existing values and parent-sample evidence
def is_empty(value) -> bool:
    if pd.isna(value):
        return True
    return str(value).strip().lower() in _EMPTY_VALUES


def extract_sample_accession(value) -> str:
    """Return one normalized INSDC sample accession, or an empty string.

    Parent-record enrichment requires a single sample record. Comma-separated
    lists, ranges, run accessions, and arbitrary text are therefore left for
    table-only evidence instead of silently selecting the wrong parent.
    """
    if is_empty(value):
        return ""
    match = _SAMPLE_ACCESSION_RE.fullmatch(str(value).strip())
    return match.group(0).upper() if match else ""


def parse_term(value):
    """Parse ``label [CURIE]`` or a bare CURIE; return None for free text."""
    if is_empty(value):
        return None
    text = str(value).strip()
    structured = _STRUCTURED_TERM_RE.fullmatch(text)
    if structured:
        return {
            "label": structured.group("label").strip(),
            "curie": structured.group("curie").upper(),
        }
    if _CURIE_RE.fullmatch(text):
        return {"label": "", "curie": text.upper()}
    return None


def format_term(label: str, curie: str) -> str:
    """Return the canonical MIxS identified-term representation."""
    return f"{str(label).strip()} [{str(curie).strip().upper()}]"


def _mention_label(raw_label: str) -> str:
    """Reduce prose before an ontology identifier to a readable term label."""
    label = " ".join(str(raw_label).split()).strip(" ,;:()[]")
    # Descriptions commonly use constructions such as "targeting a [SRF]
    # surface water layer (ENVO:...)".  Retain the part after the role cue.
    cue_parts = re.split(
        r"(?i)\b(?:targeting|within|including|collected in|consists? of|composed of)\b",
        label,
    )
    label = cue_parts[-1].strip(" ,;:()[]")
    label = re.sub(r"^\[[^]]+\]\s*", "", label)
    label = re.sub(r"^(?:in|within)\s+(?:an?|the)\s+", "", label, flags=re.IGNORECASE)
    label = re.sub(r"^(?:an?|the)\s+", "", label, flags=re.IGNORECASE)
    words = label.split()
    if len(words) > 10:
        label = " ".join(words[-10:])
    return label


def _field_for_parent_mention(label: str, context: str) -> str | None:
    """Infer a MIxS field only when a parent's prose supplies a clear role."""
    label_norm = _normalise_text(label)
    context_norm = _normalise_text(context)
    label_tokens = set(label_norm.split())

    if "biome" in label_tokens:
        return "broad-scale environmental context"

    # Spatial entities take precedence over material words: "surface water
    # layer" is a local feature, not an environmental-medium value.
    if label_tokens & _LOCAL_ROLE_WORDS:
        return "local environmental context"

    if label_tokens & _MATERIAL_ROLE_WORDS:
        return "environmental medium"

    if re.search(
        r"\b(?:broad scale environmental context|environment(?:al)? biome)\b",
        context_norm,
    ):
        return "broad-scale environmental context"
    if re.search(
        r"\b(?:local environmental context|environment(?:al)? feature|targeting)\b",
        context_norm,
    ):
        return "local environmental context"
    if re.search(
        r"\b(?:environmental medium|environment(?:al)? material|sample material)\b",
        context_norm,
    ):
        return "environmental medium"
    return None


def extract_parent_envo_evidence(title: str = "", description: str = "") -> dict[str, list[dict]]:
    """Extract explicit, field-assigned ENVO mentions from parent sample text.

    A term is emitted only when its label or surrounding prose identifies one
    of the three MIxS environmental roles.  This deliberately abstains on an
    unqualified identifier rather than offering the same term to every field.
    """
    evidence = {field: [] for field in ENVIRONMENT_FIELDS}
    seen: set[tuple[str, str]] = set()
    for source, text in (
        ("parent sample title ENVO term", title),
        ("parent sample description ENVO term", description),
    ):
        compact = " ".join(str(text or "").split())
        for match in _ENVO_MENTION_RE.finditer(compact):
            curie = (match.group("paren") or match.group("bracket")).upper()
            label = _mention_label(match.group("label"))
            start = max(0, match.start() - 100)
            end = min(len(compact), match.end() + 60)
            field = _field_for_parent_mention(label, compact[start:end])
            if not field or (field, curie) in seen:
                continue
            seen.add((field, curie))
            evidence[field].append(
                {
                    "column": source,
                    "value": format_term(label or curie, curie),
                    "weight": 1.0,
                }
            )
    return evidence


def _normalise_text(value: str) -> str:
    text = str(value).casefold().replace("_", " ")
    text = re.sub(r"\[[^]]+\]", "", text)
    text = re.sub(r"\bmetagenom(?:e|ic)\b", "", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return " ".join(text.split())


def _candidate_synonyms(candidate: dict) -> list[str]:
    raw = candidate.get("synonym") or candidate.get("synonyms") or []
    if isinstance(raw, str):
        return [raw]
    return [str(value) for value in raw if value]


def lexical_score(query: str, candidate: dict) -> tuple[float, str]:
    """Score an OLS candidate and explain the strongest lexical relationship."""
    query_norm = _normalise_text(query)
    if not query_norm:
        return 0.0, "No searchable text"

    label = str(candidate.get("label", ""))
    alternatives = [(label, "preferred label")]
    alternatives.extend((synonym, "synonym") for synonym in _candidate_synonyms(candidate))

    best_score, best_reason = 0.0, "Weak lexical similarity"
    for text, kind in alternatives:
        candidate_norm = _normalise_text(text)
        if not candidate_norm:
            continue
        if candidate_norm == query_norm:
            score = 0.99 if kind == "preferred label" else 0.96
            reason = f"Exact {kind} match"
        elif candidate_norm in query_norm or query_norm in candidate_norm:
            coverage = min(len(candidate_norm), len(query_norm)) / max(
                len(candidate_norm), len(query_norm)
            )
            score = 0.76 + (0.14 * coverage)
            reason = f"Phrase match in {kind}"
        else:
            ratio = SequenceMatcher(None, query_norm, candidate_norm).ratio()
            query_tokens = set(query_norm.split())
            candidate_tokens = set(candidate_norm.split())
            union = query_tokens | candidate_tokens
            token_score = len(query_tokens & candidate_tokens) / len(union) if union else 0.0
            score = min(0.89, (0.62 * ratio) + (0.38 * token_score))
            reason = f"Fuzzy {kind} match"
        if score > best_score:
            best_score, best_reason = score, reason
    return round(best_score, 4), best_reason


def _semantic_role_adjustment(query: str, candidate: dict, field: str) -> tuple[float, str]:
    """Penalise lexical matches that have the wrong MIxS semantic role.

    ENVO is polyhierarchical: a term can occur below a permitted ontology root
    while still being an unsuitable answer for the MIxS question.  For example,
    a tidal *zone* is not the material immediately surrounding a sample, and
    marine *mud* is not a broad environmental system.  These safeguards affect
    confidence rather than silently deleting candidates from curator review.
    """
    label_norm = _normalise_text(candidate.get("label", ""))
    description_norm = _normalise_text(candidate.get("description", ""))
    query_norm = _normalise_text(query)
    label_tokens = set(label_norm.split())
    query_tokens = set(query_norm.split())

    if field == "broad-scale environmental context":
        if "biome" not in label_tokens and label_tokens & _MATERIAL_ROLE_WORDS:
            return 0.35, "material-like term is not a broad environmental system"

    elif field == "environmental medium":
        if label_tokens & _LOCAL_ROLE_WORDS:
            return 0.35, "spatial feature/zone is not an environmental material"
        if any(marker in label_norm or marker in description_norm for marker in _PROCESS_WORDS):
            return 0.35, "process term is not an environmental material"
        generic_adjectives = {"aquatic", "marine", "oceanic", "volcanic"}
        if len(query_tokens) == 1 and query_tokens & generic_adjectives:
            return 0.55, "generic context does not support the proposed material specificity"

    elif field == "local environmental context":
        if any(marker in label_norm or marker in description_norm for marker in _PROCESS_WORDS):
            return 0.35, "process term is not a local environmental entity"
        generic_context = {"aquatic", "marine", "oceanic", "volcanic"}
        if len(query_tokens) == 1 and query_tokens & generic_context and len(label_tokens) > 1:
            return 0.55, "generic context does not support the proposed local entity"

    return 1.0, ""


# OLS queries are cached during each resolver run
class OlsEnvoClient:
    """Small resilient client for field-constrained ENVO searches in OLS4."""

    base_url = "https://www.ebi.ac.uk/ols4/api"

    def __init__(self, *, timeout: float = 8, retries: int = 1):
        self.timeout = timeout
        self.retries = retries
        self.session = requests.Session()
        self._cache: dict[tuple[str, str, bool], list[dict]] = {}
        self._failures: dict[tuple[str, str, bool], str] = {}
        self.ontology_version = self._fetch_ontology_version()

    def _get(self, path: str, *, params: dict | None = None) -> requests.Response:
        last_error: Exception | None = None
        for attempt in range(self.retries + 1):
            try:
                response = self.session.get(
                    f"{self.base_url}{path}",
                    params=params,
                    timeout=self.timeout,
                )
                response.raise_for_status()
                return response
            except requests.RequestException as exc:
                last_error = exc
                if attempt < self.retries:
                    time.sleep(0.25 * (2 ** attempt))
        raise ResolverServiceError(f"OLS request failed: {last_error}")

    def _fetch_ontology_version(self) -> str:
        try:
            payload = self._get("/ontologies/envo").json()
        except (ResolverServiceError, ValueError):
            return "unknown"
        config = payload.get("config") or {}
        version = (
            config.get("version")
            or payload.get("version")
            or config.get("versionIri")
            or config.get("version_iri")
        )
        return str(version or "unknown")

    def search(self, query: str, field: str, *, exact: bool = False) -> list[dict]:
        if field not in ENVO_ROOTS:
            raise ValueError(f"Unsupported environmental field: {field}")
        query = str(query).strip()
        if not query:
            return []
        cache_key = (field, query.casefold(), exact)
        if cache_key in self._cache:
            return self._cache[cache_key]
        if cache_key in self._failures:
            raise ResolverServiceError(self._failures[cache_key])

        params = {
            "q": query,
            "type": "class",
            "queryFields": "obo_id,label,synonym",
            "fieldList": "iri,label,obo_id,description,synonym,type",
            "local": "true",
            "obsoletes": "false",
            "groupField": "true",
            "rows": 20,
        }
        parsed = parse_term(query)
        if (
            field != "local environmental context"
            or not parsed
            or parsed["curie"].startswith("ENVO:")
        ):
            params["ontology"] = "envo"
        root_iri = ENVO_ROOTS[field]["iri"]
        if root_iri:
            params["allChildrenOf"] = root_iri
        if exact:
            params["exact"] = "true"

        try:
            payload = self._get("/search", params=params).json()
        except ResolverServiceError as exc:
            self._failures[cache_key] = str(exc)
            raise
        except ValueError as exc:
            message = f"OLS returned invalid JSON: {exc}"
            self._failures[cache_key] = message
            raise ResolverServiceError(message) from exc
        docs = payload.get("response", {}).get("docs", [])
        results = []
        for doc in docs:
            label, obo_id = doc.get("label"), doc.get("obo_id")
            if not label or not obo_id:
                continue
            description = doc.get("description") or []
            if isinstance(description, list):
                description = description[0] if description else ""
            results.append(
                {
                    "label": str(label),
                    "obo_id": str(obo_id).upper(),
                    "iri": str(doc.get("iri") or ""),
                    "description": str(description),
                    "synonym": _candidate_synonyms(doc),
                }
            )
        self._cache[cache_key] = results
        return results


def _parent_text_for_row(
    row: pd.Series,
    parent_metadata_by_accession: dict[str, dict] | None,
) -> tuple[str, str]:
    if not parent_metadata_by_accession:
        return "", ""
    accession = str(row.get("sample derived from", "")).strip()
    record = parent_metadata_by_accession.get(accession) or {}
    title = (
        record.get("_sample_title")
        or record.get("sample_title")
        or record.get("title")
        or ""
    )
    description = (
        record.get("_sample_description")
        or record.get("sample_description")
        or record.get("description")
        or ""
    )
    return str(title), str(description)


def _evidence_for_row(
    row: pd.Series,
    field: str,
    parent_metadata_by_accession: dict[str, dict] | None = None,
) -> list[dict]:
    evidence: list[dict] = []

    current = row.get(field)
    if not is_empty(current):
        evidence.append({"column": field, "value": str(current).strip(), "weight": 1.0})

    if field == "broad-scale environmental context":
        supporting = (("metagenomic source", 0.78), ("isolation_source", 0.68))
    elif field == "local environmental context":
        supporting = (("isolation_source", 0.82), ("metagenomic source", 0.62))
    else:
        supporting = (("isolation_source", 0.80), ("metagenomic source", 0.58))

    for column, weight in supporting:
        value = row.get(column)
        if not is_empty(value):
            evidence.append({"column": column, "value": str(value).strip(), "weight": weight})

    title, description = _parent_text_for_row(row, parent_metadata_by_accession)
    parent_evidence = extract_parent_envo_evidence(title, description)
    evidence.extend(parent_evidence[field])
    return evidence


def _evidence_signature(
    row: pd.Series,
    field: str,
    parent_metadata_by_accession: dict[str, dict] | None = None,
) -> tuple:
    return (
        field,
        str(row.get("sample derived from", "")).strip(),
        tuple(
            (item["column"], _normalise_text(item["value"]))
            for item in _evidence_for_row(row, field, parent_metadata_by_accession)
        ),
    )


def build_resolution_groups(
    df: pd.DataFrame,
    *,
    empty_only: bool = True,
    parent_metadata_by_accession: dict[str, dict] | None = None,
) -> list[dict]:
    """Group rows sharing a parent and the same environmental evidence."""
    groups: dict[tuple, dict] = {}
    for idx, row in df.iterrows():
        for field in ENVIRONMENT_FIELDS:
            current = row.get(field)
            if empty_only and not is_empty(current):
                continue
            signature = _evidence_signature(row, field, parent_metadata_by_accession)
            if signature not in groups:
                groups[signature] = {
                    "field": field,
                    "row_indices": [],
                    "sample_names": [],
                    "parent_accession": str(row.get("sample derived from", "")).strip(),
                    "evidence": _evidence_for_row(
                        row,
                        field,
                        parent_metadata_by_accession,
                    ),
                    "original_value": "" if is_empty(current) else str(current).strip(),
                }
            groups[signature]["row_indices"].append(idx)
            groups[signature]["sample_names"].append(str(row.get("sample_name", idx + 1)))
    return list(groups.values())


# Candidate ranking and batch proposals
def _merge_candidates(client, group):
    field = group["field"]
    combined: dict[str, dict] = {}
    errors: list[str] = []

    for evidence in group["evidence"]:
        parsed = parse_term(evidence["value"])
        query = parsed["curie"] if parsed else evidence["value"]
        try:
            hits = client.search(query, field, exact=bool(parsed))
        except ResolverServiceError as exc:
            errors.append(str(exc))
            continue

        excluded_local_ids: set[str] = set()
        if field == "local environmental context":
            # A lexical hit can be a valid ENVO class but still describe the
            # wrong MIxS level.  Detect candidates that OLS classifies as a
            # biome or environmental material and lower their confidence.
            for other_field in (
                "broad-scale environmental context",
                "environmental medium",
            ):
                try:
                    excluded_local_ids.update(
                        str(item.get("obo_id", "")).upper()
                        for item in client.search(query, other_field, exact=bool(parsed))
                    )
                except ResolverServiceError:
                    pass

        for hit in hits:
            curie = str(hit.get("obo_id", "")).upper()
            if not curie:
                continue
            if parsed and curie == parsed["curie"]:
                score, lexical_reason = 1.0, "Existing ontology identifier verified"
            else:
                lexical, lexical_reason = lexical_score(evidence["value"], hit)
                score = lexical * evidence["weight"]
            if field == "local environmental context" and curie in excluded_local_ids:
                score *= 0.6
                lexical_reason += "; candidate also belongs to a biome/material branch"
            role_multiplier, role_reason = _semantic_role_adjustment(
                evidence["value"],
                hit,
                field,
            )
            if role_multiplier < 1.0:
                score *= role_multiplier
                lexical_reason += f"; semantic-role warning: {role_reason}"
            scope_name = ENVO_ROOTS[field]["name"]
            scope_description = (
                f"{scope_name} branch" if ENVO_ROOTS[field]["iri"] else f"{scope_name} scope"
            )
            candidate = {
                **hit,
                "value": format_term(hit["label"], curie),
                "score": round(score, 4),
                "reason": (
                    f"{lexical_reason}; {scope_description}; "
                    f"evidence: {evidence['column']}"
                ),
                "evidence_column": evidence["column"],
                "evidence_value": evidence["value"],
            }
            previous = combined.get(curie)
            if previous is None or candidate["score"] > previous["score"]:
                combined[curie] = candidate

    candidates = sorted(combined.values(), key=lambda item: (-item["score"], item["label"]))
    return candidates, errors


def _confidence(candidates: list[dict], group: dict) -> tuple[str, str]:
    if not candidates:
        return "Unresolved", "No suitable term found in the field-appropriate OLS search"

    top = candidates[0]
    second_score = candidates[1]["score"] if len(candidates) > 1 else 0.0
    margin = top["score"] - second_score
    parsed_original = parse_term(group.get("original_value", ""))

    if (
        parsed_original
        and top["obo_id"] == parsed_original["curie"]
        and top["score"] == 1.0
    ):
        return "Verified", top["reason"]
    if top["score"] >= 0.98 and top["reason"].startswith("Exact preferred label match"):
        return "High", top["reason"]
    if top["score"] >= 0.93 and margin >= 0.08:
        return "High", top["reason"]
    if top["score"] >= 0.70 and margin < 0.08:
        return "Ambiguous", f"Close alternatives; {top['reason']}"
    if top["score"] >= 0.70:
        return "Review", top["reason"]
    return "Low", f"Weak contextual evidence; {top['reason']}"


def resolve_environment_batch(
    df,
    client,
    empty_only=True,
    parent_metadata_by_accession=None,
):
    """Build non-mutating, curator-reviewable environmental term proposals."""
    proposals = []
    groups = build_resolution_groups(
        df,
        empty_only=empty_only,
        parent_metadata_by_accession=parent_metadata_by_accession,
    )
    for number, group in enumerate(groups, start=1):
        candidates, errors = _merge_candidates(client, group)
        confidence, reason = _confidence(candidates, group)
        if not candidates and errors:
            confidence = "Service error"
            reason = errors[0]
        top = candidates[0] if candidates else {}
        evidence_text = "; ".join(
            f"{item['column']}: {item['value']}" for item in group["evidence"]
        ) or "No contextual evidence available"
        proposals.append(
            {
                "proposal_id": f"ENV-{number:04d}",
                "field": group["field"],
                "row_indices": group["row_indices"],
                "sample_names": group["sample_names"],
                "row_count": len(group["row_indices"]),
                "parent_accession": group["parent_accession"],
                "original_value": group["original_value"],
                "evidence": evidence_text,
                "selected_value": top.get("value", ""),
                "selected_label": top.get("label", ""),
                "selected_id": top.get("obo_id", ""),
                "score": top.get("score", 0.0),
                "confidence": confidence,
                "reason": reason,
                "accept": False,
                "candidates": candidates,
                "ontology_version": getattr(client, "ontology_version", "unknown"),
            }
        )
    return proposals


def proposals_to_review_dataframe(proposals: Iterable[dict]) -> pd.DataFrame:
    rows = []
    for proposal in proposals:
        rows.append(
            {
                "accept": bool(proposal.get("accept", False)),
                "proposal_id": proposal["proposal_id"],
                "rows": proposal["row_count"],
                "field": proposal["field"],
                "source evidence": proposal["evidence"],
                "proposed value": proposal.get("selected_value", ""),
                "confidence": proposal["confidence"],
                "score": round(float(proposal.get("score", 0.0)), 3),
                "reason": proposal["reason"],
            }
        )
    return pd.DataFrame(rows)


def apply_environment_proposals(
    df: pd.DataFrame,
    proposals: Iterable[dict],
    accepted_ids: Iterable[str],
    *,
    overwrite: bool = False,
) -> tuple[pd.DataFrame, list[dict]]:
    """Apply accepted proposals and return the modified dataframe and audit rows."""
    result = df.copy()
    accepted = set(accepted_ids)
    audit = []
    for proposal in proposals:
        if proposal["proposal_id"] not in accepted or not proposal.get("selected_value"):
            continue
        changed_rows = []
        for idx in proposal["row_indices"]:
            current = result.at[idx, proposal["field"]]
            if overwrite or is_empty(current):
                result.at[idx, proposal["field"]] = proposal["selected_value"]
                changed_rows.append(idx + 1 if isinstance(idx, int) else str(idx))
        if changed_rows:
            audit.append(
                {
                    "proposal_id": proposal["proposal_id"],
                    "rows": ", ".join(map(str, changed_rows)),
                    "sample_names": ", ".join(map(str, proposal.get("sample_names", []))),
                    "parent_accession": proposal.get("parent_accession", ""),
                    "field": proposal["field"],
                    "original_value": proposal.get("original_value", ""),
                    "applied_value": proposal["selected_value"],
                    "selected_label": proposal.get("selected_label", ""),
                    "selected_id": proposal.get("selected_id", ""),
                    "source_evidence": proposal["evidence"],
                    "confidence": proposal["confidence"],
                    "score": proposal.get("score", 0.0),
                    "reason": proposal["reason"],
                    "ontology_version": proposal.get("ontology_version", "unknown"),
                    "alternatives": " | ".join(
                        item.get("value", "") for item in proposal.get("candidates", [])[1:6]
                    ),
                    "decision": "user approved",
                    "applied_at": datetime.now(timezone.utc).isoformat(),
                }
            )
    return result, audit
