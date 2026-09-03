"""Geographic summaries for the Home and Update metadata editors."""

from __future__ import annotations

import html
import math
import re

import pandas as pd
import pydeck as pdk
import streamlit as st


LATITUDE_FIELD = "geographic location (latitude)"
LONGITUDE_FIELD = "geographic location (longitude)"
PARENT_FIELD = "sample derived from"

_EMPTY_VALUES = {"", "<na>", "nan", "none", "missing"}
_PARENT_ACCESSION = re.compile(
    r"(?:[ESD]RS\d{6,}|SAM[END][AG]?\d+|EGAN\d{11})",
    re.IGNORECASE,
)
_COORDINATE = re.compile(
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?)\s*°?\s*([NSEW])?",
    re.IGNORECASE,
)


def _text(value) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in _EMPTY_VALUES else text


def parse_coordinate(value, axis: str) -> float | None:
    """Parse a decimal or cardinal coordinate and enforce geographic bounds."""
    text = _text(value)
    match = _COORDINATE.fullmatch(text)
    if not match or axis not in {"latitude", "longitude"}:
        return None

    number = float(match.group(1))
    direction = str(match.group(2) or "").upper()
    allowed = {"N", "S"} if axis == "latitude" else {"E", "W"}
    if direction and direction not in allowed:
        return None
    if direction:
        number = abs(number) * (-1 if direction in {"S", "W"} else 1)

    limit = 90 if axis == "latitude" else 180
    if not math.isfinite(number) or not -limit <= number <= limit:
        return None
    return round(number, 6)


def _parent_accession(value) -> str:
    text = _text(value)
    match = _PARENT_ACCESSION.search(text)
    return match.group(0).upper() if match else text


def _unique_values(rows: list[dict], *fields: str) -> list[str]:
    values = []
    for row in rows:
        for field in fields:
            value = _text(row.get(field, ""))
            if value and value not in values:
                values.append(value)
    return values


def _short_list(values: list[str], limit: int = 4) -> str:
    if not values:
        return "Not provided"
    shown = values[:limit]
    text = "; ".join(shown)
    if len(values) > limit:
        text += f"; +{len(values) - limit} more"
    return html.escape(text)


def build_geographic_overview(metadata: pd.DataFrame) -> tuple[list[dict], dict]:
    """Group MAG rows by parent sample and aggregate identical map coordinates."""
    summary = {
        "metadata_rows": len(metadata),
        "parent_samples": 0,
        "mapped_parent_samples": 0,
        "mapped_mags": 0,
        "map_locations": 0,
        "missing_coordinate_parents": 0,
        "invalid_coordinate_parents": 0,
        "conflicting_coordinate_parents": 0,
        "rows_without_parent": 0,
    }
    if metadata.empty or PARENT_FIELD not in metadata.columns:
        summary["rows_without_parent"] = len(metadata)
        return [], summary

    parents: dict[str, list[dict]] = {}
    for row in metadata.fillna("").to_dict(orient="records"):
        parent = _parent_accession(row.get(PARENT_FIELD, ""))
        if not parent:
            summary["rows_without_parent"] += 1
            continue
        parents.setdefault(parent, []).append(row)
    summary["parent_samples"] = len(parents)

    mapped_parents = []
    for parent, rows in parents.items():
        coordinates = set()
        coordinate_text_present = False
        for row in rows:
            raw_latitude = _text(row.get(LATITUDE_FIELD, ""))
            raw_longitude = _text(row.get(LONGITUDE_FIELD, ""))
            coordinate_text_present |= bool(raw_latitude or raw_longitude)
            latitude = parse_coordinate(raw_latitude, "latitude")
            longitude = parse_coordinate(raw_longitude, "longitude")
            if latitude is not None and longitude is not None:
                coordinates.add((latitude, longitude))

        if not coordinates:
            key = (
                "invalid_coordinate_parents"
                if coordinate_text_present
                else "missing_coordinate_parents"
            )
            summary[key] += 1
            continue
        if len(coordinates) > 1:
            summary["conflicting_coordinate_parents"] += 1
            continue

        latitude, longitude = next(iter(coordinates))
        names = _unique_values(rows, "sample_name")
        mag_count = len(names) if names else len(rows)
        mapped_parents.append(
            {
                "parent": parent,
                "latitude": latitude,
                "longitude": longitude,
                "mag_count": mag_count,
                "mag_names": names,
                "metagenomic_sources": _unique_values(rows, "metagenomic source"),
                "countries": _unique_values(
                    rows, "geographic location (country and/or sea)"
                ),
                "isolation_sources": _unique_values(
                    rows, "isolation_source", "isolation source"
                ),
                "collection_dates": _unique_values(rows, "collection date"),
                "broad_contexts": _unique_values(
                    rows, "broad-scale environmental context"
                ),
                "local_contexts": _unique_values(
                    rows, "local environmental context"
                ),
                "environmental_media": _unique_values(rows, "environmental medium"),
            }
        )

    locations: dict[tuple[float, float], dict] = {}
    for parent in mapped_parents:
        key = (parent["latitude"], parent["longitude"])
        location = locations.setdefault(
            key,
            {
                "latitude": parent["latitude"],
                "longitude": parent["longitude"],
                "parents": [],
                "mag_count": 0,
                "mag_names": [],
                "metagenomic_sources": [],
                "countries": [],
                "isolation_sources": [],
                "collection_dates": [],
                "broad_contexts": [],
                "local_contexts": [],
                "environmental_media": [],
            },
        )
        location["parents"].append((parent["parent"], parent["mag_count"]))
        location["mag_count"] += parent["mag_count"]
        for field in (
            "mag_names",
            "metagenomic_sources",
            "countries",
            "isolation_sources",
            "collection_dates",
            "broad_contexts",
            "local_contexts",
            "environmental_media",
        ):
            for value in parent[field]:
                if value not in location[field]:
                    location[field].append(value)

    points = []
    for location in locations.values():
        shown_parents = location["parents"][:8]
        parent_details = "; ".join(
            f"{accession} ({count} MAG{'s' if count != 1 else ''})"
            for accession, count in shown_parents
        )
        if len(location["parents"]) > len(shown_parents):
            parent_details += (
                f"; +{len(location['parents']) - len(shown_parents)} more parent samples"
            )
        points.append(
            {
                "latitude": location["latitude"],
                "longitude": location["longitude"],
                "coordinates": (
                    f"{location['latitude']:.6f}, {location['longitude']:.6f}"
                ),
                "parent_count": len(location["parents"]),
                "parent_details": html.escape(parent_details),
                "mag_count": location["mag_count"],
                "mag_names": _short_list(location["mag_names"], limit=6),
                "metagenomic_sources": _short_list(location["metagenomic_sources"]),
                "countries": _short_list(location["countries"]),
                "isolation_sources": _short_list(location["isolation_sources"]),
                "collection_dates": _short_list(location["collection_dates"]),
                "broad_contexts": _short_list(location["broad_contexts"]),
                "local_contexts": _short_list(location["local_contexts"]),
                "environmental_media": _short_list(location["environmental_media"]),
                "radius": min(34, 12 + 4 * math.sqrt(location["mag_count"])),
            }
        )

    summary["mapped_parent_samples"] = len(mapped_parents)
    summary["mapped_mags"] = sum(item["mag_count"] for item in mapped_parents)
    summary["map_locations"] = len(points)
    return points, summary


def _view_state(points: list[dict]) -> pdk.ViewState:
    latitudes = [point["latitude"] for point in points]
    longitudes = [point["longitude"] for point in points]
    latitude = sum(latitudes) / len(latitudes)

    sin_mean = sum(math.sin(math.radians(value)) for value in longitudes)
    cos_mean = sum(math.cos(math.radians(value)) for value in longitudes)
    longitude = math.degrees(math.atan2(sin_mean, cos_mean))
    longitude_span = max(
        abs((value - longitude + 180) % 360 - 180) for value in longitudes
    ) * 2
    span = max(max(latitudes) - min(latitudes), longitude_span)
    if span < 0.02:
        zoom = 10
    elif span < 0.2:
        zoom = 8
    elif span < 2:
        zoom = 6
    elif span < 12:
        zoom = 4
    elif span < 50:
        zoom = 2
    else:
        zoom = 1
    return pdk.ViewState(latitude=latitude, longitude=longitude, zoom=zoom)


def render_geographic_overview(metadata: pd.DataFrame, *, key: str) -> None:
    """Render a collapsed, read-only geographic summary of the current table."""
    with st.expander("Geographic overview", expanded=False):
        st.caption(
            "One point represents a parent metagenomic sample. Samples with identical "
            "coordinates are aggregated at their exact location; point size reflects the "
            "number of associated MAGs. This view does not alter the metadata table."
        )
        points, summary = build_geographic_overview(metadata)

        metric_columns = st.columns(3)
        metric_columns[0].metric(
            "Mapped parent samples",
            f"{summary['mapped_parent_samples']:,} / {summary['parent_samples']:,}",
        )
        metric_columns[1].metric("Associated MAGs", f"{summary['mapped_mags']:,}")
        metric_columns[2].metric("Map locations", f"{summary['map_locations']:,}")

        notices = []
        if summary["missing_coordinate_parents"]:
            notices.append(
                f"{summary['missing_coordinate_parents']:,} parent sample(s) have no coordinates"
            )
        if summary["invalid_coordinate_parents"]:
            notices.append(
                f"{summary['invalid_coordinate_parents']:,} have incomplete or invalid coordinates"
            )
        if summary["conflicting_coordinate_parents"]:
            notices.append(
                f"{summary['conflicting_coordinate_parents']:,} have conflicting coordinates across MAG rows"
            )
        if summary["rows_without_parent"]:
            notices.append(
                f"{summary['rows_without_parent']:,} MAG row(s) have no parent accession"
            )
        if notices:
            st.warning("Not mapped: " + "; ".join(notices) + ".")
        if not points:
            st.info(
                "No parent samples currently have a complete, valid coordinate pair. "
                "Enter coordinates or use ENA parent-record autofill to populate them."
            )
            return

        layer = pdk.Layer(
            "ScatterplotLayer",
            data=points,
            get_position="[longitude, latitude]",
            get_radius="radius",
            radius_units="pixels",
            radius_min_pixels=8,
            radius_max_pixels=34,
            get_fill_color=[31, 119, 180, 185],
            get_line_color=[12, 67, 105, 230],
            line_width_min_pixels=1,
            pickable=True,
            auto_highlight=True,
            stroked=True,
        )
        tooltip = {
            "html": (
                "<b>Parent sample(s)</b>: {parent_details}<br/>"
                "<b>Associated MAGs</b>: {mag_count}<br/>"
                "<b>Coordinates</b>: {coordinates}<br/>"
                "<b>Country or sea</b>: {countries}<br/>"
                "<b>Metagenomic source</b>: {metagenomic_sources}<br/>"
                "<b>Isolation source</b>: {isolation_sources}<br/>"
                "<b>Collection date</b>: {collection_dates}<br/>"
                "<b>Broad-scale context</b>: {broad_contexts}<br/>"
                "<b>Local context</b>: {local_contexts}<br/>"
                "<b>Environmental medium</b>: {environmental_media}"
            ),
            "style": {
                "backgroundColor": "#17324d",
                "color": "white",
                "fontSize": "12px",
                "maxWidth": "460px",
            },
        }
        deck = pdk.Deck(
            layers=[layer],
            initial_view_state=_view_state(points),
            map_style="https://basemaps.cartocdn.com/gl/positron-gl-style/style.json",
            tooltip=tooltip,
        )
        st.pydeck_chart(deck, height=500, key=key, on_select="ignore")
