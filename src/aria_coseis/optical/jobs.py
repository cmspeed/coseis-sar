"""Optical (autoRIFT) job JSON helpers."""

from __future__ import annotations

import re
from typing import Any


def get_utm_zone(granule_id: str) -> str:
    """Extracts UTM zone from Sentinel-2 Granule ID (e.g., 'T46QHK' -> '46').
    :param granule_id: The granule ID string from which to extract the UTM zone.
    :return: The UTM zone as a string, or 'Unknown' if it cannot be extracted.
    """
    match = re.search(r"_T(\d{2})[A-Z]{3}_", granule_id)
    return match.group(1) if match else "Unknown"


def make_optical_job_json(
    title: str,
    event_id: str,
    orbit_id: str,
    pre_date: str,
    post_date: str,
    reference_ids: list[str],
    secondary_ids: list[str],
    pre_cc: float | None = None,
    post_cc: float | None = None,
    status: str = "COMPLETE",
    zone_suffix: str | None = None,
) -> dict[str, Any]:
    """
    Helper function to create a JSON object for an AUTORIFT job. Matches the schema defined in ARIA_AUTORIFT.yml.
    :param title: Title of the job (usually the earthquake event name)
    :param event_id: The earthquake event ID (e.g., "us7000dflf")
    :param orbit_id: The orbit number (e.g., "047")
    :param pre_date: Date of the Pre-Event image (YYYY-MM-DD)
    :param post_date: Date of the Post-Event image (YYYY-MM-DD)
    :param reference_ids: List of granule IDs for the reference (Pre-Event) images
    :param secondary_ids: List of granule IDs for the secondary (Post-Event) images
    :param status: Status of the job ("COMPLETE", "PARTIAL", "FAILED") - used for internal tracking
    :param zone_suffix: Optional suffix for the orbit key if using "Orbit_Zone" grouping (e.g., "Z46")
    :return: A dictionary representing the job configuration for AUTORIFT.
    """
    # Create a descriptive job name
    orbit_str = f"R{orbit_id}"
    if zone_suffix:
        orbit_str += f"-{zone_suffix}"

    job_name = f"{title}-S2-{orbit_str}-{pre_date}_{post_date}"

    job_json = {
        "name": job_name,
        "job_type": "AUTORIFT",
        "event_id": event_id,
        "pre_cloud_cover": round(pre_cc, 2)
        if pre_cc is not None
        else None,  # Temporary debug field
        "post_cloud_cover": round(post_cc, 2)
        if post_cc is not None
        else None,  # Temporary debug field
        "job_parameters": {
            "reference": reference_ids,
            "secondary": secondary_ids,
            "chip_size": 24,
            "search_range": 64,
        },
    }

    # Internal tracking for partial jobs (optional)
    if status != "COMPLETE":
        job_json["status_note"] = status

    return job_json
