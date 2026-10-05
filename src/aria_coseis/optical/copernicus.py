"""Sentinel-2 scene search and pairing via the Copernicus Data Space OData API."""

from __future__ import annotations

import re
from collections import defaultdict
from datetime import datetime, timezone
from typing import Any

import requests
from shapely import wkt
from shapely.geometry import mapping
from shapely.geometry.base import BaseGeometry

from aria_coseis.config import OPTICAL_CLOUD_THRESHOLD
from aria_coseis.optical.jobs import make_optical_job_json
from aria_coseis.utils import convert_time


def search_copernicus_public(
    aoi_polygon: BaseGeometry, start_date: str, end_date: str
) -> list[dict[str, Any]]:
    """
    Searches CDSE Public OData.
    Fetches 'Footprint' to calculate true coverage area.
    """
    # OData requires strict WKT: "POLYGON((...))" (no space)
    wkt_aoi = aoi_polygon.wkt.replace("POLYGON ((", "POLYGON((")

    base_url = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"

    # Filter: Sentinel-2 L1C, Intersects AOI, Date Range, Cloud Cover < 20.0
    filter_query = (
        f"Collection/Name eq 'SENTINEL-2' and "
        f"contains(Name,'MSIL1C') and "
        f"OData.CSC.Intersects(area=geography'SRID=4326;{wkt_aoi}') and "
        f"ContentDate/Start ge {start_date} and "
        f"ContentDate/Start le {end_date} and "
        f"Attributes/OData.CSC.DoubleAttribute/any(att:att/Name eq 'cloudCover' and att/Value le {OPTICAL_CLOUD_THRESHOLD})"
    )

    # Removing $select ensures we get the 'Footprint' field
    params = {"$filter": filter_query, "$orderby": "ContentDate/Start asc", "$top": 1000}

    print("Searching Copernicus (Public OData) for Sentinel-2 L1C...")

    granules = []
    next_link = base_url
    session = requests.Session()

    while next_link:
        try:
            if next_link == base_url:
                response = session.get(next_link, params=params)
            else:
                response = session.get(next_link)

            if response.status_code != 200:
                print(f"Error URL: {response.url}")
                print(f"Response: {response.text}")
                response.raise_for_status()

            data = response.json()
            products = data.get("value", [])

            for prod in products:
                name = prod.get("Name")
                start = prod.get("ContentDate", {}).get("Start")

                # --- PARSE GEOMETRY ---
                # Format is: geography'SRID=4326;POLYGON ((...))'
                raw_footprint = prod.get("Footprint")
                footprint = None

                if raw_footprint:
                    try:
                        # Extract just the WKT part (POLYGON...)
                        # Split by semicolon to remove SRID
                        # Remove trailing quote if present
                        clean_wkt = raw_footprint.split(";")[-1].replace("'", "")
                        footprint = wkt.loads(clean_wkt)
                    except Exception:
                        pass

                # --- EXTRACT CLOUD COVER ---
                cloud_cover = 0.0
                attrs = prod.get("Attributes", [])
                for attr in attrs:
                    if attr.get("Name") == "cloudCover":
                        cloud_cover = attr.get("Value", 0.0)
                        break

                # --- DETERMINE PLATFORM ---
                platform = "sentinel-2"
                if name.startswith("S2A"):
                    platform = "sentinel-2a"
                elif name.startswith("S2B"):
                    platform = "sentinel-2b"
                elif name.startswith("S2C"):
                    platform = "sentinel-2c"

                granules.append(
                    {
                        "granule_id": name,
                        "date": start.split("T")[0],
                        "datetime": start,
                        "cloud_cover": cloud_cover,
                        "platform": platform,
                        "footprint": footprint,  # Successfully parsed geometry
                    }
                )

            next_link = data.get("@odata.nextLink")

        except Exception as e:
            print(f"Error searching Public OData: {e}")
            break

    print(f"  Found {len(granules)} Sentinel-2 products.")
    return granules


def find_optical_pairs_copernicus(
    optical_scenes: list[dict[str, Any]],
    rupture_time: int,
    title: str,
    event_id: str,
    aoi_polygon: BaseGeometry,
    job_list: bool = True,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """
    Generate Optical Pairs grouped by Relative Orbit.
    Prioritizes: 1. True Coverage Area %, 2. Cloud Cover, 3. Time.
    Uses the Copernicus data nomenclature

    Updates:
    - Uses make_optical_job_json for ARIA_AUTORIFT formatting.
    - Merges S2A/S2B/S2C (no platform filtering).
    - Creates 'partial' job entries if only one side (Pre/Post) is found.
    """
    rupture_dt = convert_time(rupture_time)
    aoi_area = aoi_polygon.area

    # Group by Orbit -> Date
    orbit_groups = defaultdict(lambda: defaultdict(list))

    for scene in optical_scenes:
        granule_id = scene["granule_id"]
        match = re.search(r"_R(\d{3})_", granule_id)
        if match:
            orbit_id = match.group(1)
            orbit_groups[orbit_id][scene["date"]].append(scene)

    jobs = []
    scene_features = []

    print(f"Found {len(orbit_groups)} unique satellite tracks.")

    for orbit_id, dates_dict in orbit_groups.items():
        pre_candidates = []
        post_candidates = []

        for date_str, scenes_on_date in dates_dict.items():
            date_dt = datetime.strptime(date_str, "%Y-%m-%d").replace(tzinfo=timezone.utc)

            # Combine scenes from same date (Mosaic logic)
            scenes = scenes_on_date
            platform = scenes[0]["platform"] if scenes else "sentinel-2"

            unique_scenes = []
            seen_tiles = set()
            scenes.sort(key=lambda x: x["cloud_cover"])

            combined_poly = None

            for s in scenes:
                t_match = re.search(r"_T(\w{5})_", s["granule_id"])
                if t_match:
                    tile_id = t_match.group(1)
                    if tile_id not in seen_tiles:
                        unique_scenes.append(s)
                        seen_tiles.add(tile_id)

                        if s.get("footprint"):
                            if combined_poly is None:
                                combined_poly = s["footprint"]
                            else:
                                combined_poly = combined_poly.union(s["footprint"])

            # Calculate Coverage
            coverage_pct = 0.0
            if combined_poly and aoi_polygon:
                try:
                    clean_poly = combined_poly.buffer(0)
                    intersection = clean_poly.intersection(aoi_polygon)
                    coverage_pct = (intersection.area / aoi_area) * 100.0
                except Exception:
                    pass

            tile_count = len(unique_scenes)
            avg_cc = (
                sum(s["cloud_cover"] for s in unique_scenes) / len(unique_scenes)
                if unique_scenes
                else 100
            )
            delta_days = abs((date_dt - rupture_dt).days)

            candidate = {
                "date": date_str,
                "scenes": unique_scenes,
                "platform": platform,
                "coverage": coverage_pct,
                "tile_count": tile_count,
                "cc": avg_cc,
                "delta": delta_days,
            }

            if date_dt < rupture_dt:
                pre_candidates.append(candidate)
            elif date_dt > rupture_dt:
                post_candidates.append(candidate)

        # --- HANDLE INCOMPLETE PAIRS ---
        if not pre_candidates or not post_candidates:
            print(f"  Orbit {orbit_id}: [PARTIAL] Incomplete pair.")

            pre_date_str = "MISSING"
            post_date_str = "MISSING"
            primary_ids = []
            secondary_ids = []

            if pre_candidates:
                best_pre = min(pre_candidates, key=lambda x: x["cc"])
                print(f"    - Found Pre-event: {best_pre['date']} (CC: {best_pre['cc']:.1f}%)")
                pre_date_str = best_pre["date"]
                primary_ids = [s["granule_id"].replace(".SAFE", "") for s in best_pre["scenes"]]
            else:
                print("    - Missing PRE-event coverage.")

            if post_candidates:
                best_post = min(post_candidates, key=lambda x: x["cc"])
                print(f"    - Found Post-event: {best_post['date']} (CC: {best_post['cc']:.1f}%)")
                post_date_str = best_post["date"]
                secondary_ids = [s["granule_id"].replace(".SAFE", "") for s in best_post["scenes"]]
            else:
                print("    - Missing POST-event coverage.")

            # Generate the partial job and add to list
            if job_list:
                job = make_optical_job_json(
                    title,
                    event_id,
                    orbit_id,
                    pre_date_str,
                    post_date_str,
                    primary_ids,
                    secondary_ids,
                    status="PARTIAL",
                )
                jobs.append(job)
            continue

        # --- STANDARD SELECTION ---
        # Sort: Max Coverage > Low Cloud > Low Delta
        pre_candidates.sort(key=lambda x: (-x["coverage"], x["cc"], x["delta"]))
        best_pre = pre_candidates[0]

        post_candidates.sort(key=lambda x: (-x["coverage"], x["cc"], x["delta"]))
        best_post = post_candidates[0]

        print(
            f"  Orbit {orbit_id}: Pre={best_pre['date']} ({best_pre['coverage']:.1f}% Cov), Post={best_post['date']} ({best_post['coverage']:.1f}% Cov)"
        )

        for stage, candidate in [("Pre-Event", best_pre), ("Post-Event", best_post)]:
            for scene in candidate["scenes"]:
                if scene.get("footprint"):
                    feature = {
                        "type": "Feature",
                        "geometry": mapping(scene["footprint"]),
                        "properties": {
                            "earthquake": title,
                            "role": stage,
                            "date": candidate["date"],
                            "orbit": orbit_id,
                            "platform": candidate["platform"],
                            "granule_id": scene["granule_id"],
                            "cloud_cover": scene["cloud_cover"],
                            "coverage_pct": round(candidate["coverage"], 1),
                        },
                    }
                    scene_features.append(feature)

        if job_list:
            primary_ids = [s["granule_id"].replace(".SAFE", "") for s in best_pre["scenes"]]
            secondary_ids = [s["granule_id"].replace(".SAFE", "") for s in best_post["scenes"]]

            job = make_optical_job_json(
                title,
                event_id,
                orbit_id,
                best_pre["date"],
                best_post["date"],
                primary_ids,
                secondary_ids,
            )
            jobs.append(job)

    return jobs, scene_features
