"""Per-earthquake processing: AOI, SAR pairing or optical backends."""

from __future__ import annotations

import json
import math
import os
from datetime import datetime, timedelta
from typing import Any

import geojson
import geopandas as gpd
from shapely.geometry import mapping

from aria_coseis import config
from aria_coseis.aoi import load_aoi_from_json, make_aoi
from aria_coseis.notify import make_interactive_map
from aria_coseis.optical.copernicus import find_optical_pairs_copernicus, search_copernicus_public
from aria_coseis.optical.element84 import find_optical_pairs_element84, search_element84_stac
from aria_coseis.optical.gee import (
    assign_nodata,
    download_from_gcs,
    export_gee_landsat_composite,
    export_gee_sentinel2_composite,
    merge_and_compress_chips,
    wait_for_gee_tasks,
)
from aria_coseis.sar.pairing import find_reference_and_secondary_pairs
from aria_coseis.sar.search import get_path_and_frame_numbers, get_SLCs
from aria_coseis.usgs import get_ffm_geojson_url
from aria_coseis.utils import convert_time, to_snake_case


def process_earthquake(
    eq: dict[str, Any],
    aoi: str | None,
    pairing_mode: str | None,
    job_list: bool,
    resolution: int = 90,
    sensor: str = "sar",
    optical_backend: str = "copernicus",
    optical_level: str = "toa",
) -> tuple[list[list[dict[str, Any]]], list[dict[str, Any]]]:
    """
    Process earthquake event and generate the necessary SLC pairs for InSAR processing.
    :param eq: dictionary containing earthquake data
    :param aoi: Area of Interest (AOI) as a GeoJSON file
    :param pairing_mode: 'all', 'sequential', or 'coseismic' for specifying desired SLC pairing
    :param job_list: True if the JSON objects are for HYP3 job submission, False otherwise
    :param resolution: Output resolution for the topsApp processing, default is 90m
    :param sensor: 'sar' for SAR processing, 'sentinel-2' or 'landsat' for optical processing
    :optical_backend: 'copernicus', 'element84' for source data file nomenclature (only applicable if sensor is 'optical')
    :optical_level: 'toa', 'sr', or 'raw' for the desired optical product level (only applicable if sensor is 'optical')
    :return: List of JSON objects containing the parameters for each pair of SLCs
    """
    title = eq.get("title", "")
    title = to_snake_case(title)
    print(f"title: {title}")
    coords = eq.get("coordinates", [])
    event_id = eq.get("id", "")

    # Get the FFM geometry (if it exists)
    ffm_url = get_ffm_geojson_url(event_id)

    if ffm_url:
        print("FFM URL found. Loading AOI from FFM GeoJSON...")
        # Load the FFM geometry
        aoi = load_aoi_from_json(ffm_url)

        # Buffer the FFM to capture the full deformation field (optical only; SAR frame
        # selection on main is tuned to the unbuffered FFM). .envelope forces a clean rectangle.
        if sensor != "sar":
            buffer_deg = 0.15
            aoi = aoi.buffer(buffer_deg).envelope
            print(f"  -> Buffered FFM bounds by {buffer_deg} degrees to ensure coverage.")

    elif aoi:
        print("AOI provided. Using the provided AOI...")
        # Load the AOI from the provided JSON file
        aoi = load_aoi_from_json(aoi)

    # Generate AOI or use the user-provided AOI
    else:
        aoi = make_aoi(coords)  # Create AOI if not provided

    # Write AOI to a geojson file
    with open(config.output_path(f"{title}_{sensor}_{optical_level}_AOI.geojson"), "w") as f:
        geojson.dump(mapping(aoi), f, indent=2)

    all_jobs = []
    all_features = []

    if sensor == "sar":
        path_frame_numbers, frame_dataframe = get_path_and_frame_numbers(aoi, eq.get("time"))

        # Frame visualization logic (SAR specific)
        if not job_list:
            frame_gdf = gpd.GeoDataFrame(frame_dataframe, geometry="geometry", crs="EPSG:4326")
            for col in frame_gdf.columns:
                if frame_gdf[col].apply(lambda x: isinstance(x, list)).any():
                    frame_gdf[col] = frame_gdf[col].astype(str)
            frame_gdf.to_file(config.output_path(f"{title}_frames.geojson"), driver="GeoJSON")
            make_interactive_map(
                frame_dataframe, eq.get("title", ""), eq.get("coordinates", []), eq.get("url", "")
            )

        for (flight_direction, path_number), frame_numbers in path_frame_numbers.items():
            frame_numbers = list(set(fn[0] for fn in frame_numbers))
            SLCs = get_SLCs(
                flight_direction, path_number, aoi.wkt, eq.get("time"), processing_mode="historic"
            )
            isce_jobs = find_reference_and_secondary_pairs(
                SLCs,
                eq.get("time"),
                flight_direction,
                path_number,
                title,
                aoi,
                event_id,
                pairing_mode,
                job_list,
                resolution,
            )
            all_jobs.append(isce_jobs)

    elif sensor in ["sentinel-2", "landsat"]:
        rupture_time = eq.get("time")
        rupture_dt = convert_time(rupture_time)

        start_search = (rupture_dt - timedelta(days=90)).strftime("%Y-%m-%dT%H:%M:%SZ")
        end_search = (rupture_dt + timedelta(days=90)).strftime("%Y-%m-%dT%H:%M:%SZ")

        if optical_backend == "copernicus":
            print("Routing to Copernicus Public OData backend...")
            s2_scenes = search_copernicus_public(aoi, start_search, end_search)
            s2_jobs, s2_features = find_optical_pairs_copernicus(
                s2_scenes, rupture_time, title, event_id, aoi, job_list
            )

        elif optical_backend == "element84":
            print("Routing to Element84 STAC backend...")
            s2_scenes = search_element84_stac(aoi, start_search, end_search)
            s2_jobs, f_dom, f_split = find_optical_pairs_element84(
                s2_scenes, rupture_time, title, event_id, aoi
            )
            s2_features = f_dom + f_split

        elif optical_backend == "gee":
            print("Routing to Google Earth Engine backend...")

            # Define 90-day pre- and post-event temporal windows
            rupture_dt = convert_time(rupture_time).replace(tzinfo=None)
            pre_start = (rupture_dt - timedelta(days=90)).strftime("%Y-%m-%d")

            # Keep exact UTC time for the rupture boundaries
            pre_end = rupture_dt.strftime("%Y-%m-%dT%H:%M:%S")
            post_start = rupture_dt.strftime("%Y-%m-%dT%H:%M:%S")
            post_end = (rupture_dt + timedelta(days=90)).strftime("%Y-%m-%d")

            lon, lat = coords[0], coords[1]
            utm_zone = math.floor((lon + 180) / 6) + 1
            epsg_base = 32600 if lat >= 0 else 32700
            target_crs = f"EPSG:{epsg_base + utm_zone}"

            # --- MISSION DECISION LOGIC ---
            if sensor == "sentinel-2":
                if optical_level.lower() == "sr":
                    s2_collection = "COPERNICUS/S2_SR_HARMONIZED"
                elif optical_level.lower() == "toa":
                    s2_collection = "COPERNICUS/S2_HARMONIZED"
                else:
                    print("  [Warning] Sentinel-2 does not support 'raw'. Defaulting to 'toa'.")
                    s2_collection = "COPERNICUS/S2_HARMONIZED"
                    optical_level = "toa"
                print(
                    f"  [Sentinel-2 Setup] Level: {optical_level.upper()} | Collection: {s2_collection} | Band: B8"
                )

            elif sensor == "landsat":
                pre_dt = datetime.strptime(pre_start, "%Y-%m-%d")
                post_dt = datetime.strptime(post_end, "%Y-%m-%d")

                mission = "L5"
                if post_dt < datetime(2012, 1, 1):
                    mission = "L5"
                elif pre_dt >= datetime(2013, 4, 15):
                    mission = "L8"
                else:
                    mission = "L7"

                if optical_level.lower() == "raw":
                    if mission == "L5":
                        landsat_collection, landsat_band, landsat_scale = (
                            "LANDSAT/LT05/C02/T1",
                            "B2",
                            30,
                        )
                    if mission == "L7":
                        landsat_collection, landsat_band, landsat_scale = (
                            "LANDSAT/LE07/C02/T1",
                            "B8",
                            15,
                        )
                    if mission == "L8":
                        landsat_collection, landsat_band, landsat_scale = (
                            "LANDSAT/LC08/C02/T1",
                            "B8",
                            15,
                        )
                elif optical_level.lower() == "sr":
                    # SR doesn't have Pan. Switching to Red Band (30m). AutoRIFT will dynamically adjust to 30m.
                    if mission == "L5":
                        landsat_collection, landsat_band, landsat_scale = (
                            "LANDSAT/LT05/C02/T1_L2",
                            "SR_B3",
                            30,
                        )
                    if mission == "L7":
                        landsat_collection, landsat_band, landsat_scale = (
                            "LANDSAT/LE07/C02/T1_L2",
                            "SR_B3",
                            30,
                        )
                    if mission == "L8":
                        landsat_collection, landsat_band, landsat_scale = (
                            "LANDSAT/LC08/C02/T1_L2",
                            "SR_B4",
                            30,
                        )
                else:
                    optical_level = "toa"
                    if mission == "L5":
                        landsat_collection, landsat_band, landsat_scale = (
                            "LANDSAT/LT05/C02/T1_TOA",
                            "B2",
                            30,
                        )
                    if mission == "L7":
                        landsat_collection, landsat_band, landsat_scale = (
                            "LANDSAT/LE07/C02/T1_TOA",
                            "B8",
                            15,
                        )
                    if mission == "L8":
                        landsat_collection, landsat_band, landsat_scale = (
                            "LANDSAT/LC08/C02/T1_TOA",
                            "B8",
                            15,
                        )

                print(
                    f"  [Landsat Setup] Level: {optical_level.upper()} | Collection: {landsat_collection} | Band: {landsat_band}"
                )

            # Create job_list, if applicable
            if job_list:
                print(f"  Generating GEE Job payload for {title} (Skipping computation).")
                gee_job = {
                    "name": f"{title}-GEE_{sensor.upper()}",
                    "job_type": "GEE_OPTICAL_COSEIS",
                    "event_id": event_id,
                    "job_parameters": {
                        "event_title": title,
                        "sensor": sensor,
                        "target_crs": target_crs,
                        "pre_start": pre_start,
                        "pre_end": pre_end,
                        "post_start": post_start,
                        "post_end": post_end,
                        # Inject Landsat parameters if applicable
                        "landsat_collection": landsat_collection if sensor == "landsat" else None,
                        "landsat_band": landsat_band if sensor == "landsat" else None,
                        "landsat_scale": landsat_scale if sensor == "landsat" else None,
                    },
                }
                return [[gee_job]], [
                    {
                        "type": "Feature",
                        "geometry": mapping(aoi),
                        "properties": {"title": title, "crs": target_crs},
                    }
                ]

            # Execution logic
            local_dir = os.path.join(config.root_dir, "GEE_Optical_Downloads", title)
            manifest_path = os.path.join(
                local_dir, f"{title}_{sensor}_{optical_level.lower()}_autorift_manifest.json"
            )
            if os.path.exists(manifest_path):
                print(f"  Data already downloaded for {title}. Skipping GEE computation.")
                return [], []

            gcs_bucket = os.getenv("COSEIS_GCS_BUCKET")
            if not gcs_bucket:
                print("Error: COSEIS_GCS_BUCKET environment variable is not set.")
                return [], []

            import ee

            try:
                print("initializing with coseis-1")
                ee.Initialize(project="coseis-1")
            except Exception as e:
                print(
                    "Earth Engine not authenticated. Run 'earthengine authenticate --auth_mode=notebook' in your terminal."
                )
                raise e

            if sensor == "sentinel-2":
                print("Generating Pre-Event Sentinel-2 Composites...")
                pre_exports, pre_dates = export_gee_sentinel2_composite(
                    aoi,
                    pre_start,
                    pre_end,
                    title,
                    "PRE",
                    gcs_bucket,
                    s2_collection,
                    optical_level,
                    crs_epsg=target_crs,
                )
                print("Generating Post-Event Sentinel-2 Composites...")
                post_exports, post_dates = export_gee_sentinel2_composite(
                    aoi,
                    post_start,
                    post_end,
                    title,
                    "POST",
                    gcs_bucket,
                    s2_collection,
                    optical_level,
                    crs_epsg=target_crs,
                )

            elif sensor == "landsat":
                print("Generating Pre-Event Landsat Composites...")
                pre_exports, pre_dates = export_gee_landsat_composite(
                    aoi,
                    pre_start,
                    pre_end,
                    title,
                    "PRE",
                    gcs_bucket,
                    landsat_collection,
                    landsat_band,
                    landsat_scale,
                    optical_level,
                    crs_epsg=target_crs,
                )
                print("Generating Post-Event Landsat Composites...")
                post_exports, post_dates = export_gee_landsat_composite(
                    aoi,
                    post_start,
                    post_end,
                    title,
                    "POST",
                    gcs_bucket,
                    landsat_collection,
                    landsat_band,
                    landsat_scale,
                    optical_level,
                    crs_epsg=target_crs,
                )

            # Extract task objects from the dictionaries to monitor them
            all_tasks = [v["task"] for v in pre_exports.values()] + [
                v["task"] for v in post_exports.values()
            ]
            wait_for_gee_tasks(all_tasks)

            # Download the files using a Prefix Search
            print("\nDownloading independent track composites from Google Cloud Storage...")
            track_pairs = {}

            # Find the overlapping tracks that have both Pre and Post data
            valid_tracks = set(pre_exports.keys()).intersection(set(post_exports.keys()))

            for track in valid_tracks:
                print(f"Processing Track {track}...")

                # Download Pre-event for this specific track
                pre_local_paths = download_from_gcs(
                    gcs_bucket, pre_exports[track]["prefix"], local_dir
                )
                full_pre_filename = os.path.basename(pre_exports[track]["prefix"])
                final_pre_path = os.path.join(local_dir, f"{full_pre_filename}.tif")

                if len(pre_local_paths) > 1:
                    merge_and_compress_chips(pre_local_paths, final_pre_path)
                elif len(pre_local_paths) == 1:
                    os.rename(pre_local_paths[0], final_pre_path)

                # Download Post-event for this specific track
                post_local_paths = download_from_gcs(
                    gcs_bucket, post_exports[track]["prefix"], local_dir
                )
                full_post_filename = os.path.basename(post_exports[track]["prefix"])
                final_post_path = os.path.join(local_dir, f"{full_post_filename}.tif")

                if len(post_local_paths) > 1:
                    merge_and_compress_chips(post_local_paths, final_post_path)
                elif len(post_local_paths) == 1:
                    os.rename(post_local_paths[0], final_post_path)

                # --- NEW VALIDATION LOGIC ---
                # Verify files exist before assigning nodata or adding to manifest
                if not os.path.exists(final_pre_path) or not os.path.exists(final_post_path):
                    print(f"  Warning: Missing data for track {track}. Skipping pair in manifest.")
                    continue
                # ----------------------------

                assign_nodata(final_pre_path, nodata_val=0)
                assign_nodata(final_post_path, nodata_val=0)

                track_pairs[track] = {"pre_image": final_pre_path, "post_image": final_post_path}

            # Ensure the directory exists (in case download_from_gcs was never triggered)
            os.makedirs(local_dir, exist_ok=True)

            # --- HANDLE NO DATA SCENARIO ---
            if not track_pairs:
                print(
                    f"\n  No valid tracks downloaded for {title}. Writing failed manifest to prevent retries."
                )
                manifest_payload = {
                    "event_title": title,
                    "event_id": event_id,
                    "sensor": sensor,
                    "optical_level": optical_level,
                    "backend": "Google Earth Engine",
                    "track_pairs": {},
                    "status": "FAILED_NO_DATA",
                }
                with open(manifest_path, "w") as f:
                    json.dump(manifest_payload, f, indent=4)
                return [], []
            # -------------------------------

            # Create the local manifest for AutoRIFT (Success Case)
            manifest_payload = {
                "event_title": title,
                "event_id": event_id,
                "sensor": sensor,
                "optical_level": optical_level,
                "backend": "Google Earth Engine",
                "track_pairs": track_pairs,
                "pre_dates_used": pre_dates,
                "post_dates_used": post_dates,
                "status": "DOWNLOADED_READY_FOR_AUTORIFT",
            }

            with open(manifest_path, "w") as f:
                json.dump(manifest_payload, f, indent=4)

            print(f"\nManifest written to: {manifest_path}")

            return [], []

        if s2_jobs:
            print(f"Generated {len(s2_jobs)} Sentinel-2 jobs using {optical_backend}.")
            all_jobs.append(s2_jobs)
            all_features.extend(s2_features)

            if s2_features:
                fc = {"type": "FeatureCollection", "features": s2_features}
                scene_filename = config.output_path(f"{title}_selected_scenes.geojson")
                with open(scene_filename, "w") as f:
                    json.dump(fc, f, indent=2)
                print(f"Saved selected scene footprints to {scene_filename}")
        else:
            print(f"No optical pairs found using {optical_backend}.")

    return all_jobs, all_features
