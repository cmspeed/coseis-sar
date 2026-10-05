import re
import unicodedata
import os
import argparse
import asf_search as asf
from dateutil import parser as dateparser
from dateutil.parser import isoparse
from pathlib import Path
import requests
import json
import folium
import geojson
import geopandas as gpd
import glob
import csv
import math
from shapely import wkt
from shapely.geometry import mapping, shape, box, Point, Polygon, LineString, MultiLineString, MultiPolygon
from shapely.ops import unary_union, linemerge
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from collections import defaultdict
from itertools import combinations
import logging
import yagmail
import time
from time import sleep
from types import SimpleNamespace
from urllib.parse import urlparse
from typing import List, Dict, Any, Optional

from aria_coseis import config
from aria_coseis.config import (
    ASF_DAAC_API,
    OPTICAL_CLOUD_THRESHOLD,
    USGS_api_alltime,
    coastline_api,
)
from aria_coseis.utils import convert_time, to_snake_case
from aria_coseis.notify import ascii_table_to_html, get_next_pass, make_interactive_map, send_email
from aria_coseis.usgs import (
    get_event_rake,
    get_ffm_geojson_url,
    get_historic_earthquake_data_date_range,
    get_historic_earthquake_data_single_date,
    parse_custom_eq_list,
    parse_geojson,
)
from aria_coseis.aoi import load_aoi_from_json, make_aoi
from aria_coseis.significance import check_significance
from aria_coseis.sar.search import get_SLCs, get_path_and_frame_numbers
from aria_coseis.sar.pairing import find_reference_and_secondary_pairs, make_job_json
from aria_coseis.sar.topsapp import run_dockerized_topsApp
from aria_coseis.optical.jobs import get_utm_zone, make_optical_job_json
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
from aria_coseis.tracking import add_to_tracker, check_tracker_for_updates, load_tracker


def process_earthquake(eq, aoi, pairing_mode, job_list, resolution=90, sensor='sar', optical_backend='copernicus', optical_level='toa'):
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
    title = eq.get('title', '')
    title = to_snake_case(title)
    print(f"title: {title}")
    coords = eq.get('coordinates', [])
    event_id = eq.get('id', '')

    # Get the FFM geometry (if it exists)
    ffm_url = get_ffm_geojson_url(event_id)

    if ffm_url:
        print("FFM URL found. Loading AOI from FFM GeoJSON...")
        # Load the FFM geometry
        aoi = load_aoi_from_json(ffm_url)
        
        # Buffer the FFM to capture the full deformation field (optical only; SAR frame
        # selection on main is tuned to the unbuffered FFM). .envelope forces a clean rectangle.
        if sensor != 'sar':
            buffer_deg = 0.15
            aoi = aoi.buffer(buffer_deg).envelope
            print(f"  -> Buffered FFM bounds by {buffer_deg} degrees to ensure coverage.")

    elif aoi:
        print("AOI provided. Using the provided AOI...")
        # Load the AOI from the provided JSON file
        aoi = load_aoi_from_json(aoi)

    # Generate AOI or use the user-provided AOI
    else:
        aoi = make_aoi(coords) # Create AOI if not provided

    # Write AOI to a geojson file
    with open(f'{title}_{sensor}_{optical_level}_AOI.geojson', 'w') as f:
        geojson.dump(mapping(aoi), f, indent=2)

    all_jobs = []
    all_features = []

    if sensor == 'sar':
        path_frame_numbers, frame_dataframe = get_path_and_frame_numbers(aoi, eq.get('time'))
        
        # Frame visualization logic (SAR specific)
        if not job_list:
            frame_gdf = gpd.GeoDataFrame(frame_dataframe, geometry="geometry", crs="EPSG:4326")
            for col in frame_gdf.columns:
                if frame_gdf[col].apply(lambda x: isinstance(x, list)).any():
                    frame_gdf[col] = frame_gdf[col].astype(str)
            frame_gdf.to_file(f"{title}_frames.geojson", driver="GeoJSON")
            make_interactive_map(frame_dataframe, eq.get('title', ''), eq.get('coordinates', []), eq.get('url', ''))

        for (flight_direction, path_number), frame_numbers in path_frame_numbers.items():
            frame_numbers = list(set(fn[0] for fn in frame_numbers))
            SLCs = get_SLCs(flight_direction, path_number, aoi.wkt, eq.get('time'), processing_mode='historic')
            isce_jobs = find_reference_and_secondary_pairs(SLCs, eq.get('time'), flight_direction, path_number, 
                                                           title, aoi, event_id, pairing_mode, job_list, resolution)
            all_jobs.append(isce_jobs)

    elif sensor in ['sentinel-2', 'landsat']:
        rupture_time = eq.get('time')
        rupture_dt = convert_time(rupture_time)
        
        start_search = (rupture_dt - timedelta(days=90)).strftime('%Y-%m-%dT%H:%M:%SZ')
        end_search = (rupture_dt + timedelta(days=90)).strftime('%Y-%m-%dT%H:%M:%SZ')

        if optical_backend == 'copernicus':
            print("Routing to Copernicus Public OData backend...")
            s2_scenes = search_copernicus_public(aoi, start_search, end_search)
            s2_jobs, s2_features = find_optical_pairs_copernicus(s2_scenes, rupture_time, title, event_id, aoi, job_list)
            
        elif optical_backend == 'element84':
            print("Routing to Element84 STAC backend...")
            s2_scenes = search_element84_stac(aoi, start_search, end_search)
            s2_jobs, f_dom, f_split = find_optical_pairs_element84(s2_scenes, rupture_time, title, event_id, aoi)
            s2_features = f_dom + f_split
            
        elif optical_backend == 'gee':
            print("Routing to Google Earth Engine backend...")
            
            # Define 60-day temporal windows
            rupture_dt = convert_time(rupture_time).replace(tzinfo=None)
            pre_start = (rupture_dt - timedelta(days=90)).strftime('%Y-%m-%d')
            
            # Keep exact UTC time for the rupture boundaries
            pre_end = rupture_dt.strftime('%Y-%m-%dT%H:%M:%S')
            post_start = rupture_dt.strftime('%Y-%m-%dT%H:%M:%S')
            post_end = (rupture_dt + timedelta(days=90)).strftime('%Y-%m-%d')

            lon, lat = coords[0], coords[1]
            utm_zone = math.floor((lon + 180) / 6) + 1
            epsg_base = 32600 if lat >= 0 else 32700
            target_crs = f"EPSG:{epsg_base + utm_zone}"
            
            # --- MISSION DECISION LOGIC ---
            if sensor == 'sentinel-2':
                if optical_level.lower() == 'sr':
                    s2_collection = 'COPERNICUS/S2_SR_HARMONIZED'
                elif optical_level.lower() == 'toa':
                    s2_collection = 'COPERNICUS/S2_HARMONIZED'
                else:
                    print("  [Warning] Sentinel-2 does not support 'raw'. Defaulting to 'toa'.")
                    s2_collection = 'COPERNICUS/S2_HARMONIZED'
                    optical_level = 'toa'
                print(f"  [Sentinel-2 Setup] Level: {optical_level.upper()} | Collection: {s2_collection} | Band: B8")
            
            elif sensor == 'landsat':
                pre_dt = datetime.strptime(pre_start, '%Y-%m-%d')
                post_dt = datetime.strptime(post_end, '%Y-%m-%d')
                
                mission = 'L5'
                if post_dt < datetime(2012, 1, 1):
                    mission = 'L5'
                elif pre_dt >= datetime(2013, 4, 15):
                    mission = 'L8'
                else:
                    mission = 'L7'

                if optical_level.lower() == 'raw':
                    if mission == 'L5': landsat_collection, landsat_band, landsat_scale = 'LANDSAT/LT05/C02/T1', 'B2', 30
                    if mission == 'L7': landsat_collection, landsat_band, landsat_scale = 'LANDSAT/LE07/C02/T1', 'B8', 15
                    if mission == 'L8': landsat_collection, landsat_band, landsat_scale = 'LANDSAT/LC08/C02/T1', 'B8', 15
                elif optical_level.lower() == 'sr':
                    # SR doesn't have Pan. Switching to Red Band (30m). AutoRIFT will dynamically adjust to 30m.
                    if mission == 'L5': landsat_collection, landsat_band, landsat_scale = 'LANDSAT/LT05/C02/T1_L2', 'SR_B3', 30
                    if mission == 'L7': landsat_collection, landsat_band, landsat_scale = 'LANDSAT/LE07/C02/T1_L2', 'SR_B3', 30
                    if mission == 'L8': landsat_collection, landsat_band, landsat_scale = 'LANDSAT/LC08/C02/T1_L2', 'SR_B4', 30
                else:
                    optical_level = 'toa'
                    if mission == 'L5': landsat_collection, landsat_band, landsat_scale = 'LANDSAT/LT05/C02/T1_TOA', 'B2', 30
                    if mission == 'L7': landsat_collection, landsat_band, landsat_scale = 'LANDSAT/LE07/C02/T1_TOA', 'B8', 15
                    if mission == 'L8': landsat_collection, landsat_band, landsat_scale = 'LANDSAT/LC08/C02/T1_TOA', 'B8', 15
                
                print(f"  [Landsat Setup] Level: {optical_level.upper()} | Collection: {landsat_collection} | Band: {landsat_band}")

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
                        "landsat_collection": landsat_collection if sensor == 'landsat' else None,
                        "landsat_band": landsat_band if sensor == 'landsat' else None,
                        "landsat_scale": landsat_scale if sensor == 'landsat' else None
                    }
                }
                return [[gee_job]], [{"type": "Feature", "geometry": mapping(aoi), "properties": {"title": title, "crs": target_crs}}]

            # Execution logic
            local_dir = os.path.join(config.root_dir, "GEE_Optical_Downloads", title)
            manifest_path = os.path.join(local_dir, f"{title}_{sensor}_{optical_level.lower()}_autorift_manifest.json")
            if os.path.exists(manifest_path):
                print(f"  Data already downloaded for {title}. Skipping GEE computation.")
                return [], []
            
            gcs_bucket = os.getenv('COSEIS_GCS_BUCKET')
            if not gcs_bucket:
                print("Error: COSEIS_GCS_BUCKET environment variable is not set.")
                return [], []

            import ee
            try:
                print('initializing with coseis-1')
                ee.Initialize(project='coseis-1')
            except Exception as e:
                print("Earth Engine not authenticated. Run 'earthengine authenticate --auth_mode=notebook' in your terminal.")
                raise e

            if sensor == 'sentinel-2':
                print(f"Generating Pre-Event Sentinel-2 Composites...")
                pre_exports, pre_dates = export_gee_sentinel2_composite(
                    aoi, pre_start, pre_end, title, "PRE", gcs_bucket, s2_collection, optical_level, crs_epsg=target_crs
                )
                print(f"Generating Post-Event Sentinel-2 Composites...")
                post_exports, post_dates = export_gee_sentinel2_composite(
                    aoi, post_start, post_end, title, "POST", gcs_bucket, s2_collection, optical_level, crs_epsg=target_crs
                )
            
            elif sensor == 'landsat':
                print(f"Generating Pre-Event Landsat Composites...")
                pre_exports, pre_dates = export_gee_landsat_composite(
                    aoi, pre_start, pre_end, title, "PRE", gcs_bucket, landsat_collection,
                    landsat_band, landsat_scale, optical_level, crs_epsg=target_crs
                )
                print(f"Generating Post-Event Landsat Composites...")
                post_exports, post_dates = export_gee_landsat_composite(
                    aoi, post_start, post_end, title, "POST", gcs_bucket, landsat_collection,
                    landsat_band, landsat_scale, optical_level, crs_epsg=target_crs
                )

            # Extract task objects from the dictionaries to monitor them
            all_tasks = [v['task'] for v in pre_exports.values()] + [v['task'] for v in post_exports.values()]
            wait_for_gee_tasks(all_tasks)

            # Download the files using a Prefix Search
            print("\nDownloading independent track composites from Google Cloud Storage...")
            track_pairs = {}
            
            # Find the overlapping tracks that have both Pre and Post data
            valid_tracks = set(pre_exports.keys()).intersection(set(post_exports.keys()))
            
            for track in valid_tracks:
                print(f"Processing Track {track}...")
                
                # Download Pre-event for this specific track
                pre_local_paths = download_from_gcs(gcs_bucket, pre_exports[track]['prefix'], local_dir)
                full_pre_filename = os.path.basename(pre_exports[track]['prefix'])
                final_pre_path = os.path.join(local_dir, f"{full_pre_filename}.tif")
                
                if len(pre_local_paths) > 1:
                    merge_and_compress_chips(pre_local_paths, final_pre_path)
                elif len(pre_local_paths) == 1:
                    os.rename(pre_local_paths[0], final_pre_path)
                    
                # Download Post-event for this specific track
                post_local_paths = download_from_gcs(gcs_bucket, post_exports[track]['prefix'], local_dir)
                full_post_filename = os.path.basename(post_exports[track]['prefix'])
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
                
                track_pairs[track] = {
                    "pre_image": final_pre_path,
                    "post_image": final_post_path
                }

            # Ensure the directory exists (in case download_from_gcs was never triggered)
            os.makedirs(local_dir, exist_ok=True)

            # --- HANDLE NO DATA SCENARIO ---
            if not track_pairs:
                print(f"\n  No valid tracks downloaded for {title}. Writing failed manifest to prevent retries.")
                manifest_payload = {
                    "event_title": title,
                    "event_id": event_id,
                    "sensor": sensor,
                    "optical_level": optical_level,
                    "backend": "Google Earth Engine",
                    "track_pairs": {},
                    "status": "FAILED_NO_DATA"
                }
                with open(manifest_path, 'w') as f:
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
                "status": "DOWNLOADED_READY_FOR_AUTORIFT"
            }
            
            with open(manifest_path, 'w') as f:
                json.dump(manifest_payload, f, indent=4)
                
            print(f"\nManifest written to: {manifest_path}")

            return [], []
        
        if s2_jobs:
            print(f"Generated {len(s2_jobs)} Sentinel-2 jobs using {optical_backend}.")
            all_jobs.append(s2_jobs)
            all_features.extend(s2_features)
        
            if s2_features:
                fc = {"type": "FeatureCollection", "features": s2_features}
                scene_filename = f"{title}_selected_scenes.geojson"
                with open(scene_filename, 'w') as f:
                    json.dump(fc, f, indent=2)
                print(f"Saved selected scene footprints to {scene_filename}")
        else:
            print(f"No optical pairs found using {optical_backend}.")

    return all_jobs, all_features


def main_forward(pairing_mode=None, resolution=30, do_processing=False, send_email_flag=False, process_only=False):
    """
    Runs the main query and processing workflow in forward processing mode.
    Used to produce co-seismic product for new earthquakes when new SLC data becomes available.
    :param pairing_mode: 'all', 'sequential', or 'coseismic' for specifying desired SLC pairing
    :param resolution: Output resolution for the topsApp processing, default is 30m
    :param do_processing: If True, runs the dockerized topsApp processing workflow after generating the JSONs. Default is False.
    :param send_email_flag: If True, sends an email alert after processing. Default is False.
    :param process_only: If True, only runs the processing workflow without generating new JSONs or sending emails. Default is False.
    """
    import shutil

    # A lock file to prevent overlapping runs
    lock_file = "/tmp/coseis_processing.lock"

    if os.path.exists(lock_file):
        print("Previous processing run still active. Exiting.")
        return
    try:
        with open(lock_file, 'w') as f:
            f.write("running")

        print('=========================================')
        print("Running cronjob to check for new earthquakes...")
        print('=========================================')
        
        # Initialize the tracking directory if it doesn't exist
        if not os.path.exists(config.TRACKING_DIR):
            os.makedirs(config.TRACKING_DIR, exist_ok=True)

        if not process_only:
            # Check for New Earthquakes over 48-hour window to ensre no events are missed due to API delays
            two_days_ago = (datetime.now(timezone.utc) - timedelta(days=2)).strftime('%Y-%m-%dT%H:%M:%S')
            
            # Define parameters for a custom search on the USGS 'alltime' endpoint
            params = {
                "format": "geojson",
                "starttime": two_days_ago,
                "minmagnitude": 5.5
            }
            print(f"Checking for earthquakes since {two_days_ago}...")

            # Use the query endpoint instead of the static summary feeds
            response = requests.get(USGS_api_alltime, params=params)
            response.raise_for_status()
            geojson_data = response.json()

            start_date = datetime.now().strftime('%Y-%m-%d')
            current_time = datetime.now(timezone.utc).strftime("%Y-%m-%d at %H:%M:%S UTC")

            if geojson_data:
                # Parse GeoJSON and create variables for each feature's properties
                earthquakes = parse_geojson(geojson_data)
                eq_sig = check_significance(earthquakes, start_date, end_date=None, mode = 'forward')

                if eq_sig is not None:
                    for eq in eq_sig:
                        # Check for duplicate entry in the pending queue
                        tracker = load_tracker()
                        if eq.get('id') in tracker:
                            print(f"Earthquake with ID {eq.get('id')} is already in the pending queue. Skipping.")
                            continue

                        title = eq.get('title', '')
                        title_snake = to_snake_case(title)
                        print(f"title: {title_snake}")
                        coords = eq.get('coordinates', [])
                        
                        # Initial AOI creation (a 1-degree box)
                        aoi = make_aoi(coords)

                        # Write AOI to a geojson file
                        with open(f'{title}_AOI.geojson', 'w') as f:
                            geojson.dump(aoi, f, indent=2)

                        # Get path/frame numbers for the initial AOI
                        path_frame_numbers, frame_dataframe = get_path_and_frame_numbers(aoi, eq.get('time'))
                        
                        # Create a timestamp string
                        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

                        # Create the output directory
                        timestamp_dir = Path(f"nextpass_outputs_{timestamp}")
                        timestamp_dir.mkdir(parents=True, exist_ok=True)

                        # Run next_pass to get the next overpasses
                        s1_info, nisar_info, overpass_map = get_next_pass(aoi, timestamp_dir)
                        
                        # Setup Github pages directory
                        docs_maps_dir = Path(os.getcwd()).parent / "docs" / "maps"
                        docs_maps_dir.mkdir(parents=True, exist_ok=True)
                        GITHUB_PAGES_BASE_URL = "https://cmspeed.github.io/coseis-sar"
                        
                        # Create a unique ID for the filenames so they aren't overwritten
                        unique_id = eq.get('id', datetime.now().strftime('%Y%m%d%H%M%S'))

                        # Route Joint Map
                        map_url = ""
                        if overpass_map and os.path.exists(overpass_map):
                            new_map_name = f"{title_snake}_{unique_id}_overpass_map.html"
                            shutil.copy(overpass_map, docs_maps_dir / new_map_name)
                            map_url = f"{GITHUB_PAGES_BASE_URL}/maps/{new_map_name}"

                        # Construct and send the initial email alert
                        message_dict = {
                            "title": eq.get('title', ''),
                            "time": convert_time(eq['time']).strftime('%Y-%m-%d %H:%M:%S'),
                            "coordinates": [round(coord, 3) for coord in eq.get('coordinates', [])],
                            "magnitude": eq.get('mag', ''),
                            "depth": round(eq.get('coordinates', [])[2], 1),
                            "alert": eq.get('alert', ''),
                            "url": eq.get('url', '')
                        }

                        # Convert the raw text table to an HTML table
                        html_s1_table = ascii_table_to_html(s1_info)
                        html_nisar_table = ascii_table_to_html(nisar_info)

                        # Construct and send the email
                        header_html = f"""
                        <div style="font-family: Arial, sans-serif; color: #333; max-width: 850px; margin: auto; border: 1px solid #e0e0e0; border-radius: 8px; overflow: hidden;">
                            <div style="background-color: #003366; color: white; padding: 20px;">
                                <h2 style="margin: 0; font-size: 22px;">{message_dict['title']}</h2>
                                <p style="margin: 5px 0 0; font-size: 14px; color: #b3d4fc;">{message_dict['time']} UTC</p>
                            </div>
                            <div style="padding: 20px;">
                                <h3 style="margin: 0 0 10px 0; border-bottom: 2px solid #f0f0f0; padding-bottom: 8px; color: #003366;">Event Details</h3>
                                <table style="width: 100%; text-align: left; margin-bottom: 25px; border-collapse: collapse;">
                                    <tr>
                                        <th style="width: 150px; padding: 4px 0;">Epicenter (Lat, Lon):</th>
                                        <td style="padding: 4px 0;">{message_dict['coordinates'][1]}, {message_dict['coordinates'][0]}</td>
                                    </tr>
                                    <tr>
                                        <th style="padding: 4px 0;">Depth:</th>
                                        <td style="padding: 4px 0;">{message_dict['depth']} km</td>
                                    </tr>
                                </table>
                                <a href="{message_dict['url']}" style="display: inline-block; padding: 10px 18px; background-color: #0055a4; color: white; text-decoration: none; border-radius: 5px; font-weight: bold; margin-bottom: 30px;">View on USGS Hazard Portal</a>
                        """
                        
                        s1_section = f"""
                                <h3 style="margin: 0 0 10px 0; border-bottom: 2px solid #f0f0f0; padding-bottom: 8px; color: #003366;">Sentinel-1 Acquisitions</h3>
                                {html_s1_table}
                        """
                        
                        nisar_section = f"""
                                <h3 style="margin: 25px 0 10px 0; border-bottom: 2px solid #f0f0f0; padding-bottom: 8px; color: #003366;">NISAR Acquisitions</h3>
                                {html_nisar_table}
                        """

                        footer_html = """
                            </div>
                            <div style="background-color: #f9f9f9; padding: 15px; text-align: center; font-size: 12px; color: #888; border-top: 1px solid #e0e0e0;">
                                This is an automated message. Please do not reply.<br>
                                For product-specific inquiries, contact Dr. Cole Speed (<a href="mailto:cole.speed@jpl.nasa.gov">cole.speed@jpl.nasa.gov</a>) and Dr. Grace Bato (<a href="mailto:bato@jpl.nasa.gov">bato@jpl.nasa.gov</a>).
                            </div>
                        </div>
                        """

                        # Helper function to generate the clickable button to route to HTML map
                        def get_button_html(url):
                            if not url: return ""
                            return f"""
                            <div style="text-align: center; margin: 30px 0;">
                                <a href="{url}" style="background-color: #003366; color: white; padding: 15px 40px; text-decoration: none; border-radius: 50px; display: inline-block; font-family: Arial, sans-serif; font-size: 18px;">
                                    Click here for interactive overpass map
                                </a>
                            </div>
                            """

                        if send_email_flag:
                            subject_text = f"New Event: {message_dict['title']}"
                            
                            # Send joint S1 + NISAR email to PRIMARY_RECIPIENTS
                            if config.PRIMARY_RECIPIENTS:
                                primary_body = (header_html + s1_section + nisar_section + get_button_html(map_url) + footer_html).replace('\n', '')
                                send_email(subject=subject_text, body=primary_body, recipients=config.PRIMARY_RECIPIENTS)
                                print('=========================================')
                                print('Joint S1 and NISAR email sent to primary recipients.')
                                print('=========================================')
                        else:
                            print('=========================================')
                            print('Email sending is disabled (--send_email not provided).')
                            print('=========================================')

                        # START TRACKING FOR THIS EVENT
                        # Finds pre-seismic SLCs, creates partial job list, and saves to tracking file
                        add_to_tracker(eq, aoi, resolution)
                        
                else:
                    print(f"No new significant earthquakes found as of {current_time}.")
        else:
            print("Running in --process_only mode. Skipping USGS earthquake discovery.")
        
        # Check ASF DAAC for available SLCs for pending earthquakes
        check_tracker_for_updates(do_processing, send_email_flag)
        
    finally:
        if os.path.exists(lock_file):
            os.remove(lock_file)


def main_historic(start_date=None, end_date=None, eq_list_path=None, aoi=None, pairing_mode=None, job_list=False, resolution=90, sensor='sar', optical_backend='copernicus', optical_level='toa'):
    """
    Runs the main query and processing workflow in historic processing mode.
    Used to produce 'pre-seismic', 'co-seismic', and 'post-seismic' displacement products for historic earthquakes.
    Pre-seismic and post-seismic data are generated for 90 days before and 30 days after the event.
    A single date or date range can be provided for processing.
    :param start_date: The query start date in YYYY-MM-DD format
    :param end_date: The query end date in YYYY-MM-DD format (Optional)
    :param eq_list_path: The path to a custom JSON file containing a list of earthquakes to process (Optional)
    :param aoi: The path to a JSON file representing the area of interest (AOI) (Optional)
    :param pairing_mode: 'all', 'sequential', or 'coseismic' for specifying desired SLC pairing.
    :param job_list: If True, create a list of jobs in HYP3 format for cloud processing.
    :param resolution: Output resolution for the topsApp processing, default is 90m
    :param sensor: 'sar' for SAR processing, 'optical' for optical processing
    :optical_backend: 'copernicus', 'element84' for source data file nomenclature (only applicable if sensor is 'optical')
    :optical_level: 'toa' or 'sr' for the desired optical product level (only applicable if sensor is 'optical')
    """
    # Generate the list of earthquakes
    geojson_data = None
    
    # If a custom earthquake list is provided, use that instead of querying the USGS API, else use the provided dates to query the USGS API for earthquakes in that time range
    if eq_list_path:
        eq_sig = parse_custom_eq_list(eq_list_path)
        if not eq_sig:
            print("No valid earthquakes found in the provided list.")
            return
    else:    
        if start_date and not end_date:
            print('=========================================')
            print(f"Running historic processing in single-date mode for date: {start_date}")
            print('=========================================')
            geojson_data = get_historic_earthquake_data_single_date(USGS_api_alltime, str(start_date))

        elif start_date and end_date:
            print('=========================================')
            print(f"Running historic processing in date range mode for dates: {start_date} to {end_date}")
            print('=========================================')
            geojson_data = get_historic_earthquake_data_date_range(USGS_api_alltime, str(start_date), str(end_date))

        if geojson_data:
            earthquakes = parse_geojson(geojson_data)
            eq_sig = check_significance(earthquakes, start_date, end_date, sensor=sensor, mode='historic')
        else:
            eq_sig = None
            
    # Process the list of earthquakes
    if eq_sig is not None:
        jobs_dict = []
        master_scene_features = []
        earthquake_infos = [] 

        for eq in eq_sig:
            try:
                eq_jsons, eq_features = process_earthquake(eq, aoi, pairing_mode, job_list, resolution, sensor, optical_backend, optical_level)

                if eq_features:
                    master_scene_features.extend(eq_features)
                
                if eq_jsons: 
                    event_dt = convert_time(eq['time'])
                    eq_info = {
                        "title": eq.get('title'),
                        "epicenter": {
                            "latitude": eq['coordinates'][1],
                            "longitude": eq['coordinates'][0],
                            "depth_km": eq['coordinates'][2],
                            "usgs_event_url": eq.get('url', '')
                        },
                        "time": event_dt.strftime("%Y-%m-%d %H:%M:%S UTC")
                    }
                    earthquake_infos.append(eq_info)

            except Exception as e:
                print(f"Error processing {eq['title']}: {e}")
                continue

            if eq_jsons:
                for i, eq_json in enumerate(eq_jsons):
                    for j, json_data in enumerate(eq_json):
                        jobs_dict.append(json_data)

        if jobs_dict:
            current_time = datetime.now(timezone.utc).strftime("%Y-%m-%d_%H-%M-%S_UTC")
            
            with open(f'jobs_list_{current_time}.json', 'w') as f:
                json.dump(jobs_dict, f, indent=4)
            
            with open(f'earthquake_info_{current_time}.json', 'w', encoding='utf-8') as f:
                json.dump(earthquake_infos, f, indent=4, ensure_ascii=False)

        if master_scene_features:
            current_time = datetime.now(timezone.utc).strftime("%Y-%m-%d_%H-%M-%S")
            feature_filename = f'all_selected_scenes_{sensor}_{optical_backend}_{current_time}.geojson'
            fc = {"type": "FeatureCollection", "features": master_scene_features}
            with open(feature_filename, 'w') as f:
                json.dump(fc, f, indent=2)
            print(f"Saved master scene footprints to {feature_filename}")
    else:
        if eq_list_path:
            print("No significant earthquakes found in the provided list.")
        else:
            print(f"No significant earthquakes found between {start_date} and {end_date}.")
