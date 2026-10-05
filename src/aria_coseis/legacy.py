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

def load_tracker():
    """Loads all active jobs from the tracking directory."""
    tracker = {}
    if not os.path.exists(config.TRACKING_DIR):
        os.makedirs(config.TRACKING_DIR, exist_ok=True)
        return tracker
    
    for file in glob.glob(os.path.join(config.TRACKING_DIR, "*.json")):
        try:
            with open(file, "r") as f:
                event_data = json.load(f)
                # The filename (minus .json) is the event_id
                event_id = os.path.basename(file).replace('.json', '')
                tracker[event_id] = event_data
        except json.JSONDecodeError:
            continue
    return tracker

def save_tracker(data):
    """Saves the tracking data back to individual files."""
    if not os.path.exists(config.TRACKING_DIR):
        os.makedirs(config.TRACKING_DIR, exist_ok=True)
        
    # Save current tracker state to individual files
    for event_id, event_data in data.items():
        file_path = os.path.join(config.TRACKING_DIR, f"{event_id}.json")
        with open(file_path, "w") as f:
            json.dump(event_data, f, indent=4)
            
    # Remove files for events that are no longer in the tracker dictionary
    for file in glob.glob(os.path.join(config.TRACKING_DIR, "*.json")):
        event_id = os.path.basename(file).replace('.json', '')
        if event_id not in data:
            os.remove(file)


def add_to_tracker(eq, aoi, resolution=30):
    """
    Initializes tracking for a new earthquake.
    Identifies intersecting tracks, finds pre-seismic SLCs for each track, 
    creates a partial job file (granules=Empty, secondary_granules=Filled), and adds entry to tracking file.
    :param eq: dictionary containing earthquake data
    :param aoi: shapely Polygon object representing the area of interest
    :param resolution: desired output resolution for processing (default=30)
    """
    tracker = load_tracker()
    event_id = eq.get('id')
    title = to_snake_case(eq.get('title'))
    event_time = eq.get('time')
    
    if event_id in tracker:
        print(f"Event {title} is already being tracked.")
        return

    print('=========================================')
    print(f"Initializing tracking for {title}...")
    print('=========================================')

    # Get intersecting tracks
    path_frame_numbers, _ = get_path_and_frame_numbers(aoi, event_time)
    
    tracks_info = {}

    for (flight_direction, path_number), frame_numbers_set in path_frame_numbers.items():
        frame_numbers = list(set(fn[0] for fn in frame_numbers_set))
        
        # Unique key for this track
        track_key = f"{flight_direction}_{path_number}"
        
        # Fetch SLCs intersecting the AOI (Historical search relative to event time)
        slcs = get_SLCs(flight_direction, path_number, aoi.wkt, event_time, processing_mode='historic')
        
        rupture_dt = convert_time(event_time).replace(tzinfo=None)
        pre_slcs = []
        reference_date = None
        
        if slcs:
            # Determine the maximum footprint this specific track has over the AOI
            scenes_by_date = defaultdict(list)
            for s in slcs:
                scenes_by_date[s['date'][:10]].append(s)
                
            max_track_area = 0
            for date_str, scenes in scenes_by_date.items():
                union_geom = unary_union([s['geometry'] for s in scenes])
                intersection_area = union_geom.intersection(aoi).area
                if intersection_area > max_track_area:
                    max_track_area = intersection_area

            # Filter for pre-seismic scenes and select the closest valid date
            valid_pre_scenes = [s for s in slcs if datetime.strptime(s['date'], "%Y-%m-%dT%H:%M:%SZ") < rupture_dt]
            
            if valid_pre_scenes:
                # Sort dates newest to oldest (closest to earthquake first)
                dates = sorted(list(set(s['date'][:10] for s in valid_pre_scenes)), reverse=True)
                
                for ref_date in dates:
                    candidate_scenes = [s for s in valid_pre_scenes if s['date'][:10] == ref_date]
                    union_geom = unary_union([s['geometry'] for s in candidate_scenes])
                    intersection_area = union_geom.intersection(aoi).area
                    
                    # Ensure coverage is at least 95% of the track's max expected footprint
                    if max_track_area > 0 and (intersection_area / max_track_area) > 0.95:
                        pre_slcs = [s['fileID'].removesuffix("-SLC") for s in candidate_scenes]
                        reference_date = ref_date
                        break # Successfully found valid coverage
        
        if not pre_slcs:
             print(f"No fully overlapping pre-seismic coverage found for {track_key}. Skipping track.")
             continue

        # Create Partial Job List
        # We leave 'granules' (post-seismic) empty for now
        # We fill 'secondary_granules' (pre-seismic)
        job_filename = f"job_{title}_{track_key}_partial.json"
        
        # Create the standard HYP3 structure
        job_json = make_job_json(title, event_id, flight_direction, path_number, [], pre_slcs, resolution)
        
        # Save Partial File
        with open(job_filename, "w") as f:
            json.dump([job_json], f, indent=4) # List of 1 job

        # Add to tracks info
        tracks_info[track_key] = {
            "flight_direction": flight_direction,
            "path_number": path_number,
            "frame_numbers": frame_numbers,
            "partial_job_file": job_filename,
            "reference_date": reference_date,
            "status": "AWAITING_POST_SEISMIC"
        }
        print(f"  Initialized track {track_key}. Pre-seismic date: {reference_date}. Waiting for post-seismic.")

    if tracks_info:
        tracker[event_id] = {
            "title": title,
            "time": event_time,
            "aoi": mapping(aoi),
            "tracks": tracks_info
        }
        save_tracker(tracker)
        print(f"Added {title} to tracking file with {len(tracks_info)} tracks.")
    else:
        print(f"No valid tracks initialized for {title}.")


def check_tracker_for_updates(do_processing=False, send_email_flag=False):
    """
    Iterates through the tracking file using a unified state machine.
    - Local Machine (do_processing=True): Looks for 'AWAITING_POST_SEISMIC', runs topsApp, updates state to 'READY_FOR_EMAIL'.
    - GitHub Action (send_email_flag=True): Looks for 'READY_FOR_EMAIL', sends email, removes event.

    Args:
        do_processing (bool): If True, this machine will perform local processing for tracks awaiting post-seismic data.
        send_email_flag (bool): If True, this machine will send emails for tracks that are ready for email.
    """
    tracker = load_tracker()
    if not tracker:
        return

    events_to_remove = []
    completed_jobs_summary = []

    for event_id, event_data in tracker.items():
        title = event_data['title']
        event_time = event_data['time']
        tracks = event_data['tracks']
        
        tracks_to_remove = []

        for track_key, track_info in tracks.items():
            
            # Local Processing
            if track_info['status'] == "AWAITING_POST_SEISMIC":
                if not do_processing:
                    continue 

                flight_dir = track_info['flight_direction']
                path_num = track_info['path_number']
                
                # Reconstruct the AOI geometry from the tracker dictionary
                aoi_geom = shape(event_data['aoi'])
                
                slcs = get_SLCs(flight_dir, path_num, aoi_geom.wkt, event_time, processing_mode='forward')
                
                rupture_dt = convert_time(event_time).replace(tzinfo=None)
                post_slcs = []
                secondary_date = None

                if slcs:
                    # Determine the maximum footprint this specific track has over the AOI
                    scenes_by_date = defaultdict(list)
                    for s in slcs:
                        scenes_by_date[s['date'][:10]].append(s)
                        
                    max_track_area = 0
                    for date_str, scenes in scenes_by_date.items():
                        union_geom = unary_union([s['geometry'] for s in scenes])
                        intersection_area = union_geom.intersection(aoi_geom).area
                        if intersection_area > max_track_area:
                            max_track_area = intersection_area

                    # Filter for post-seismic scenes and select the closest valid date
                    slcs.sort(key=lambda x: x['date'])
                    valid_post_scenes = [s for s in slcs if datetime.strptime(s['date'], "%Y-%m-%dT%H:%M:%SZ") > rupture_dt]                    
                    
                    if valid_post_scenes:
                        # Sort dates oldest to newest (closest to earthquake first)
                        post_dates = sorted(list(set(s['date'][:10] for s in valid_post_scenes)))
                        
                        for sec_date in post_dates:
                            candidate_scenes = [s for s in valid_post_scenes if s['date'][:10] == sec_date]
                            union_geom = unary_union([s['geometry'] for s in candidate_scenes])
                            intersection_area = union_geom.intersection(aoi_geom).area
                            
                            # Ensure coverage is at least 95% of the track's max expected footprint
                            if max_track_area > 0 and (intersection_area / max_track_area) > 0.95:
                                post_slcs = [s['fileID'].removesuffix("-SLC") for s in candidate_scenes]
                                secondary_date = sec_date
                                break # Successfully found valid coverage
                
                if post_slcs:
                    pre_seismic_date = track_info['reference_date']
                    post_seismic_date = secondary_date
                    
                    partial_file = track_info['partial_job_file']
                    try:
                        with open(partial_file, 'r') as f:
                            job_list = json.load(f)
                            job = job_list[0]
                        
                        job['job_parameters']['granules'] = post_slcs

                        older_date_str = pre_seismic_date.split('T')[0].replace("-", "")
                        newer_date_str = post_seismic_date.split('T')[0].replace("-", "")
                        pair_folder_name = f"{flight_dir}{int(path_num):03d}_{older_date_str}_{newer_date_str}"
                        processing_dir = os.path.join(
                            config.root_dir, title, f"{flight_dir}{int(path_num):03d}", "coseismic", pair_folder_name
                        )
                        
                        print(f"    Starting automatic processing for {pair_folder_name}")
                        try:
                            run_dockerized_topsApp(job, processing_dir)
                            track_info['processing_status'] = "Success"

                            # Delete raw SLCs and intermediate files if final .nc product exists
                            nc_files = glob.glob(os.path.join(processing_dir, "**", "*.nc"), recursive=True)
                            
                            if nc_files:
                                print(f"    Output .nc file found. Cleaning up SLCs and heavy intermediate files in {pair_folder_name}...")
                                import shutil
                                
                                # Delete raw .zip and .SAFE files
                                for slc_zip in glob.glob(os.path.join(processing_dir, "S1[A-D]*.zip")):
                                    os.remove(slc_zip)
                                for slc_safe in glob.glob(os.path.join(processing_dir, "S1[A-D]*.SAFE")):
                                    shutil.rmtree(slc_safe, ignore_errors=True)
                                    
                                # Delete intermediate ISCE2 folders
                                intermediate_dirs = [
                                    "geom_reference", "ion", "fine_interferogram", 
                                    "fine_offsets", "fine_coreg", "mask", 
                                    "reference", "secondary", "PICKLE", "aux_cal", "orbits"
                                ]
                                for idir in intermediate_dirs:
                                    dir_path = os.path.join(processing_dir, idir)
                                    if os.path.exists(dir_path):
                                        shutil.rmtree(dir_path, ignore_errors=True)
                            else:
                                print(f"    WARNING: No final .nc product found in {pair_folder_name}. Retaining raw and intermediate files for debugging.")

                        except Exception as e:
                            print(f"    Processing failed for {pair_folder_name}: {e}")
                            track_info['processing_status'] = f"Failed: {str(e)}"

                        completed_filename = os.path.join(processing_dir, f"job_{title}_{track_key}_COMPLETED.json")
                        with open(completed_filename, 'w') as f:
                            json.dump([job], f, indent=4)
                        
                        # Notify GitHub Actions that the job is done and ready for email
                        track_info['status'] = "READY_FOR_EMAIL"
                        track_info['dates'] = f"{pre_seismic_date} (Pre) - {post_seismic_date} (Post)"
                        track_info['location'] = processing_dir
                        
                        if os.path.exists(partial_file):
                            os.remove(partial_file)

                    except Exception as e:
                        print(f"    Error processing partial file {partial_file}: {e}")
                else:
                    print(f"    No post-seismic data yet.")
                    
            # Email with Github Actions
            elif track_info['status'] == "READY_FOR_EMAIL":
                if not send_email_flag:
                    continue 
                    
                completed_jobs_summary.append({
                    "title": title,
                    "track": track_key,
                    "dates": track_info.get('dates', 'Unknown'),
                    "status": track_info.get('processing_status', 'Unknown'),
                    "location": track_info.get('location', 'Unknown')
                })
                
                # Determine whether to delete or quarantine based on success/failure
                if track_info.get('processing_status', '').startswith("Failed"):
                    print(f"    Job {track_key} failed. Moving to FAILED_NEEDS_ATTENTION state.")
                    track_info['status'] = "FAILED_NEEDS_ATTENTION"
                else:
                    # Mark successful track for removal
                    tracks_to_remove.append(track_key)

        # Remove completed tracks
        for tk in tracks_to_remove:
            del tracks[tk]
        
        # If no tracks left, mark event for removal
        if not tracks:
            events_to_remove.append(event_id)

    # Clean up tracker
    for eid in events_to_remove:
        print(f"Event {tracker[eid]['title']} fully processed. Removing from tracker.")
        del tracker[eid]
    
    save_tracker(tracker)

    # Send "PROCESSING COMPLETED" email
    if completed_jobs_summary and send_email_flag:
        has_failures = any(item['status'].startswith("Failed") for item in completed_jobs_summary)
        status_tag = "WITH FAILURES" if has_failures else "SUCCESS"
        subject = f"PROCESSING COMPLETED ({status_tag}): {len(completed_jobs_summary)} Jobs Processed"
        
        # Start HTML body with inline CSS for universal email client support
        body = """
        <html>
            <body style="font-family: Arial, sans-serif; color: #333; line-height: 1.6; margin: 0; padding: 20px;">
            <div style="max-width: 600px; margin: 0 auto; border: 1px solid #e0e0e0; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05);">
                <h2 style="color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; margin-top: 0;">Automated Processing Update</h2>
                <p>The following SAR processing jobs have concluded:</p>
        """
        
        # Generate a styled card for each completed job
        for item in completed_jobs_summary:
            status_color = "#e74c3c" if item['status'].startswith("Failed") else "#27ae60"
            
            body += f"""
                <div style="background-color: #f8f9fa; padding: 15px; margin-bottom: 15px; border-left: 5px solid {status_color}; border-radius: 4px;">
                <p style="margin: 0 0 5px;"><strong>Event:</strong> {item['title']}</p>
                <p style="margin: 0 0 5px;"><strong>Track:</strong> {item['track']}</p>
                <p style="margin: 0 0 5px;"><strong>Dates:</strong> {item['dates']}</p>
                <p style="margin: 0 0 5px;"><strong>Status:</strong> <span style="color: {status_color}; font-weight: bold;">{item['status']}</span></p>
                <p style="margin: 0;"><strong>Location:</strong> <br>
                    <code style="background: #e9ecef; padding: 4px; display: block; margin-top: 5px; word-wrap: break-word; font-size: 12px; color: #c0392b;">
                    {item['location']}
                    </code>
                </p>
                </div>
            """

        # Close the HTML body
        body += """
                <p style="font-size: 12px; color: #95a5a6; border-top: 1px solid #e0e0e0; padding-top: 15px; margin-top: 20px; text-align: center;">
                This is an automated message from the COSEIS pipeline.
                </p>
            </div>
            </body>
        </html>
        """
        
        print("Sending completion email to secondary recipients...")
        send_email(subject, body, recipients=config.SECONDARY_RECIPIENTS)


def export_gee_sentinel2_composite(aoi_polygon, start_date, end_date, title, stage, gcs_bucket, collection_id, optical_level, crs_epsg='EPSG:4326'):
    """
    Generates a cloud-free median composite in GEE and exports to Google Cloud Storage.
    :param aoi_polygon: Shapely Polygon representing the Area of Interest
    :param start_date: Start date for the image collection filter (YYYY-MM-DD)
    :param end_date: End date for the image collection filter (YYYY-MM-DD)
    :param title: Title of the earthquake event (used for file naming)
    :param stage: 'pre-event' or 'post-event' to indicate the timing of the composite
    :param gcs_bucket: Name of the Google Cloud Storage bucket to export the composite
    :param collection_id: Sentinel-2 collection ID (e.g., 'COPERNICUS/S2_HARMONIZED')
    :param optical_level: Optical processing level (e.g., TOA, SR) for file naming
    :param crs_epsg: EPSG code for the coordinate reference system to use in the export (default is 'EPSG:4326')
    :return: The GEE export task object
    """
    import ee

    # Convert Shapely polygon to GEE Geometry
    bounds = aoi_polygon.bounds
    ee_roi = ee.Geometry.Rectangle([bounds[0], bounds[1], bounds[2], bounds[3]])

    # Query the parameterized Harmonized Collection
    s2_col = ee.ImageCollection(collection_id) \
        .filterBounds(ee_roi) \
        .filterDate(start_date, end_date)

    # Apply Cloud Score+ Masking
    cs_plus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')
    s2_linked = s2_col.linkCollection(cs_plus, ['cs_cdf'])

    def mask_clouds(img):
        mask = img.select('cs_cdf').gte(0.65)
        return img.updateMask(mask)

    s2_masked = s2_linked.map(mask_clouds)

    # Fetch distinct orbits intersecting the AOI to the local Python environment
    distinct_orbits = ee.List(s2_masked.aggregate_array('SENSING_ORBIT_NUMBER')).distinct().getInfo()
    
    if not distinct_orbits:
        print(f"  No Sentinel-2 data found for {stage} between {start_date} and {end_date}.")
        return {}, []

    print(f"  Found {len(distinct_orbits)} distinct SENSING_ORBIT_NUMBER(s): {distinct_orbits}")
    
    orbit_exports = {}
    all_unique_dates = []

    # Iterate over each distinct orbit
    for orbit in distinct_orbits:
        orbit_col = s2_masked.filter(ee.Filter.eq('SENSING_ORBIT_NUMBER', orbit))
        
        def get_date(img):
            return ee.Feature(None, {'date': img.date().format("yyyy-MM-dd'T'HH:mm:ss")})
        
        raw_dates = orbit_col.map(get_date).aggregate_array('date').getInfo()
        orbit_dates = sorted(list(set(raw_dates)))
        all_unique_dates.extend(orbit_dates)
        
        # Calculate median for JUST this orbit
        composite = orbit_col.median().select('B8').toUint16().clip(ee_roi)

        # Clean up the start and end dates for file naming
        clean_start = start_date.replace('T', '_').replace(':', '') + 'UTC' if 'T' in start_date else start_date + '_000000UTC'
        clean_end = end_date.replace('T', '_').replace(':', '') + 'UTC' if 'T' in end_date else end_date + '_000000UTC'

        # Create Export Task
        file_name = f"{title}_S2_Path{orbit}_{optical_level.upper()}_B8_{stage}_{clean_start}_to_{clean_end}"
        prefix = f"COSEIS_Composites/{title}/{file_name}"
        
        task = ee.batch.Export.image.toCloudStorage(
            image=composite,
            description=file_name[:100],
            bucket=gcs_bucket,
            fileNamePrefix=prefix,
            region=ee_roi,
            scale=10, 
            crs=crs_epsg,
            maxPixels=1e13,
            formatOptions={'cloudOptimized': False}
        )
        
        task.start()
        print(f"  Started GCS Export Task for Orbit {orbit}: {file_name}")
        
        orbit_exports[str(orbit)] = {
            'task': task,
            'gcs_uri': f"gs://{gcs_bucket}/{prefix}.tif",
            'prefix': prefix,
            'dates_included': orbit_dates
        }

    return orbit_exports, sorted(list(set(all_unique_dates)))


def export_gee_landsat_composite(aoi_polygon, start_date, end_date, title, stage, gcs_bucket, collection_id, band_name, scale, optical_level, crs_epsg='EPSG:4326'):
    """Generates a cloud-free median TOA composite for the explicitly provided Landsat mission.
    :param aoi_polygon: Shapely Polygon representing the Area of Interest
    :param start_date: Start date for the image collection filter (YYYY-MM-DD)
    :param end_date: End date for the image collection filter (YYYY-MM-DD)
    :param title: Title of the earthquake event (used for file naming)
    :param stage: 'pre-event' or 'post-event' to indicate the timing of the composite
    :param gcs_bucket: Name of the Google Cloud Storage bucket to export the composite
    :param collection_id: Landsat collection ID (e.g., 'LANDSAT/LC08/C02/T1_L2')
    :param band_name: Name of the optical band to export (e.g., 'SR_B4' for Landsat 8 Red)
    :param scale: Scale in meters for the export (e.g., 30 for Landsat)
    :param optical_level: Optical processing level (e.g., TOA, SR, RAW) for file naming
    :param crs_epsg: EPSG code for the coordinate reference system to use in the export (default is 'EPSG:4326')
    :return: The GEE export task object and a list of unique acquisition dates
    """
    import ee

    bounds = aoi_polygon.bounds
    ee_roi = ee.Geometry.Rectangle([bounds[0], bounds[1], bounds[2], bounds[3]])

    # Determine Landsat Mission from collection_id for dynamic file naming
    if 'LT05' in collection_id:
        mission_name = 'Landsat5'
    elif 'LE07' in collection_id:
        mission_name = 'Landsat7'
    elif 'LC08' in collection_id:
        mission_name = 'Landsat8'
    elif 'LC09' in collection_id:
        mission_name = 'Landsat9'
    else:
        mission_name = 'Landsat'

    # Query the explicit collection passed into the function
    l_col = ee.ImageCollection(collection_id) \
        .filterBounds(ee_roi) \
        .filterDate(start_date, end_date)

    # Apply USGS QA_PIXEL bitmask for clouds and shadows
    def mask_clouds(img):
        qa = img.select('QA_PIXEL')
        cloud_shadow_bitmask = (1 << 4)
        clouds_bitmask = (1 << 3)
        dilated_cloud_bitmask = (1 << 1)
        
        mask = qa.bitwiseAnd(cloud_shadow_bitmask).eq(0) \
            .And(qa.bitwiseAnd(clouds_bitmask).eq(0)) \
            .And(qa.bitwiseAnd(dilated_cloud_bitmask).eq(0))
            
        return img.updateMask(mask).select([band_name], ['OPTICAL_BAND'])

    l_masked = l_col.map(mask_clouds)

    # Fetch distinct paths intersecting the AOI to the local Python environment
    # .getInfo() pulls the list from GEE servers to local execution
    distinct_paths = ee.List(l_masked.aggregate_array('WRS_PATH')).distinct().getInfo()
    
    if not distinct_paths:
        print(f"  No Landsat data found for {stage} between {start_date} and {end_date}.")
        return {}, []

    print(f"  Found {len(distinct_paths)} distinct WRS_PATH(s): {distinct_paths}")
    
    path_exports = {}
    all_unique_dates = []

    # Iterate over each distinct path to create independent composites and export tasks
    for path in distinct_paths:
        path_col = l_masked.filter(ee.Filter.eq('WRS_PATH', path))
        
        # Map over the specific path collection to get acquisition dates
        def get_date(img):
            return ee.Feature(None, {'date': img.date().format("yyyy-MM-dd'T'HH:mm:ss")})
        
        raw_dates = path_col.map(get_date).aggregate_array('date').getInfo()
        path_dates = sorted(list(set(raw_dates)))
        all_unique_dates.extend(path_dates)
        
        # Calculate median for JUST this path to avoid cross-track blending
        composite_raw = path_col.median().select('OPTICAL_BAND').clip(ee_roi)
        
        # Scale based on the collection type to ensure autoRIFT gets clean uint16
        if 'TOA' in collection_id:
            composite = composite_raw.multiply(10000).toUint16()
        elif '_L2' in collection_id:
            # GEE C02 L2 Surface Reflectance scaling: (pixel * 0.0000275) - 0.2
            composite = composite_raw.multiply(0.0000275).subtract(0.2).multiply(10000).max(0).toUint16()
        else:
            composite = composite_raw.unmask(0).toUint16()
        
        # Clean up the start and end dates for file naming
        clean_start = start_date.replace('T', '_').replace(':', '') + 'UTC' if 'T' in start_date else start_date + '_000000UTC'
        clean_end = end_date.replace('T', '_').replace(':', '') + 'UTC' if 'T' in end_date else end_date + '_000000UTC'

        # Include the specific WRS_PATH and optical_level in the file name
        file_name = f"{title}_{mission_name}_{optical_level.upper()}_Path{path}_{band_name}_{stage}_{clean_start}_to_{clean_end}"
        prefix = f"COSEIS_Composites/{title}/{file_name}"
        
        task = ee.batch.Export.image.toCloudStorage(
            image=composite,
            description=file_name[:100],
            bucket=gcs_bucket,
            fileNamePrefix=prefix,
            region=ee_roi,
            scale=scale, 
            crs=crs_epsg,
            maxPixels=1e13,
            formatOptions={'cloudOptimized': False}
        )
        
        task.start()
        print(f"  Started GCS Export Task for Path {path}: {file_name}")
        
        # Store the task and metadata keyed by path string to build the manifest later
        path_exports[str(path)] = {
            'task': task,
            'gcs_uri': f"gs://{gcs_bucket}/{prefix}.tif",
            'prefix': prefix,
            'dates_included': path_dates
        }

    return path_exports, sorted(list(set(all_unique_dates)))


def wait_for_gee_tasks(tasks, timeout_mins=60):
    """
    Polls GEE until all provided tasks are either COMPLETED or FAILED.
    Includes a timeout to prevent infinite hangs on stuck GEE backend tasks.
    :param tasks: List of GEE export task objects to monitor
    :param timeout_mins: Maximum minutes to wait before canceling tasks.
    """
    print(f"Waiting for Google Earth Engine exports to complete (Timeout: {timeout_mins} mins)...")
    start_time = time.time()
    timeout_seconds = timeout_mins * 60

    for task in tasks:
        while task.active():
            elapsed = time.time() - start_time
            if elapsed > timeout_seconds:
                print(f"  Task {task.id} TIMED OUT after {timeout_mins} minutes. Canceling task.")
                try:
                    task.cancel()
                except Exception as e:
                    print(f"  Could not cancel task {task.id}: {e}")
                break
                
            print(f"  Task {task.id} is {task.status()['state']}... waiting 30 seconds.")
            time.sleep(30)
        
        status = task.status()
        if status['state'] == 'COMPLETED':
            print(f"  Task {task.id} COMPLETED.")
        elif status['state'] in ['CANCELLED', 'CANCELED']:
            print(f"  Task {task.id} CANCELLED due to timeout.")
        else:
            print(f"  Task {task.id} FAILED: {status.get('error_message', 'Unknown error')}")


def download_from_gcs(bucket_name, prefix, local_dir):
    """Downloads a blob from the GCS bucket to user's local machine.
    :param bucket_name: Name of the GCS bucket
    :param prefix: Prefix of the blob in the GCS bucket (including path) to identify the file(s) to download
    :param local_dir: Local directory where the file(s) should be downloaded
    """
    from google.cloud import storage

    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blobs = bucket.list_blobs(prefix=prefix)

    os.makedirs(local_dir, exist_ok=True)
    downloaded_files = []

    for blob in blobs:
        # Extract the actual filename GEE generated (including any chip suffixes)
        filename = os.path.basename(blob.name)
        dest_path = os.path.join(local_dir, filename)
        
        blob.download_to_filename(dest_path)
        print(f"  Successfully downloaded chip: {dest_path}")
        downloaded_files.append(dest_path)
        
    return downloaded_files


def merge_and_compress_chips(chip_paths, output_path):
    """Uses GDAL to merge GEE image chips into a single compressed GeoTIFF.
    :param chip_paths: List of file paths to the individual image chips downloaded from GCS
    :param output_path: Desired file path for the final merged and compressed GeoTIFF
    """
    print(f"  Stitching and compressing {len(chip_paths)} file(s)...")
    
    # Create temporary output name to prevent GDAL from overwriting input
    tmp_output = output_path + ".tmp.tif"
    
    # gdalwarp merges the chips and applies DEFLATE compression
    cmd = [
        "gdalwarp",
        "-co", "COMPRESS=DEFLATE",
        "-co", "PREDICTOR=2",
        "-co", "TILED=YES",
        "-co", "NUM_THREADS=ALL_CPUS",
        *chip_paths,
        tmp_output
    ]
    
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Delete the individual chips
    for chip in chip_paths:
        try:
            os.remove(chip)
        except OSError:
            pass
            
    # Rename the clean, temporary file to your final desired filename
    os.rename(tmp_output, output_path)
        
    print(f"  Successfully created unified composite: {output_path}")
    return output_path


def assign_nodata(filepath, nodata_val=0):
    """
    Opens the specified GeoTIFF and explicitly writes the NoData 
    value into the metadata header.
    """
    from osgeo import gdal

    # Open the file in update mode, explicitly passing the open option to break COG layout
    ds = gdal.OpenEx(filepath, gdal.OF_UPDATE, open_options=["IGNORE_COG_LAYOUT_BREAK=YES"])
    
    if ds is not None:
        band = ds.GetRasterBand(1)
        if band.GetNoDataValue() is None:
            band.SetNoDataValue(nodata_val)
        # FlushCache ensures the metadata is written securely to disk
        ds.FlushCache()
        ds = None
    else:
        print(f"Warning: Could not open {filepath} to assign NoData value.")


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
