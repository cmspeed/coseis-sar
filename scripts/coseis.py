import re
import unicodedata
import os
import argparse
import asf_search as asf
from dateutil import parser as dateparser
from dateutil.parser import isoparse
import ee
from pathlib import Path
from pystac_client import Client
import requests
import json
import folium
import geojson
import geopandas as gpd
import csv
from shapely import wkt
from shapely.geometry import mapping, shape, box, Point, Polygon, LineString, MultiLineString, MultiPolygon
from shapely.ops import unary_union
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from collections import defaultdict
from itertools import combinations
import logging
import yagmail
import time
from datetime import timedelta
from google.cloud import storage
from time import sleep
from types import SimpleNamespace
from urllib.parse import urlparse
from typing import List, Dict, Any, Optional

# Set logging level to WARNING to suppress DEBUG and INFO logs
logging.basicConfig(level=logging.WARNING)

# API endpoints
USGS_api_hourly = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_hour.geojson"  # USGS Earthquake API - Hourly
USGS_api_daily = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_day.geojson"  # USGS Earthquake API - Daily
USGS_api_30day = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_month.geojson"  # USGS Earthquake API - Monthly
USGS_api_alltime = "https://earthquake.usgs.gov/fdsnws/event/1/query" # USGS Earthquake API - All Time
coastline_api = "https://raw.githubusercontent.com/OSGeo/PROJ/refs/heads/master/docs/plot/data/coastline.geojson" # Coastline API
ASF_DAAC_API = "https://api.daac.asf.alaska.edu/services/search/param" # ASF DAAC API endpoint
CMR_API_URL = "https://cmr.earthdata.nasa.gov/search/granules.json" # NASA CMR API endpoint
root_dir = os.path.join(os.getcwd(), "data")

# Global variables
OPTICAL_CLOUD_THRESHOLD = 20.0  # Maximum cloud cover percentage for optical data

TRACKING_FILE = "active_job_tracking.json"

def get_recipients_from_env(var_name):
    """
    Retrieves a list of emails from an environment variable.
    """
    env_val = os.getenv(var_name, "")
    # Split by comma and strip whitespace
    return [email.strip() for email in env_val.split(',') if email.strip()]

# Load recipients from environment variables
PRIMARY_RECIPIENTS = get_recipients_from_env('COSEIS_PRIMARY_RECIPIENTS')
SECONDARY_RECIPIENTS = get_recipients_from_env('COSEIS_SECONDARY_RECIPIENTS')
TERTIARY_RECIPIENTS = get_recipients_from_env('COSEIS_TERTIARY_RECIPIENTS')

def load_tracker():
    """Loads the active job tracking file."""
    if not os.path.exists(TRACKING_FILE):
        return {}
    try:
        with open(TRACKING_FILE, "r") as f:
            return json.load(f)
    except json.JSONDecodeError:
        return {}


def save_tracker(data):
    """Saves the active job tracking file."""
    with open(TRACKING_FILE, "w") as f:
        json.dump(data, f, indent=4)


def add_to_tracker(eq, aoi, resolution=90):
    """
    Initializes tracking for a new earthquake.
    Identifies intersecting tracks, finds pre-seismic SLCs for each track, 
    creates a partial job file (granules=Empty, secondary_granules=Filled), and adds entry to tracking file.
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
        
        # Find Pre-seismic SLCs (Historical search relative to event time)
        # We look back 24 days to ensure we get the latest coverage
        slcs = get_SLCs(flight_direction, path_number, frame_numbers, event_time, processing_mode='historic')
        
        # Filter for pre-seismic only (closest to event)
        rupture_dt = convert_time(event_time).replace(tzinfo=None)
        pre_slcs = []
        if slcs:
            # Sort by date
            slcs.sort(key=lambda x: x['date'])
            # Find the SLCs immediately preceding the rupture
            # We want the single latest acquisition date before rupture
            dates = sorted(list(set(s['date'] for s in slcs)))
            pre_dates = [d for d in dates if datetime.strptime(d[:10], "%Y-%m-%d") < rupture_dt]
            
            if pre_dates:
                reference_date = pre_dates[-1] # The last date before the earthquake
                pre_slcs = [s['fileID'].removesuffix("-SLC") for s in slcs if s['date'] == reference_date]
            else:
                print(f"No pre-seismic data found for {track_key}. Skipping track.")
                continue
        else:
             print(f"No SLCs found for {track_key}. Skipping track.")
             continue

        if not pre_slcs:
            continue

        # Create Partial Job List
        # We leave 'granules' (post-seismic) empty for now
        # We fill 'secondary_granules' (pre-seismic)
        job_filename = f"job_{title}_{track_key}_partial.json"
        
        # Create the standard HYP3 structure
        job_json = make_job_json(title, flight_direction, path_number, [], pre_slcs, resolution)
        
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
                frames = track_info['frame_numbers']
                
                slcs = get_SLCs(flight_dir, path_num, frames, event_time, processing_mode='forward')
                
                rupture_dt = convert_time(event_time).replace(tzinfo=None)
                post_slcs = []
                secondary_date = None

                if slcs:
                    slcs.sort(key=lambda x: x['date'])
                    post_dates = sorted(list(set(s['date'] for s in slcs if datetime.strptime(s['date'][:10], "%Y-%m-%d") > rupture_dt)))                    
                    if post_dates:
                        secondary_date = post_dates[0]
                        post_slcs = [s['fileID'].removesuffix("-SLC") for s in slcs if s['date'] == secondary_date]
                
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
                            root_dir, title, f"{flight_dir}{int(path_num):03d}", "coseismic", pair_folder_name
                        )
                        
                        print(f"    Starting automatic processing for {pair_folder_name}")
                        try:
                            run_dockerized_topsApp(job, processing_dir)
                            track_info['processing_status'] = "Success"
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
        send_email(subject, body, recipients=SECONDARY_RECIPIENTS)


def get_historic_earthquake_data_single_date(eq_api, input_date):
    """
    Fetch data from the USGS Earthquake Portal for a single date and returns it as a GeoJSON object.
    The data returned will depend on the parameters included with the API request.
    :param eq_api: USGS API endpoint
    :param input_date: date in the format 'YYYY-MM-DD'
    :return: GeoJSON object containing earthquake data
    """
    print('=========================================')
    print(f"Fetching historic earthquake data from {input_date}...")
    print('=========================================')
    try:

        # Parameters for the API request
        params = {
            "format": "geojson",
            "starttime": input_date + "00:00:00",
            "endtime": input_date + "23:59:59",
            "minmagnitude": 6.0,
            "maxdepth": 30.0
        }

        # Fetch data from the USGS Earthquake API
        response = requests.get(eq_api, params=params)
        response.raise_for_status()  # Raise error if request fails
        
        # Parse the response as GeoJSON
        earthquakes = geojson.loads(response.text)

        return earthquakes

    except requests.RequestException as e:
        print(f"Error accessing primary API: {e}")
        return None
    except geojson.GeoJSONDecodeError as e:
        print(f"Error parsing GeoJSON data: {e}")
        return None


def get_historic_earthquake_data_date_range(eq_api, start_date, end_date):
    """
    Fetch data from the USGS Earthquake Portal over the date range and returns it as a GeoJSON object.
    The data returned will depend on the parameters included with the API request.
    :param eq_api: USGS API endpoint
    :param start_date: start date in the format 'YYYY-MM-DD'
    :param end_date: end date in the format 'YYYY-MM-DD'
    :return: GeoJSON object containing earthquake data over the date range requested
    """
    start_date = start_date + "T00:00:00"
    end_date = end_date + "T23:59:59"

    print('=========================================')
    print(f"Fetching historic earthquake data from {start_date} to {end_date}...")
    print('=========================================')
    try:

        # Parameters for the API request
        params = {
            "format": "geojson",
            "starttime": start_date,
            "endtime": end_date,
            "minmagnitude": 5.0,
            "maxdepth": 40.0
        }

        # Fetch data from the USGS Earthquake API
        response = requests.get(eq_api, params=params)
        response.raise_for_status()  # Raise error if request fails
        
        # Parse the response as GeoJSON
        earthquakes = geojson.loads(response.text)
        return earthquakes

    except requests.RequestException as e:
        print(f"Error accessing primary API: {e}")
        return None
    except geojson.GeoJSONDecodeError as e:
        print(f"Error parsing GeoJSON data: {e}")
        return None


def get_ffm_geojson_url(event_id):
    """
    Retrieves the URL to the FFM.geojson for a given earthquake event ID.
    """
    detail_url = f"https://earthquake.usgs.gov/fdsnws/event/1/query"
    params = {
        "eventid": event_id,
        "format": "geojson"
    }

    print(f"Fetching event detail for {event_id}...")
    response = requests.get(detail_url, params=params)
    response.raise_for_status()
    data = response.json()

    products = data.get("properties", {}).get("products", {})
    finite_faults = products.get("finite-fault", [])

    if not finite_faults:
        print("No finite-fault product available.")
        return None

    ff_product = finite_faults[0]
    contents = ff_product.get("contents", {})

    for key, info in contents.items():
        if key.endswith("FFM.geojson"):
            return info.get("url")

    print("FFM.geojson not found in finite-fault contents.")
    return None


def check_for_new_data(eq_api):
    """
    Fetch data from the USGS Earthquake Portal and returns it as a GeoJSON object.
    The data returned is updated hourly and includes all earthquakes occuring during that period.
    :param eq_api: USGS API endpoint
    :return: GeoJSON object containing earthquake data
    """
    print('=========================================')
    print("Checking for new earthquake data...")
    print('=========================================')

    try:
        # Fetch data from the USGS Earthquake API
        response = requests.get(eq_api)
        response.raise_for_status()  # Raise error if request fails
        
        # Parse the response as GeoJSON
        earthquakes = geojson.loads(response.text)
    
        # Save to a GeoJSON file
        output_file = "earthquakes.geojson"
        with open(output_file, "w") as f:
            json.dump(earthquakes, f, indent=2)
        return earthquakes

    except requests.RequestException as e:
        print(f"Error accessing primary API: {e}")
        return None
    except geojson.GeoJSONDecodeError as e:
        print(f"Error parsing GeoJSON data: {e}")
        return None


def get_coastline(coastline_api):
    """
    Fetch coastline data from OSGEO/PROJ Github repo and return it as a GeoJSON object.
    The data returned will be in the form of a MultiPolygon. The data are also written to a file: "coastline_buffered.geojson".
    :param coastline_api: API endpoint for the coastline data (https://raw.githubusercontent.com/OSGeo/PROJ/refs/heads/master/docs/plot/data/coastline.geojson)
    :return: GeoJSON object containing coastline data
    """
    try:
        # Fetch data from the specified API
        response = requests.get(coastline_api)
        response.raise_for_status()  # Raise error if request fails
        
        # Parse the response as GeoJSON
        coastline = geojson.loads(response.text)
        
        # Extract the LineString features from the GeoJSON data
        features = []
        for feature in coastline["features"]:
            if feature["geometry"]["type"] == "LineString":
                features.append(LineString(feature["geometry"]["coordinates"]))
            else: 
                print("Coastline data is not in LineString format.") # Ensure each feature is a LineString
        
        # Convert LineStrings to Polygons
        polygons = []
        for line in features:
            if line.is_ring:  # Check if the LineString is closed
                polygons.append(Polygon(line))
            else:
                # Close the LineString and create a Polygon
                closed_line = LineString(list(line.coords) + [line.coords[0]])
                polygons.append(Polygon(closed_line))

        # Combine all polygons into a MultiPolygon
        coastline = MultiPolygon(polygons)

        # Buffer the coastline by 0.5 degrees        
        coastline = coastline.buffer(0.5)

        # Convert the MultiLineString to GeoJSON format
        geojson_data = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": mapping(coastline),
                    "properties": {}
                }
            ]
        }

        # Save to a GeoJSON file
        output_file = "coastline_buffered.geojson"
        with open(output_file, "w") as f:
            json.dump(geojson_data, f, indent=2)
        return coastline
    
    except requests.RequestException as e:
        print(f"Error accessing coastline API: {e}")
        return None
    except geojson.GeoJSONDecodeError as e:
        print(f"Error parsing GeoJSON data: {e}")
        return None
    

def parse_geojson(geojson_data):
    """
    Parse the features of a GeoJSON object and create a dictionary for each earthquake (feature),
    with property names as the keys and property values as the values.
    :param geojson_data: GeoJSON object containing earthquake data
    :return: List of dictionaries containing earthquake data
    """
    earthquakes = []

    # Loop through each feature in the GeoJSON data
    for feature in geojson_data['features']:
        # Extract the properties of the feature
        properties = feature['properties']
        
        # Extract the geometry (coordinates) of the feature
        geometry = feature['geometry']
        coordinates = geometry['coordinates'] if geometry and 'coordinates' in geometry else None
        
        # Create a dictionary for the current feature with property names as keys
        feature_dict = {key: value for key, value in properties.items()}
        
        # Add geometry coordinates to the dictionary
        feature_dict['coordinates'] = coordinates
        
        # Add the USGS ID from the GeoJSON data
        feature_dict['id'] = feature['id']

        # Append the dictionary to the list
        earthquakes.append(feature_dict)
    return earthquakes


def withinCoastline(earthquake, coastline):
    """
    Determine if earthquake epicenter is within 0.5 decimal degrees (~55 km) of the coastline.
    This is one filtering parameter to determine if an earthquake is "significant" within the scope of this project.
    :param earthquake: dictionary containing earthquake data
    :param coastline: shapely Polygon object representing the coastline
    :return: True if the epicenter is within the coastline, False otherwise
    """
    # Extract the coordinates of the earthquake's epicenter
    coords = earthquake.get('coordinates', [])
    if not coords:
        return None

    # Create a Point object from the earthquake's coordinates
    epicenter = Point(coords[:2])
    
    # Buffer the coastline by 0.5 degrees
    coastline_buffer = coastline.buffer(0.5)

    # Determine if the epicenter is within the coastline
    within_coastline_buffer = coastline_buffer.contains(epicenter)
    return within_coastline_buffer


def check_significance(earthquakes, start_date, end_date=None, mode='historic'):
    """
    Check the significance of each earthquake based on its 
    (1) magnitude (>=6.0), (2) USGS alert level (['green','yellow','orange','red]),
    (3) depth (<=30.0 km), and (4) distance from land (within 0.5 degrees, ~55 km of the coastline).
    :param earthquakes: list of dictionaries containing earthquake data
    :return: List of dictionaries containing significant earthquakes
    """
    print('=========================================')
    print("Checking for significant earthquakes...")
    print('=========================================')

    significant_earthquakes = []
    alert_list = ['green','yellow', 'orange', 'red']
    coastline = get_coastline(coastline_api)

    # Include alert criteria for historic data
    if mode == 'historic':
        for earthquake in earthquakes:
            magnitude = earthquake.get('mag')
            alert = earthquake.get('alert')
            depth = earthquake.get('coordinates', [])[2] if earthquake.get('coordinates') else None
            within_Coastline_buffer = withinCoastline(earthquake, coastline)
            if all(var is not None for var in (magnitude, alert, depth)):
                if (magnitude >= 5.5) and (alert in alert_list) and (depth <= 40.0) and within_Coastline_buffer:
                    significant_earthquakes.append(earthquake)

    # Base significance on magnitude, depth, and distance from land for forward-looking data
    if mode =='forward':
        for earthquake in earthquakes:
            magnitude = earthquake.get('mag')
            depth = earthquake.get('coordinates', [])[2] if earthquake.get('coordinates') else None
            within_Coastline_buffer = withinCoastline(earthquake, coastline)
            if all(var is not None for var in (magnitude, depth)):
                if (magnitude >= 5.5) and (depth <= 15.0) and within_Coastline_buffer:
                    significant_earthquakes.append(earthquake)

    # Write significant earthquakes to a GeoJSON file
    if len(significant_earthquakes) > 0:
        print('=========================================')
        print(f"Found {len(significant_earthquakes)} significant earthquakes.")
        print('=========================================')
        for eq in significant_earthquakes:
            print(f"Name: {eq['title']}")
            print(f"Rupture Date/Time: {convert_time(eq['time'])} UTC")
            print(f"Magnitude: {eq['mag']}")
            print(f"Depth: {eq['coordinates'][2]} km")
            print(f"Alert Level: {eq['alert']}")
            print('=========================================')
        print('=========================================')
        significant_earthquakes_to_geojson_and_csv(significant_earthquakes, start_date, end_date)
        return significant_earthquakes
    else:
        return None


def significant_earthquakes_to_geojson_and_csv(significant_earthquakes, start_date, end_date=None):
    """
    Write the significant earthquakes to a GeoJSON file, namely: "significant_earthquakes_full_record.geojson"
    :param significant_earthquakes: list of dictionaries containing significant earthquake metadata
    """
    geojson_features = []
    for eq in significant_earthquakes:
        # Convert time to a datetime object
        dt = convert_time(eq['time'])

        # Create a GeoJSON feature for each earthquake
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Point",
                "coordinates": eq["coordinates"][:2]
            },
            "properties": {
                "place": eq["place"],
                "magnitude": eq["mag"],
                "date": dt.strftime("%Y-%m-%d"),
                "time_utc": dt.strftime("%H:%M:%S"),
                "longitude": eq["coordinates"][0],
                "latitude": eq["coordinates"][1],
                "depth_km": eq["coordinates"][2],
                "alert": eq["alert"],
                "url": eq["url"]
            }
        }
        geojson_features.append(feature)

    # Create the final GeoJSON structure
    geojson_data = {
        "type": "FeatureCollection",
        "features": geojson_features
    }

    # Output the data to GeoJSON and CSV files
    if start_date and end_date:
        with open(f'significant_earthquakes_{start_date}_to_{end_date}_M6.geojson', 'w') as f:
            geojson.dump(geojson_data, f)

        with open(f'significant_earthquakes_{start_date}_to_{end_date}.csv', 'w', newline='') as f:
            writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["Place", "Magnitude", "Date", "Time_utc", "Longitude", "Latitude", "Depth_km", "Alert", "URL"])
            
            for eq in significant_earthquakes:
                dt = convert_time(eq['time'])
                writer.writerow([
                    eq["place"], eq["mag"], dt.strftime("%Y-%m-%d"), dt.strftime("%H:%M:%S"),
                    eq["coordinates"][0], eq["coordinates"][1], eq["coordinates"][2], eq["alert"], eq["url"]
                ])
    else:
        with open(f'significant_earthquakes_{start_date}.geojson', 'w') as f:
            geojson.dump(geojson_data, f)

        with open(f'significant_earthquakes_{start_date}.csv', 'w', newline='') as f:
            writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["Place", "Magnitude", "Date", "Time_utc", "Longitude", "Latitude", "Depth_km", "Alert", "URL"])
            
            for eq in significant_earthquakes:
                dt = convert_time(eq['time'])
                writer.writerow([
                    eq["place"], eq["mag"], dt.strftime("%Y-%m-%d"), dt.strftime("%H:%M:%S"),
                    eq["coordinates"][0], eq["coordinates"][1], eq["coordinates"][2], eq["alert"], eq["url"]
            ])
    return


def make_aoi(coordinates):
    """
    Create an Area of Interest (AOI) polygon based on the given coordinates.
    The AOI is a square with a side length of 1 degree (~111 km) centered on the earthquake's epicenter.
    The AOI is written to a GeoJSON file, "AOI.geojson".
    :param coordinates: list containing the longitude and latitude of the earthquake's epicenter
    :return: Shapely Polygon object representing the AOI
    """
    print('=========================================')
    print("Creating Area of Interest (AOI) polygon...")
    print('=========================================')

    # Extract the X and Y coordinates of the earthquake's epicenter
    X = coordinates[0]
    Y = coordinates[1]

    # Define the side length of the square AOI in decimal degrees
    side_length = 1.0 # 1 degree is ~111 km at the equator

    # Calculate half side length
    half_side = side_length / 2
    
    # Define the square's vertices relative to the center point
    square_coords = [
        (X - half_side, Y - half_side),  # LL
        (X + half_side, Y - half_side),  # LR
        (X + half_side, Y + half_side),  # UR
        (X - half_side, Y + half_side)   # UL
    ]
    
    # Create the square polygon
    AOI = Polygon(square_coords)

    print(f"Area of Interest (AOI) created: {AOI}")
    print('=========================================')
    return AOI


def load_aoi_from_json(aoi_path_or_url):
    """
    Load an AOI from a given GeoJSON file or URL and return a Shapely Polygon or MultiPolygon.

    :param aoi_path_or_url: Path to the AOI GeoJSON file, or URL.
    :return: Shapely Polygon or MultiPolygon representing the AOI.
    """
    try:
        # Determine if input is a URL
        parsed = urlparse(aoi_path_or_url)
        if parsed.scheme in ("http", "https"):
            print(f"Loading AOI from URL: {aoi_path_or_url}")
            response = requests.get(aoi_path_or_url)
            response.raise_for_status()
            aoi_data = response.json()
        else:
            print(f"Loading AOI from file: {aoi_path_or_url}")
            with open(aoi_path_or_url, 'r') as f:
                aoi_data = json.load(f)

        # Ensure the GeoJSON contains features
        features = aoi_data.get("features", [])
        if not features:
            raise ValueError("Invalid AOI file: No features found.")

        # Convert all geometries and union them
        geometries = [shape(f["geometry"]) for f in features if "geometry" in f]
        if not geometries:
            raise ValueError("No valid geometries in features.")

        combined = unary_union(geometries)
        minx, miny, maxx, maxy = combined.bounds
        return box(minx, miny, maxx, maxy)

    except Exception as e:
        print(f"Error loading AOI from JSON: {e}")
        exit(1)


def convert_time(time):
    """
    Convert the given Unix timestamp in milliseconds to a UTC datetime object.
    :param time_ms: Unix timestamp in milliseconds (int or float)
    :return: Datetime object in UTC in this format: 'YYYY-MM-DDTHH:MM:SS'
    """
    timestamp_s = time / 1000 # Convert to milliseconds
    dt = datetime.fromtimestamp(timestamp_s, tz=timezone.utc) # Convert to datetime object in UTC
    dt = dt.replace(microsecond=0) # Remove microseconds
    return dt


def make_interactive_map(frame_dataframe, title, coords, url):
        
    # Extract latitude and longitude from coords
    lon, lat = coords[0], coords[1]
    
    # Generate an interactive map
    map_object = frame_dataframe.explore(column="pathNumber", 
                                            cmap="viridis", 
                                            tooltip=["flightDirection", "pathNumber", "frameNumber", "startTime"],
                                            )

    # Add a basemap explicitly with Folium (Esri World Imagery)
    folium.TileLayer('Esri World Imagery').add_to(map_object)

    # Create a popup with both the title and a clickable URL
    popup_content = f"""
    <b>{title}</b><br>
    <a href="{url}" target="_blank">{url}</a>
    """
    # Add AOI centroid marker with a popup that includes the title and URL
    folium.Marker(
        location=[lat, lon],
        popup=folium.Popup(popup_content, max_width=300),
        icon=folium.Icon(color='red', icon='info-sign')
    ).add_to(map_object)

    # Save to an HTML file
    map_filename = f"{title}_SLC_Map.html"
    map_object.save(map_filename)
    
    return map_filename


def get_path_and_frame_numbers(AOI, time):
    """
    Query the ASF DAAC API for SLC data intersecting the Area of Interest (AOI) over the previous 24 days.
    This ensures all possible intersecting SLC fileIDs are returned for the given AOI.
    :param AOI: Shapely Polygon object representing the Area of Interest
    :param time: Unix timestamp representing the earthquake's origin time
    :return: Dictionary containing the path and frame numbers for each *unique* intersecting SLC.
    """
    # Establish the date range for the query
    rupture_date = convert_time(time)
    start_date = rupture_date - timedelta(days=24) # 24 days before the earthquake
    start_date = start_date.replace(hour=0, minute=0, second=0)
    end_date = rupture_date.replace(hour=23, minute=59, second=59) # the day of the earthquake

    # Format the datetime object into a string
    start_date = start_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    end_date = end_date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # Define the ASF query parameters
    params = {
        'intersectsWith': AOI.wkt,
        'dataset': 'SENTINEL-1',
        'processingLevel': 'SLC',
        'beamSwath': 'IW',
        'start': start_date,
        'end': end_date,
        'output': 'geojson'
    }

    print('Performing ASF DAAC API query to return path and frame numbers for SLCs intersecting AOI over the preceding 24 days...')

    # Sometimes the request to ASF DAAC times out for various reasons. This logic is meant to reduce that.
    MAX_RETRIES = 10
    WAIT_SECONDS = 30

    for attempt in range(MAX_RETRIES):
        try:
            # Fetch data from the ASF DAAC API
            print(f"Attempt {attempt + 1} of {MAX_RETRIES}")
            response = requests.get(ASF_DAAC_API, params=params, timeout=160)
            response.raise_for_status()

            # Parse the response as GeoJSON
            data = geojson.loads(response.text)

            # Check if features are returned
            if not data.get('features'):
                print("No SLC features found intersecting the AOI.")
                return {}, gpd.GeoDataFrame()

            # Convert the GeoJSON data to a GeoDataFrame for visualization
            # Wrapped in try/except to handle cases where features lack geometry
            try:
                frame_dataframe = gpd.GeoDataFrame.from_features(data['features'], crs='EPSG:4326')
            except Exception as e:
                print(f"Warning: Could not create GeoDataFrame (likely missing geometry): {e}")
                return {}, gpd.GeoDataFrame()

            # Initialize an empty dictionary to store the path and frame numbers as sets
            path_frame_numbers = defaultdict(lambda: defaultdict(set))

            # Extract the path and frame numbers from the GeoJSON data
            for feature in data['features']:
                flight_direction = feature['properties']['flightDirection']
                start_time = dateparser.isoparse(
                    feature['properties']['startTime']
                    ).strftime("%Y-%m-%d %H:%M:%S") + ' UTC'
                path_number = feature['properties']['pathNumber']
                frame_number = feature['properties']['frameNumber']
                path_frame_numbers[flight_direction][path_number].add((frame_number,start_time))  # Use a set to avoid duplicates

            # Reformat the dictionary into usable format with tuples as keys and lists as values
            reformatted = {
                (flight_direction, path_number): sorted(frame_numbers)
                for flight_direction, path_frame in path_frame_numbers.items()
                for path_number, frame_numbers in path_frame.items()
            }

            print('=========================================')
            print('Reformatted Path and Frame Numbers for SLCs:')
            print('=========================================')
            for key, value in reformatted.items():
                print(f"{key}: {value}")
            return reformatted, frame_dataframe

        except requests.exceptions.RequestException as e:
            print(f"Request error from ASF DAAC API: {e}")
            if attempt < MAX_RETRIES - 1:
                print(f"Retrying in {WAIT_SECONDS} seconds...")
                sleep(WAIT_SECONDS)
            else:
                print("All retry attempts failed.")
                return {}, gpd.GeoDataFrame()

        except Exception as e:
            print(f"Unexpected error while processing ASF DAAC response: {e}")
            # IMPORTANT: Return empty structures to prevent unpacking errors in main loop
            return {}, gpd.GeoDataFrame()

    # Fallback if loop finishes without returning
    return {}, gpd.GeoDataFrame()


def get_SLCs(flight_direction, path_number, frame_numbers, time, processing_mode):
    """
    Query the ASF DAAC API for SLC data based on the given path and frame numbers.
    The data are organized by flight direction, path number, and frame numbers.
    :param processing_mode: 'historic', 'forward'
    :param flight_direction: 'ASCENDING' or 'DESCENDING'
    :param path_number: Sentinel-1 path number
    :param frame_numbers: List of Sentinel-1 frame numbers
    :param time: Unix timestamp representing the earthquake's origin time
    :return: List of dictionaries containing SLC fileIDs and their respective dates
    """
    # Establish the date range for the query
    rupture_date = convert_time(time)
    
    if processing_mode == 'historic':
        start_date = rupture_date - timedelta(days=90)  # 90 days before the earthquake
        start_date= start_date.replace(hour=0, minute=0, second=0)
        end_date = rupture_date + timedelta(days=30)    # 30 days after the earthquake
        end_date = end_date.replace(hour=23, minute=59, second=59)
    elif processing_mode == 'forward':
        # In forward mode for tracking, we want data acquired AFTER the earthquake
        start_date = rupture_date
        today = datetime.now(timezone.utc)
        end_date = today

    # Format the datetime object into a string
    start_date = start_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    end_date = end_date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # Define the query parameters
    params = {
        'flightDirection': flight_direction,
        'frame': ','.join(str(f) for f in frame_numbers),
        'relativeOrbit': path_number,
        'dataset':'SENTINEL-1',
        'processingLevel': 'SLC',
        'beamSwath': 'IW',
        'start':start_date,
        'end':end_date,
        'output': 'geojson'
    }

    print('Performing ASF DAAC API query to return SLCs for the given path and frame numbers...')
    
    # Sometimes the request to ASF DAAC times out for various reasons. This logic is meant to reduce that.
    MAX_RETRIES = 10
    WAIT_SECONDS = 30

    for attempt in range(MAX_RETRIES):
        try:
            # Fetch data from the ASF DAAC API
            print(f"Attempt {attempt + 1} of {MAX_RETRIES}")
            response = requests.get(ASF_DAAC_API, params=params, timeout=160)
            response.raise_for_status()
            # Parse the response as GeoJSON
            data = geojson.loads(response.text)

            # Extract the file IDs from the GeoJSON data
            SLCs = []
            for feature in data['features']:
                start_time = feature['properties']['startTime']
                path = feature['properties']['pathNumber']
                frame = feature['properties']['frameNumber']

                try:
                    date = isoparse(start_time).strftime("%Y-%m-%dT%H:%M:%SZ")
                except Exception:
                    print(f"Warning: Unexpected date format in startTime: {start_time}")
                    date = None
                
                SLC = {
                    'fileID': feature['properties']['fileID'],
                    'date': date,
                    'pathNumber': path,
                    'frameNumber': frame
                }

                SLCs.append(SLC)

            # Print the SLCs
            print('=========================================')
            print(f"Found {len(SLCs)} SLCs for the {flight_direction} path {path_number} and frame numbers {frame_numbers}.")
            print('=========================================')
            for SLC in SLCs:
                print(f"FileID: {SLC['fileID']}, Date: {SLC['date']}")
            return SLCs
        
        except requests.exceptions.RequestException as e:
            print(f"Request error from ASF DAAC API: {e}")
            if attempt < MAX_RETRIES - 1:
                print(f"Retrying in {WAIT_SECONDS} seconds...")
                sleep(WAIT_SECONDS)
            else:
                print("All retry attempts failed.")
                return None


def search_copernicus_public(aoi_polygon, start_date, end_date):
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
    params = {
        "$filter": filter_query,
        "$orderby": "ContentDate/Start asc",
        "$top": 1000
    }
    
    print(f"Searching Copernicus (Public OData) for Sentinel-2 L1C...")
    
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
            products = data.get('value', [])
            
            for prod in products:
                name = prod.get('Name')
                start = prod.get('ContentDate', {}).get('Start')
                
                # --- PARSE GEOMETRY ---
                # Format is: geography'SRID=4326;POLYGON ((...))'
                raw_footprint = prod.get('Footprint')
                footprint = None
                
                if raw_footprint:
                    try:
                        # Extract just the WKT part (POLYGON...)
                        # Split by semicolon to remove SRID
                        # Remove trailing quote if present
                        clean_wkt = raw_footprint.split(';')[-1].replace("'", "")
                        footprint = wkt.loads(clean_wkt)
                    except Exception as e:
                        # print(f"Geometry parse error: {e}") # Optional debug
                        pass 

                # --- EXTRACT CLOUD COVER ---
                cloud_cover = 0.0
                attrs = prod.get('Attributes', [])
                for attr in attrs:
                    if attr.get('Name') == 'cloudCover':
                        cloud_cover = attr.get('Value', 0.0)
                        break

                # --- DETERMINE PLATFORM ---
                platform = "sentinel-2"
                if name.startswith("S2A"): platform = "sentinel-2a"
                elif name.startswith("S2B"): platform = "sentinel-2b"
                elif name.startswith("S2C"): platform = "sentinel-2c"

                granules.append({
                    "granule_id": name,
                    "date": start.split("T")[0],
                    "datetime": start,
                    "cloud_cover": cloud_cover,
                    "platform": platform,
                    "footprint": footprint # Successfully parsed geometry
                })
            
            next_link = data.get('@odata.nextLink')
            
        except Exception as e:
            print(f"Error searching Public OData: {e}")
            break

    print(f"  Found {len(granules)} Sentinel-2 products.")
    return granules


def search_element84_stac(aoi_polygon, start_date, end_date):
    """
    Searches Element84 Earth Search STAC API v1 for Sentinel-2 L2A COGs using pystac_client.
    """
    print("Searching Element84 Earth Search for Sentinel-2 L2A...")
    
    # Connect to the API endpoint
    api_url = "https://earth-search.aws.element84.com/v1"
    client = Client.open(api_url)
    
    # Define parameters
    collection_sentinel_2_l2a = "sentinel-2-c1-l2a" 
    intersects_geom = mapping(aoi_polygon)
    
    try:
        # Perform the search using the client
        search = client.search(
            collections=[collection_sentinel_2_l2a],
            intersects=intersects_geom,
            datetime=f"{start_date}/{end_date}",
            query={"eo:cloud_cover": {"lte": OPTICAL_CLOUD_THRESHOLD}},
            limit=100
        )
        
        # Retrieve the metadata items
        items = list(search.items())
        print(f"  Found {len(items)} Sentinel-2 products.")
        
        granules = []
        
        # Each item contains information about the scene geometry, acquisition time, and properties
        for item in items:
            props = item.properties
            
            platform = props.get('platform', 'sentinel-2')
            cloud_cover = props.get('eo:cloud_cover', 0.0)
            
            # Extract Tile ID from grid:code (e.g., 'MGRS-35SNA' -> '35SNA')
            grid_code = props.get('grid:code') or props.get('s2:mgrs_tile', '')
            tile_id = grid_code.replace('MGRS-', '') if grid_code else ''
            
            # Extract UTM Zone directly from mgrs:utm_zone
            utm_zone = str(props.get('mgrs:utm_zone', 'Unknown'))
            
            # Extract Orbit ID from s2:product_uri (e.g., '..._R107_...')
            product_uri = props.get('s2:product_uri', '')
            orbit_match = re.search(r'_R(\d{3})_', product_uri)
            orbit_id = orbit_match.group(1) if orbit_match else "000"
            
            granules.append({
                "granule_id": item.id,
                "date": props.get('datetime', '').split('T')[0],
                "datetime": props.get('datetime'),
                "cloud_cover": cloud_cover,
                "platform": platform,
                "footprint": shape(item.geometry),
                "orbit_id": orbit_id,
                "tile_id": tile_id,
                "utm_zone": utm_zone
            })
            
        return granules

    except Exception as e:
        print(f"  !! Error searching Element84 STAC with pystac_client: {e}")
        return []


def get_utm_zone(granule_id):
    """Extracts UTM zone from Sentinel-2 Granule ID (e.g., 'T46QHK' -> '46').
    :param granule_id: The granule ID string from which to extract the UTM zone.
    :return: The UTM zone as a string, or 'Unknown' if it cannot be extracted.
    """
    match = re.search(r'_T(\d{2})[A-Z]{3}_', granule_id)
    return match.group(1) if match else "Unknown"


def make_optical_job_json(title, orbit_id, pre_date, post_date, reference_ids, secondary_ids, pre_cc=None, post_cc=None, status="COMPLETE", zone_suffix=None):
    """
    Helper function to create a JSON object for an AUTORIFT job. Matches the schema defined in ARIA_AUTORIFT.yml.
    :param title: Title of the job (usually the earthquake event name)
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
        "pre_cloud_cover": round(pre_cc, 2) if pre_cc is not None else None,    # Temporary debug field
        "post_cloud_cover": round(post_cc, 2) if post_cc is not None else None,  # Temporary debug field
        "job_parameters": {
            "reference": reference_ids,
            "secondary": secondary_ids,
            "chip_size": 24,
            "search_range": 64
        }
    }
    
    # Internal tracking for partial jobs (optional)
    if status != "COMPLETE":
        job_json["status_note"] = status
        
    return job_json


def find_optical_pairs_copernicus(optical_scenes, rupture_time, title, aoi_polygon, job_list=True):
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
        granule_id = scene['granule_id']
        match = re.search(r'_R(\d{3})_', granule_id)
        if match:
            orbit_id = match.group(1)
            orbit_groups[orbit_id][scene['date']].append(scene)

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
            platform = scenes[0]['platform'] if scenes else "sentinel-2"

            unique_scenes = []
            seen_tiles = set()
            scenes.sort(key=lambda x: x['cloud_cover'])
            
            combined_poly = None
            
            for s in scenes:
                t_match = re.search(r'_T(\w{5})_', s['granule_id'])
                if t_match:
                    tile_id = t_match.group(1)
                    if tile_id not in seen_tiles:
                        unique_scenes.append(s)
                        seen_tiles.add(tile_id)
                        
                        if s.get('footprint'):
                            if combined_poly is None:
                                combined_poly = s['footprint']
                            else:
                                combined_poly = combined_poly.union(s['footprint'])

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
            avg_cc = sum(s['cloud_cover'] for s in unique_scenes) / len(unique_scenes) if unique_scenes else 100
            delta_days = abs((date_dt - rupture_dt).days)
            
            candidate = {
                'date': date_str,
                'scenes': unique_scenes,
                'platform': platform,
                'coverage': coverage_pct, 
                'tile_count': tile_count,
                'cc': avg_cc,             
                'delta': delta_days       
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
                best_pre = min(pre_candidates, key=lambda x: x['cc'])
                print(f"    - Found Pre-event: {best_pre['date']} (CC: {best_pre['cc']:.1f}%)")
                pre_date_str = best_pre['date']
                primary_ids = [s['granule_id'].replace('.SAFE', '') for s in best_pre['scenes']]
            else:
                print(f"    - Missing PRE-event coverage.")

            if post_candidates:
                best_post = min(post_candidates, key=lambda x: x['cc'])
                print(f"    - Found Post-event: {best_post['date']} (CC: {best_post['cc']:.1f}%)")
                post_date_str = best_post['date']
                secondary_ids = [s['granule_id'].replace('.SAFE', '') for s in best_post['scenes']]
            else:
                print(f"    - Missing POST-event coverage.")
            
            # Generate the partial job and add to list
            if job_list:
                job = make_optical_job_json(
                    title, orbit_id, pre_date_str, post_date_str, 
                    primary_ids, secondary_ids, status="PARTIAL"
                )
                jobs.append(job)
            continue

        # --- STANDARD SELECTION ---
        # Sort: Max Coverage > Low Cloud > Low Delta
        pre_candidates.sort(key=lambda x: (-x['coverage'], x['cc'], x['delta']))
        best_pre = pre_candidates[0]
        
        post_candidates.sort(key=lambda x: (-x['coverage'], x['cc'], x['delta']))
        best_post = post_candidates[0]
        
        print(f"  Orbit {orbit_id}: Pre={best_pre['date']} ({best_pre['coverage']:.1f}% Cov), Post={best_post['date']} ({best_post['coverage']:.1f}% Cov)")

        for stage, candidate in [("Pre-Event", best_pre), ("Post-Event", best_post)]:
            for scene in candidate['scenes']:
                if scene.get('footprint'):
                    feature = {
                        "type": "Feature",
                        "geometry": mapping(scene['footprint']),
                        "properties": {
                            "earthquake": title,
                            "role": stage,
                            "date": candidate['date'],
                            "orbit": orbit_id,
                            "platform": candidate['platform'],
                            "granule_id": scene['granule_id'],
                            "cloud_cover": scene['cloud_cover'],
                            "coverage_pct": round(candidate['coverage'], 1)
                        }
                    }
                    scene_features.append(feature)

        if job_list:
            primary_ids = [s['granule_id'].replace('.SAFE', '') for s in best_pre['scenes']]
            secondary_ids = [s['granule_id'].replace('.SAFE', '') for s in best_post['scenes']]
            
            # CALL THE HELPER FUNCTION HERE
            job = make_optical_job_json(
                title, orbit_id, best_pre['date'], best_post['date'], 
                primary_ids, secondary_ids
            )
            jobs.append(job)
            
    return jobs, scene_features


def process_candidate_group(orbit_key, dates_dict, rupture_dt, aoi_polygon, title, aoi_area, strategy_name, role=None):
    """
    Generic logic to select best Pre/Post pair from a grouped dictionary of candidates. Used by both 'Dominant' and 'Split' strategies.
    :param orbit_key: The key representing the group (e.g., "047" for Orbit-based, "047_Z46" for Orbit_Zone-based)
    :param dates_dict: Dictionary where keys are date strings and values are lists of scene dictionaries for that date
    :param rupture_dt: Datetime object representing the earthquake's origin time
    :param aoi_polygon: Shapely Polygon representing the Area of Interest (used for coverage calculation)
    :param title: Title of the earthquake event (used for job naming)
    :param aoi_area: Area of the AOI polygon (used for coverage calculation)
    :param strategy_name: Name of the strategy ("Dominant" or "Split") for logging purposes
    :param role: Optional parameter to indicate if this group is 'dominant' or 'minority' in the context of mixed zones (used for logging and feature properties)
    :return: A tuple containing a list of job JSON objects and a list of GeoJSON features for visualization.
    """
    pre_candidates = []
    post_candidates = []
    
    # Evaluate candidates for every available date
    for date_str, scenes in dates_dict.items():
        if not scenes: continue
        
        date_dt = datetime.strptime(date_str, "%Y-%m-%d").replace(tzinfo=timezone.utc)
        platform = scenes[0]['platform']

        # Deduplicate tiles for this date/group
        unique_scenes = []
        seen_tiles = set()
        scenes.sort(key=lambda x: x['cloud_cover'])
        
        combined_poly = None
        for s in scenes:
            # Create a unique key for the tile (UTM + Triplet, e.g., '46QHK')
            tile_id = s.get('tile_id')
            # Check if tile_id exists and hasn't been processed yet
            if tile_id and tile_id not in seen_tiles:
                unique_scenes.append(s)
                seen_tiles.add(tile_id)
                if s.get('footprint'):
                    if combined_poly is None: 
                        combined_poly = s['footprint']
                    else: 
                        combined_poly = combined_poly.union(s['footprint'])

        if not unique_scenes:
            continue

        # Calc coverage
        coverage_pct = 0.0
        if combined_poly and aoi_polygon:
            try:
                clean_poly = combined_poly.buffer(0)
                intersection = clean_poly.intersection(aoi_polygon)
                coverage_pct = (intersection.area / aoi_area) * 100.0
            except: pass
        
        avg_cc = sum(s['cloud_cover'] for s in unique_scenes) / len(unique_scenes)
        delta = abs((date_dt - rupture_dt).days)
        
        candidate = {
            'date': date_str, 'scenes': unique_scenes, 'platform': platform,
            'coverage': coverage_pct, 'tile_count': len(unique_scenes), 'cc': avg_cc, 'delta': delta
        }
        
        if date_dt < rupture_dt: pre_candidates.append(candidate)
        elif date_dt > rupture_dt: post_candidates.append(candidate)

    jobs = []
    features = []

    # Select Best Pair
    if not pre_candidates or not post_candidates:
        return [], []

    # Sort: Max Coverage > Low Cloud > Low Delta
    pre_candidates.sort(key=lambda x: (-x['coverage'], x['cc'], x['delta']))
    best_pre = pre_candidates[0]
    
    post_candidates.sort(key=lambda x: (-x['coverage'], x['cc'], x['delta']))
    best_post = post_candidates[0]
    
    print(f"  [{strategy_name}] {orbit_key}: Pre={best_pre['date']} ({best_pre['coverage']:.1f}% Cov), Post={best_post['date']} ({best_post['coverage']:.1f}% Cov)")

    # Generate Job
    primary_ids = [s['granule_id'].replace('.SAFE', '') for s in best_pre['scenes']]
    secondary_ids = [s['granule_id'].replace('.SAFE', '') for s in best_post['scenes']]
    
    # Handle naming differences based on strategy
    zone_suffix = None
    orbit_id = orbit_key
    if "Z" in orbit_key:
        parts = orbit_key.split('_')
        orbit_id = parts[0]
        zone_suffix = parts[1]

    job = make_optical_job_json(title, orbit_id, best_pre['date'], best_post['date'], primary_ids, secondary_ids, pre_cc=best_pre['cc'], post_cc=best_post['cc'], zone_suffix=zone_suffix)
    jobs.append(job)

    # Generate GeoJSON Features for Visualization
    for stage, candidate in [("Pre-Event", best_pre), ("Post-Event", best_post)]:
        for scene in candidate['scenes']:
            if scene.get('footprint'):
                feat = {
                    "type": "Feature",
                    "geometry": mapping(scene['footprint']),
                    "properties": {
                        "job_name": job['name'],
                        "strategy": strategy_name,
                        "role": role,
                        "timing": stage,
                        "date": candidate['date'],
                        "orbit": orbit_id,
                        "zone": get_utm_zone(scene['granule_id']),
                        "granule_id": scene['granule_id']
                    }
                }
                features.append(feat)
                
    return jobs, features


def check_mixed_zones_in_group(dates_dict):
    """
    Checks if any single date/acquisition within an orbit group contains  tiles from multiple UTM zones. 
    This indicates a 'Mixed' case that requires Split/Dominant strategy comparison.
    :param dates_dict: Dictionary where keys are date strings and values are lists of scene dictionaries for that date
    :return: True if mixed zones are found within any date, False otherwise
    """
    for date, scenes in dates_dict.items():
        zones = set(get_utm_zone(s['granule_id']) for s in scenes)
        if len(zones) > 1:
            return True
    return False


def find_optical_pairs_element84(optical_scenes, rupture_time, title, aoi_polygon):
    """
    Generates optical pairs using two strategies concurrently for comparison:
    1. DOMINANT: Enforces one UTM zone per Orbit (drops minority tiles).
    2. SPLIT: Creates separate jobs for every UTM zone found.

    Uses the Element84 file nomenclature.
    
    Returns:
    - standard_jobs (Dominant jobs for all orbits, for general processing)
    - mixed_dom_features (Features for Dominant strategy, ONLY for mixed cases)
    - mixed_split_features (Features for Split strategy, ONLY for mixed cases)
    """
    rupture_dt = convert_time(rupture_time)
    aoi_area = aoi_polygon.area

    # Group by Orbit -> Date -> List of Scenes
    orbit_groups = defaultdict(lambda: defaultdict(list))
    for scene in optical_scenes:
        if scene.get('orbit_id'):
            orbit_groups[scene['orbit_id']][scene['date']].append(scene)

    # Standard production output (Dominant strategy applied to ALL)
    all_prod_jobs = []
    
    # Comparison/Debug output (ONLY for orbits that actually have mixed zones)
    mixed_dom_feats = []
    mixed_split_feats = []

    print(f"Processing {len(orbit_groups)} tracks with DOMINANT strategy (and SPLIT where applicable)...")

    for orbit_id, dates_dict in orbit_groups.items():
        
        # Check if this orbit even has a mixing issue
        is_mixed = check_mixed_zones_in_group(dates_dict)
        
        # --- STRATEGY A: DOMINANT ZONE (Applied to ALL orbits) ---
        # Note: Strategy A filters to keep only the dominant data per date.
        dom_dates_dict = defaultdict(list)
        for d, scenes in dates_dict.items():
            zone_counts = defaultdict(int)
            for s in scenes: zone_counts[get_utm_zone(s['granule_id'])] += 1
            if zone_counts:
                winner = max(zone_counts, key=zone_counts.get)
                dom_dates_dict[d] = [s for s in scenes if get_utm_zone(s['granule_id']) == winner]
        
        # Generate Dominant Jobs (Always added to production list)
        j_dom, f_dom = process_candidate_group(orbit_id, dom_dates_dict, rupture_dt, aoi_polygon, title, aoi_area, "DOMINANT", role="dominant")
        all_prod_jobs.extend(j_dom)

        # If this was a Mixed case, add to our specific debug lists
        if is_mixed:
            mixed_dom_feats.extend(f_dom)

            # --- STRATEGY B: SPLIT ZONES (Only needed for Mixed cases) ---
            # Determine the Global Dominant Zone (for the whole orbit, not just per date)
            # This allows us to label the split results as 'dominant' or 'minority'
            all_orbit_scenes = [s for scenes in dates_dict.values() for s in scenes]
            global_zone_counts = defaultdict(int)
            for s in all_orbit_scenes:
                global_zone_counts[get_utm_zone(s['granule_id'])] += 1
            
            global_winner = None
            if global_zone_counts:
                global_winner = max(global_zone_counts, key=global_zone_counts.get)

            split_groups = defaultdict(lambda: defaultdict(list))
            for d, scenes in dates_dict.items():
                for s in scenes:
                    z = get_utm_zone(s['granule_id'])
                    key = f"{orbit_id}_Z{z}"
                    split_groups[key][d].append(s)
            
            for key, z_dates_dict in split_groups.items():
                # Extract zone from key to determine role
                z_str = key.split('_Z')[-1]
                role = 'dominant' if z_str == global_winner else 'minority'

                _, f_split = process_candidate_group(key, z_dates_dict, rupture_dt, aoi_polygon, title, aoi_area, "SPLIT", role=role)
                mixed_split_feats.extend(f_split)

    return all_prod_jobs, mixed_dom_feats, mixed_split_feats


def generate_pairs(pairs, mode):
    """
    Generate pairs of SLCs based on the selected pairing mode.
    :param pairs: List of SLC pairs sorted by date
    :param: mode: 'sequential' for temporally consecutive pairs, 'all' for all possible pairs, 'conseismic' for pairs bounding the rupture date only
    :return: List of SLC pairs based on the mode
    """
    if mode == 'sequential':
        return [(pairs[i], pairs[i + 1]) for i in range(len(pairs) - 1)]
    elif mode == 'all':
        all_pairs = list(combinations(pairs, 2))
        return all_pairs
    elif mode == 'coseismic':
        return []


def find_reference_and_secondary_pairs(SLCs, time, flight_direction, path_number, title, pairing_mode='sequential', job_list = False, resolution=90):
    """
    Find the reference and secondary pairs of SLCs necessary to run dockerized topsApp, 
    and determine whether each pair is pre-seismic, co-seismic, or post-seismic based on the rupture date and SLC dates.
    :param SLCs: List of dictionaries containing SLC fileIDs and their respective dates
    :param time: Unix timestamp representing the earthquake's origin time
    :param flight_direction: 'ASCENDING' or 'DESCENDING'
    :param path_number: Sentinel-1 path number
    :param title: USGS title of the earthquake event, used for file organization
    :param pairing_mode: 'sequential' for temporally consecutive pairs, 'all' for all possible pairs, 'coseismic' for pairs bounding the rupture date only
    :param job_list: True if the JSON objects are for HYP3 job submission, False otherwise
    :param resolution: Output resolution for the topsApp processing, default is 90m
    :return: List of JSON objects containing the parameters for each pair of SLCs
    """
    # Get the rupture date in the format YYYY-MM-DD
    rupture_date = convert_time(time)
    rupture_date = rupture_date.strftime('%Y-%m-%d')
    rupture_date_dt = convert_time(time).replace(tzinfo=None)
    
    # Reformatting for dictionary keys for later use    
    flight_direction = 'A' if flight_direction == 'ASCENDING' else 'D'
    path_number = f"{int(path_number):03}" 

    # Pair SLCs by date 
    slc_by_date = {}
    for slc in SLCs:
        date_obj = datetime.strptime(slc['date'][:10], "%Y-%m-%d")
        key = date_obj
        if key not in slc_by_date:
            slc_by_date[key] = []
        slc_by_date[key].append(slc)
    
    sorted_dates = sorted(slc_by_date.keys())
    initial_pairs = []
    for date in sorted_dates:
        frames = slc_by_date[date]
        frame_numbers = {slc['frameNumber'] for slc in frames}
        if len(frames) == len(frame_numbers):
            initial_pairs.append((date, frames))
    
    # Split the pairs into pre-seismic, co-seismic, and post-seismic
    pre_seismic = [pair for pair in initial_pairs if pair[0] < rupture_date_dt]
    post_seismic = [pair for pair in initial_pairs if pair[0] > rupture_date_dt]
    co_seismic = []
    
    for i in range(len(initial_pairs) - 1):
        if initial_pairs[i][0] < rupture_date_dt < initial_pairs[i + 1][0]:
            co_seismic = [(initial_pairs[i], initial_pairs[i + 1])]
            break
    
    # Pair based on the selected pairing_mode
    if pairing_mode == 'coseismic':
        paired_results = {'co-seismic': co_seismic}
    else:
        paired_results = {
            'pre_seismic': generate_pairs(pre_seismic, pairing_mode),
            'co_seismic': co_seismic,
            'post_seismic': generate_pairs(post_seismic, pairing_mode)
        }
    
    # Create JSON objects for each pair
    isce_jsons = []
    for timing, pairs in paired_results.items():
        for (secondary, reference) in pairs:
            reference_date, reference_scenes = reference
            secondary_date, secondary_scenes = secondary
            reference_scenes_ids = [slc['fileID'].removesuffix("-SLC") for slc in reference_scenes]
            secondary_scenes_ids = [slc['fileID'].removesuffix("-SLC") for slc in secondary_scenes]
            if job_list:
                json_output = make_job_json(title, flight_direction, path_number, reference_scenes_ids, secondary_scenes_ids, resolution)
            else:
                json_output = make_json(title, timing, flight_direction, path_number, list(frame_numbers), 
                                        {'date': reference_date.strftime('%Y-%m-%d')}, 
                                        {'date': secondary_date.strftime('%Y-%m-%d')}, 
                                        reference_scenes_ids, secondary_scenes_ids)
            isce_jsons.append(json_output)
    return isce_jsons


def make_json(title, timing, flight_direction, path_number, frame_numbers, reference, secondary, reference_scenes, secondary_scenes):
    """Create a JSON object containing parameters for dockerized topsApp.
    Note: Not all params here are used in the final dockerized topsApp. Some are used for file organzation.
    Note: Several params are 'hardcoded', as these should not vary between individual products.
    :param title: USGS title of the earthquake event
    :param timing: 'pre-seismic', 'co-seismic', or 'post-seismic'
    :param flight_direction: 'A' or 'D' for ascending or descending
    :param path_number: Sentinel-1 path number
    :param frame_numbers: List of intersecting Sentinel-1 frame numbers
    :param reference: Dictionary containing the reference date
    :param secondary: Dictionary containing the secondary date
    :param reference_scenes: List of reference SLC fileIDs
    :param secondary_scenes: List of secondary SLC fileIDs
    :return: JSON object containing the parameters for dockerized topsApp
    """
    # Reformatting 'fight-direction' for readability in the json
    flight_direction = 'ASCENDING' if flight_direction == 'A' else 'DESCENDING'
    
    isce_json = {
        "title": title,
        "timing": timing,
        "flight-direction": flight_direction,
        "path-number": path_number,
        "frame-numbers": frame_numbers,
        "reference-date": reference['date'],
        "secondary-date": secondary['date'],
        "reference-scenes": reference_scenes,
        "secondary-scenes": secondary_scenes,
        "frame-id": -1,
        "estimate-ionosphere-delay": True,
        "esd-coherence-threshold": -1,
        "compute-solid-earth-tide": True,
        "goldstein-filter-power": 0.5,
        "output-resolution": 30,
        "unfiltered-coherence": True,
        "dense-offsets": True,
    }
    return isce_json


def make_job_json(title, flight_direction, path_number, reference_scenes, secondary_scenes, resolution):
    """
    Create a JSON object containing parameters for dockerized topsApp on HYP3.
    :param title: USGS title of the earthquake event
    :param flight_direction: 'A' or 'D' for ascending or descending
    :param path_number: Sentinel-1 path number
    :param reference_scenes: List of reference SLC fileIDs
    :param secondary_scenes: List of secondary SLC fileIDs
    :param resolution: Output resolution for the topsApp processing, default is 90m
    :return: JSON object containing the parameters for dockerized topsApp
    """
    
    job_json = {
        "name": f"{title}-{flight_direction}{path_number}",
        "job_type": "ARIA_S1_COSEIS",
        "job_parameters": {
            "granules": reference_scenes,
            "secondary_granules": secondary_scenes,
            "frame_id": -1,
            "estimate_ionosphere_delay": True,
            "esd_coherence_threshold": -1,
            "compute_solid_earth_tide": True,
            "goldstein_filter_power": 0.5,
            "unfiltered_coherence": True,
            "dense_offsets": True,
            "output_resolution": resolution
            }
    }
    return job_json


def create_directories_from_json(eq_jsons, root_dir):
    """
    Create directories for each group of SLCs based on the JSON data provided. These directories will be used to store the outputs of the dockerized topsApp.
    The directories are created in the root directory specified and will have a name following this format: 'flight_directionpath_number_secondary_date_reference_date'.
    For example 'A012_20220101_20220112' for an ascending path 12 earthquake with secondary date 2022-01-01 and reference date 2022-01-12.
    :param eq_jsons: List of JSON objects containing the parameters for each pair of SLCs
    :param root_dir: Root directory where the directories will be created
    :return: List of directory names
    """
    dirnames = []
    total = 0
    for isce_jsons in eq_jsons:
        sub_dirnames = []  # To hold the directories for each group in `isce_jsons`
        for json_data in isce_jsons:
            title = json_data['title']
            timing = json_data['timing']
            flight_direction = 'A' if json_data['flight-direction'] == 'ASCENDING' else 'D' # Reformat 'fight-direction' to shorten dirname  
            path_number = json_data['path-number']
            secondary_date = json_data['secondary-date'].replace("-", "")
            reference_date = json_data['reference-date'].replace("-", "")
            
            # Build the full path
            base_path = os.path.join(root_dir, title, flight_direction + path_number, timing)
            sub_path = f"{flight_direction}{path_number}_{secondary_date}_{reference_date}"
            full_path = os.path.join(base_path, sub_path)
            
            # Create directories, ensuring no overwriting
            os.makedirs(full_path, exist_ok=True)
            sub_dirnames.append(full_path)
            print(f"Created: {full_path}")
            total += 1
        dirnames.append(sub_dirnames)
    return dirnames, total


def run_dockerized_topsApp(json_data, working_dir):
    """
    Run dockerized topsApp InSAR processing workflow using the provided JSON data.
    Outputs are added to the root dir + an extension for each pair.
    :param json_data: JSON object containing the parameters for dockerized topsApp
    :param working_dir: Working directory where the outputs will be stored
    """
    # Extract the parameters from the JSON data
    params = json_data['job_parameters']
    reference_scenes = " ".join(params['granules'])
    secondary_scenes = " ".join(params['secondary_granules'])
    
    # Define static/dynamic parameters
    frame_id = params.get('frame_id', -1)
    estimate_ionosphere_delay = str(params.get('estimate_ionosphere_delay'))
    esd_coherence_threshold = str(params.get('esd_coherence_threshold'))
    compute_solid_earth_tide = str(params.get('compute_solid_earth_tide'))
    goldstein_filter_power = str(params.get('goldstein_filter_power'))
    output_resolution = str(params.get('output_resolution'))
    unfiltered_coherence = str(params.get('unfiltered_coherence', True))
    dense_offsets = str(params.get('dense_offsets', True))
    
    # Construct the command
    cmd = [
        "conda", "run", "-n", "topsapp_env_trappist_python11",
        "taskset", "-c", "0-16",
        "isce2_topsapp",
        "--reference-scenes", reference_scenes,
        "--secondary-scenes", secondary_scenes,
        "--frame-id", str(frame_id),
        "--estimate-ionosphere-delay", estimate_ionosphere_delay,
        "--esd-coherence-threshold", esd_coherence_threshold,
        "--compute-solid-earth-tide", compute_solid_earth_tide,
        "--goldstein-filter-power", goldstein_filter_power,
        "--output-resolution", output_resolution,
        "--unfiltered-coherence", unfiltered_coherence,
        "--dense-offsets", dense_offsets,
        "++process", "coseis_sar"
    ]
    
    # Set Environment Variables
    env = os.environ.copy()
    env["XLA_PYTHON_CLIENT_MEM_FRACTION"] = ".20"

    print('=========================================')
    print(f"Running dockerized topsApp in {working_dir}...")
    print(f"Command: {' '.join(cmd)}")
    print('=========================================')
    
    # Run the command
    try:
        # Check if working dir exists, create if not
        os.makedirs(working_dir, exist_ok=True)
        
        # Execute
        result = subprocess.run(
            cmd,
            cwd=working_dir,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        print("Processing completed successfully.")
        
        # Optional: Save a log file in that directory
        with open(os.path.join(working_dir, "topsapp_stdout.log"), "w") as log:
            log.write(result.stdout)
            
    except subprocess.CalledProcessError as e:
        print("Error occurred while running topsApp.")
        print("Stderr:\n", e.stderr)
        # Optional: Write error log
        with open(os.path.join(working_dir, "topsapp_error.log"), "w") as log:
            log.write(e.stderr)
        raise e  # Re-raise to handle it in the calling function if needed


def to_snake_case(input_string):
    """
    Convert the given string to snake_case with ASCII-safe characters.
    Strips accents and transliterates Unicode to closest ASCII.
    """
    # Normalize and transliterate Unicode characters to closest ASCII equivalent
    normalized = unicodedata.normalize('NFKD', input_string).encode('ascii', 'ignore').decode('ascii')

    # Replace non-alphanumeric characters with spaces
    cleaned_string = re.sub(r'[^\w\s]', '', normalized)

    # Replace spaces with underscores and convert to lowercase
    snake_case_string = re.sub(r'\s+', '_', cleaned_string.strip()).lower()

    return snake_case_string


def ascii_table_to_html(ascii_table):
    """Converts a raw ASCII table string into a styled HTML table."""
    if not ascii_table or "+" not in ascii_table:
        return f"<pre style='background:#2b2b2b; color:#ccc; padding:10px;'>{ascii_table}</pre>"
    
    lines = ascii_table.strip().split('\n')
    html = '<table style="width: 100%; border-collapse: collapse; font-size: 13px; text-align: left; margin-bottom: 20px;">'
    
    is_header = True
    for line in lines:
        if line.strip().startswith('+'):
            continue  # Skip the border lines
        if line.strip().startswith('|'):
            # Extract cell content, ignoring the first and last empty splits from the outer '|'
            cells = [cell.strip() for cell in line.split('|')[1:-1]]
            html += "<tr>"
            for cell in cells:
                if is_header:
                    html += f'<th style="border-bottom: 2px solid #0055a4; padding: 10px 8px; background-color: #f0f4f8; color: #003366;">{cell}</th>'
                else:
                    html += f'<td style="border-bottom: 1px solid #e0e0e0; padding: 8px; color: #444;">{cell}</td>'
            html += "</tr>"
            is_header = False
            
    html += "</table>"
    return html


def send_email(subject, body, recipients=None):
    """
    Send an email with the earthquake information to a specified list of recipients.
    :param subject: Email subject
    :param body: Email body content
    :param recipients: List of email addresses to send the email to
    """
    if recipients is None:
        recipients = PRIMARY_RECIPIENTS

    if not recipients:
        print("WARNING: Cannot send email. No recipients configured in environment variables.")
        return

    GMAIL_USER = os.environ.get('GMAIL_USER')
    GMAIL_PSWD = os.environ.get('GMAIL_APP_PSWD')

    if not GMAIL_USER or not GMAIL_PSWD:
        print("WARNING: Cannot send email. GMAIL_USER or GMAIL_APP_PSWD environment variables are missing.")
        return

    yag = yagmail.SMTP(GMAIL_USER, GMAIL_PSWD)

    yag.send(
             bcc=recipients,
             subject=subject,
             contents=[body]
             )
    return


def export_gee_composite(aoi_polygon, start_date, end_date, title, stage, gcs_bucket):
    """
    Generates a cloud-free median composite in GEE and exports to Google Cloud Storage.
    :param aoi_polygon: Shapely Polygon representing the Area of Interest
    :param start_date: Start date for the image collection filter (YYYY-MM-DD)
    :param end_date: End date for the image collection filter (YYYY-MM-DD)
    :param title: Title of the earthquake event (used for file naming)
    :param stage: 'pre-event' or 'post-event' to indicate the timing of the composite
    :param gcs_bucket: Name of the Google Cloud Storage bucket to export the composite
    :return: The GEE export task object
    """
    # Convert Shapely polygon to GEE Geometry
    bounds = aoi_polygon.bounds # (minx, miny, maxx, maxy)
    ee_roi = ee.Geometry.Rectangle([bounds[0], bounds[1], bounds[2], bounds[3]])

    # Query the S2 Harmonized Collection
    s2_col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED') \
        .filterBounds(ee_roi) \
        .filterDate(start_date, end_date)

    # Apply Cloud Score+ Masking
    cs_plus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')
    s2_linked = s2_col.linkCollection(cs_plus, ['cs_cdf'])

    def mask_clouds(img):
        # cs_cdf is the clear-sky probability (0 to 1). We keep pixels > 0.65
        mask = img.select('cs_cdf').gte(0.65)
        return img.updateMask(mask)

    s2_masked = s2_linked.map(mask_clouds)

    # Generate the composite using Band 8 (NIR)
    composite = s2_masked.select('B8').median().clip(ee_roi)

    # Create Export Task
    file_name = f"{title}_S2_B8_{stage}_{start_date}_to_{end_date}"
    
    task = ee.batch.Export.image.toCloudStorage(
        image=composite,
        description=file_name,
        bucket=gcs_bucket,
        fileNamePrefix=f"COSEIS_Composites/{title}/{file_name}",
        region=ee_roi,
        scale=10, 
        crs='EPSG:4326',
        maxPixels=1e13
    )
    
    task.start()
    print(f"  Started GCS Export Task: {file_name}")
    return task


def wait_for_gee_tasks(tasks):
    """
    Polls GEE until all provided tasks are either COMPLETED or FAILED.
    :param tasks: List of GEE export task objects to monitor
    """
    print("Waiting for Google Earth Engine exports to complete...")
    for task in tasks:
        while task.active():
            print(f"  Task {task.id} is {task.status()['state']}... waiting 30 seconds.")
            time.sleep(30)
        
        status = task.status()
        if status['state'] == 'COMPLETED':
            print(f"  Task {task.id} COMPLETED.")
        else:
            print(f"  Task {task.id} FAILED: {status.get('error_message')}")


def download_from_gcs(bucket_name, source_blob_name, destination_file_name):
    """Downloads a blob from the GCS bucket to user's local machine.
    :param bucket_name: Name of the GCS bucket
    :param source_blob_name: Name of the blob in the GCS bucket (including path)
    :param destination_file_name: Local path where the file should be downloaded
    """
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(source_blob_name)

    os.makedirs(os.path.dirname(destination_file_name), exist_ok=True)
    blob.download_to_filename(destination_file_name)
    print(f"  Successfully downloaded to: {destination_file_name}")


def process_earthquake(eq, aoi, pairing_mode, job_list, resolution=90, mode='sar', optical_backend='copernicus'):
    """
    Process earthquake event and generate the necessary SLC pairs for InSAR processing.
    :param eq: dictionary containing earthquake data
    :param aoi: Area of Interest (AOI) as a GeoJSON file
    :param pairing_mode: 'all', 'sequential', or 'coseismic' for specifying desired SLC pairing
    :param job_list: True if the JSON objects are for HYP3 job submission, False otherwise
    :param resolution: Output resolution for the topsApp processing, default is 90m
    :param mode: 'sar' for SAR processing, 'optical' for optical processing
    :optical_backend: 'copernicus', 'element84' for source data file nomenclature (only applicable if mode is 'optical')
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

    elif aoi:
        print("AOI provided. Using the provided AOI...")
        # Load the AOI from the provided JSON file
        aoi = load_aoi_from_json(aoi)

    # Generate AOI or use the user-provided AOI
    else:
        aoi = make_aoi(coords) # Create AOI if not provided

    # Write AOI to a geojson file
    with open(f'{title}_AOI.geojson', 'w') as f:
        geojson.dump(mapping(aoi), f, indent=2)

    all_jobs = []
    all_features = []

    if mode == 'sar':
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
            SLCs = get_SLCs(flight_direction, path_number, frame_numbers, eq.get('time'), processing_mode='historic')
            isce_jobs = find_reference_and_secondary_pairs(SLCs, eq.get('time'), flight_direction, path_number, 
                                                           title, pairing_mode, job_list, resolution)
            all_jobs.append(isce_jobs)

    elif mode == 'optical':
        rupture_time = eq.get('time')
        rupture_dt = convert_time(rupture_time)
        
        start_search = (rupture_dt - timedelta(days=90)).strftime('%Y-%m-%dT%H:%M:%SZ')
        end_search = (rupture_dt + timedelta(days=90)).strftime('%Y-%m-%dT%H:%M:%SZ')

        if optical_backend == 'copernicus':
            print("Routing to Copernicus Public OData backend...")
            s2_scenes = search_copernicus_public(aoi, start_search, end_search)
            s2_jobs, s2_features = find_optical_pairs_copernicus(s2_scenes, rupture_time, title, aoi, job_list)
            
        elif optical_backend == 'element84':
            print("Routing to Element84 STAC backend...")
            s2_scenes = search_element84_stac(aoi, start_search, end_search)
            s2_jobs, f_dom, f_split = find_optical_pairs_element84(s2_scenes, rupture_time, title, aoi)
            s2_features = f_dom + f_split
            
        elif optical_backend == 'gee':
            print("Routing to Google Earth Engine backend...")
            
            # Fetch the bucket name from the environment
            gcs_bucket = os.getenv('COSEIS_GCS_BUCKET')
            if not gcs_bucket:
                print("Error: COSEIS_GCS_BUCKET environment variable is not set.")
                print("Please set it to your Google Cloud Storage bucket name (e.g., export COSEIS_GCS_BUCKET='my-gee-bucket').")
                return [], []

            try:
                ee.Initialize()
            except Exception as e:
                print("Earth Engine not authenticated. Run 'earthengine authenticate --auth_mode=notebook' in your terminal.")
                raise e

            # Define 30-day temporal windows
            rupture_dt = convert_time(rupture_time).replace(tzinfo=None)
            pre_start = (rupture_dt - timedelta(days=30)).strftime('%Y-%m-%d')
            pre_end = rupture_dt.strftime('%Y-%m-%d')
            
            post_start = rupture_dt.strftime('%Y-%m-%d')
            post_end = (rupture_dt + timedelta(days=30)).strftime('%Y-%m-%d')

            # Start GEE Tasks
            print(f"Generating Pre-Event Composite ({pre_start} to {pre_end})...")
            task_pre = export_gee_composite(aoi, pre_start, pre_end, title, "PRE", gcs_bucket)
            
            print(f"Generating Post-Event Composite ({post_start} to {post_end})...")
            task_post = export_gee_composite(aoi, post_start, post_end, title, "POST", gcs_bucket)

            # Wait for them to finish in Google Cloud Storage
            wait_for_gee_tasks([task_pre, task_post])

            # Download the composites from GCS to local directory
            print("Downloading composites from Google Cloud Storage bucket...")
            local_dir = os.path.join(root_dir, "GEE_Optical_Downloads", title)
            
            pre_filename = f"{task_pre.config['description']}.tif"
            post_filename = f"{task_post.config['description']}.tif"
            
            pre_local_path = os.path.join(local_dir, pre_filename)
            post_local_path = os.path.join(local_dir, post_filename)
            
            download_from_gcs(
                gcs_bucket, 
                f"COSEIS_Composites/{title}/{pre_filename}", 
                pre_local_path
            )
            download_from_gcs(
                gcs_bucket, 
                f"COSEIS_Composites/{title}/{post_filename}", 
                post_local_path
            )

            # Create the local manifest for AutoRIFT
            manifest_path = os.path.join(local_dir, f"{title}_autorift_manifest.json")
            manifest_payload = {
                "event_title": title,
                "backend": "Google Earth Engine",
                "pre_composite_path": pre_local_path,
                "post_composite_path": post_local_path,
                "status": "DOWNLOADED_READY_FOR_AUTORIFT"
            }
            
            with open(manifest_path, 'w') as f:
                json.dump(manifest_payload, f, indent=4)
                
            print(f"\nManifest written to: {manifest_path}")

            # Return empty lists since we bypass the traditional HYP3 job queue
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


def get_next_pass(AOI, timestamp_dir, satellite="sentinel-1"):
    """
    Get the next satellite pass over the given AOI.
    Uses the next_pass.py script to determine the next overpass.
    :param AOI: Shapely Polygon object representing the AOI
    :param satellite: Satellite name ('sentinel-1', 'sentinel-2', or 'landsat')
    :return: Next overpass time or None if an error occurs
    """
    import os
    import sys
    import next_pass
    next_pass_dir = os.path.dirname(next_pass.__file__)
    if next_pass_dir not in sys.path:
        sys.path.append(next_pass_dir)
    try:
        from utils import plot_maps
    except ImportError as e:
        print(f"Could not import plot_maps: {e}")
        return None, None, None, None
    
    from datetime import date

    min_lon, min_lat, max_lon, max_lat = AOI.bounds
    bbox = [str(min_lat), str(max_lat), str(min_lon), str(max_lon)]

    print("=========================================")
    print(f"Querying next-pass for all satellites over AOI...")
    print("=========================================")
    
    args = SimpleNamespace(bbox=bbox, sat="all", event_date=date.today(), look_back=14, cloudiness=False)

    try:
        result = next_pass.find_next_overpass(args, timestamp_dir)
    except Exception as e:
        print(f"Next pass error: {e}")
        return None, None, None, None
    
    result_s1 = result.get("sentinel-1") 
    result_nisar = result.get("nisar")
    
    s1_info = result_s1.get("next_collect_info", "No S1 info available.") if result_s1 else "No S1 info."
    nisar_info = result_nisar.get("next_collect_info", "No NISAR info available.") if result_nisar else "No NISAR info."
    
    default_map_file = timestamp_dir / "satellite_overpasses_map.html"
    
    # Generate S1-Only Map
    s1_map_path = timestamp_dir / "S1_overpass_map.html"
    try:
        # Pass None for NISAR
        plot_maps.make_overpasses_map(result_s1, None, None, None, args.bbox, timestamp_dir)
        if default_map_file.exists():
            os.rename(default_map_file, s1_map_path)
    except Exception as e:
        print(f"Could not generate S1-only map: {e}")
        s1_map_path = None

    # Generate S1 + NISAR Map
    nisar_map_path = timestamp_dir / "S1_NISAR_overpass_map.html"
    try:
        # Pass result_nisar into your plot_maps signature
        plot_maps.make_overpasses_map(result_s1, None, None, result_nisar, args.bbox, timestamp_dir)
        if default_map_file.exists():
            os.rename(default_map_file, nisar_map_path)
    except Exception as e:
        print(f"Could not generate S1+NISAR map: {e}")
        nisar_map_path = None

    return s1_info, nisar_info, s1_map_path, nisar_map_path


def main_forward(pairing_mode=None, resolution=90, do_processing=False, send_email_flag=False, mode='sar', optical_backend='copernicus'):
    """
    Runs the main query and processing workflow in forward processing mode.
    Used to produce co-seismic product for new earthquakes when new SLC data becomes available.
    :param pairing_mode: 'all', 'sequential', or 'coseismic' for specifying desired SLC pairing
    :param resolution: Output resolution for the topsApp processing, default is 90m
    :param do_processing: If True, runs the dockerized topsApp processing workflow after generating the JSONs. Default is False.
    :param send_email_flag: If True, sends an email alert after processing. Default is False.
    :param mode: 'sar' for SAR processing, 'optical' for optical processing
    :optical_backend: 'copernicus', 'element84' for source data file nomenclature (only applicable if mode is 'optical')
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
        
        # Initialize the tracking file if it doesn't exist
        if not os.path.exists(TRACKING_FILE):
            with open(TRACKING_FILE, 'w') as f:
                json.dump({}, f)

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

                    # Run next_pass to get the next S1 overpasses
                    s1_info, nisar_info, s1_map, nisar_map = get_next_pass(aoi, timestamp_dir)
                    
                    # Setup Github pages directory
                    docs_maps_dir = Path(os.getcwd()).parent / "docs" / "maps"
                    docs_maps_dir.mkdir(parents=True, exist_ok=True)
                    GITHUB_PAGES_BASE_URL = "https://cmspeed.github.io/coseis-sar"
                    
                    # Create a unique ID for the filenames so they aren't overwritten
                    unique_id = eq.get('id', datetime.now().strftime('%Y%m%d%H%M%S'))

                    # Route S1 Map
                    s1_map_url = ""
                    if s1_map and os.path.exists(s1_map):
                        new_s1_name = f"{title_snake}_{unique_id}_S1_overpass_map.html"
                        shutil.copy(s1_map, docs_maps_dir / new_s1_name)
                        s1_map_url = f"{GITHUB_PAGES_BASE_URL}/maps/{new_s1_name}"

                    # Route NISAR Map
                    nisar_map_url = ""
                    if nisar_map and os.path.exists(nisar_map):
                        new_nisar_name = f"{title_snake}_{unique_id}_S1_NISAR_overpass_map.html"
                        shutil.copy(nisar_map, docs_maps_dir / new_nisar_name)
                        nisar_map_url = f"{GITHUB_PAGES_BASE_URL}/maps/{new_nisar_name}"

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
                        
                        # Send S1-only to PRIMARY_RECIPIENTS
                        if PRIMARY_RECIPIENTS:
                            primary_body = (header_html + s1_section + get_button_html(s1_map_url) + footer_html).replace('\n', '')
                            send_email(subject=subject_text, body=primary_body, recipients=PRIMARY_RECIPIENTS)
                            print('=========================================')
                            print('S1-only email sent.')
                            print('=========================================')

                        # Send S1 + NISAR to TERTIARY_RECIPIENTS
                        if TERTIARY_RECIPIENTS:
                            tertiary_body = (header_html + s1_section + nisar_section + get_button_html(nisar_map_url) + footer_html).replace('\n', '')
                            send_email(subject=subject_text, body=tertiary_body, recipients=TERTIARY_RECIPIENTS)
                            print('=========================================')
                            print('S1+NISAR email sent.')
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

        # Check ASF DAAC for available SLCs for pending earthquakes
        check_tracker_for_updates(do_processing, send_email_flag)
        
    finally:
        if os.path.exists(lock_file):
            os.remove(lock_file)


def main_historic(start_date, end_date=None, aoi=None, pairing_mode=None, job_list=False, resolution=90, mode='sar', optical_backend='copernicus'):
    """
    Runs the main query and processing workflow in historic processing mode.
    Used to produce 'pre-seismic', 'co-seismic', and 'post-seismic' displacement products for historic earthquakes.
    Pre-seismic and post-seismic data are generated for 90 days before and 30 days after the event.
    A single date or date range can be provided for processing.
    :param start_date: The query start date in YYYY-MM-DD format
    :param end_date: The query end date in YYYY-MM-DD format (Optional)
    :param aoi: The path to a JSON file representing the area of interest (AOI) (Optional)
    :param pairing_mode: 'all', 'sequential', or 'coseismic' for specifying desired SLC pairing.
    :param job_list: If True, create a list of jobs in HYP3 format for cloud processing.
    :param resolution: Output resolution for the topsApp processing, default is 90m
    :param mode: 'sar' for SAR processing, 'optical' for optical processing
    :optical_backend: 'copernicus', 'element84' for source data file nomenclature (only applicable if mode is 'optical')
    """
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
        eq_sig = check_significance(earthquakes, start_date, end_date, mode='historic')

        if eq_sig is not None:
            jobs_dict = []
            master_scene_features = []
            earthquake_infos = [] 

            for eq in eq_sig:
                try:
                    eq_jsons, eq_features = process_earthquake(eq, aoi, pairing_mode, job_list, resolution, mode, optical_backend)

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
                feature_filename = f'all_selected_scenes_{mode}_{optical_backend}_{current_time}.geojson'
                fc = {"type": "FeatureCollection", "features": master_scene_features}
                with open(feature_filename, 'w') as f:
                    json.dump(fc, f, indent=2)
                print(f"Saved master scene footprints to {feature_filename}")
        else:
            print(f"No significant earthquakes found betweeen {start_date} and {end_date}.")


if __name__ == "__main__":
    """
    Run the main function based on the input arguments provided, in either 'historic' or 'forward' processing mode.
    Historic processing can be done for a single date or range of dates, with the option to specify the SLC pairing mode.
    Forward processing is used to generate co-seismic displacement products for new earthquakes.
    Example usage for historic processing: 
      python coseis_sar.py --historic --dates 2021-08-14 --pairing all
      python coseis_sar.py --historic --dates 2021-08-14 2021-09-07 --pairing all
      python coseis_sar.py --historic --dates 2014-06-14 2025-02-12 --pairing coseismic (all coseismic pairs from beginning of S1 data to 2025-02-12)
      python coseis_sar.py --historic --dates 2014-06-14 2025-02-12 --pairing coseismic --job_list --resolution 30 (only produce the job list for HYP3 processing, don't run any jobs locally) 
    Example usage for forward processing: 
      python coseis.py --forward
    """
    parser = argparse.ArgumentParser(
        description="Run historic or forward processing based on input arguments."
    )

    parser.add_argument("--historic", action="store_true", help="Run historic processing.")
    parser.add_argument("--forward", action="store_true", help="Run forward processing.")
    parser.add_argument("--dates", nargs="+", help="Provide one or two dates in YYYY-MM-DD format for historic processing.")
    parser.add_argument("--aoi", help="Specify a path to a json file representing the area of interest (AOI).")
    parser.add_argument("--pairing", choices=["all", "sequential", "coseismic"], help="Specify the SLC pairing mode. Required for historic processing.")
    parser.add_argument("--job_list", action="store_true", help="Create a list of jobs in HYP3 format for cloud processing.")
    parser.add_argument("--resolution", type=int, default=90, help="Output resolution for topsApp processing in meters. Default is 90m.")
    parser.add_argument("--mode", choices=["sar", "optical"], default="sar", help="Processing mode: 'sar' (Sentinel-1) or 'optical' (Landsat/Sentinel-2). Default is sar.")
    parser.add_argument("--do_processing", action="store_true", help="Execute local topsApp processing.")
    parser.add_argument("--send_email", action="store_true", help="Send email notifications.")
    parser.add_argument("--optical_backend", choices=["copernicus", "element84", "gee"], default="copernicus", help="Specify the optical data provider if mode is optical. Default is copernicus.")

    args = parser.parse_args()

    if args.mode == 'sar' and '--optical_backend' in sys.argv:
        print("Error: --optical_backend can only be used when --mode is 'optical'.")
        parser.print_help()
        exit(1)

    if args.historic:
        if args.send_email or args.do_processing:
            print("Error: --send_email and --do_processing are only supported in --forward mode.")
            parser.print_help()
            exit(1)
            
        if not args.dates:
            print("Error: --dates is required when using --historic mode.")
            parser.print_help()
            exit(1)

        if len(args.dates) > 2:
            print("Error: --dates should have at most two values (start_date [end_date]).")
            parser.print_help()
            exit(1)

        start_date = args.dates[0]
        end_date = args.dates[1] if len(args.dates) == 2 else None

        aoi = args.aoi if args.aoi else None

        if not args.pairing:
            print("Error: --pairing is required when using --historic mode. Options: 'all', 'sequential', 'coseismic'.")
            parser.print_help()
            exit(1)

        main_historic(start_date, end_date, aoi=aoi, pairing_mode=args.pairing, job_list=args.job_list, resolution=args.resolution, mode=args.mode, optical_backend=args.optical_backend)

    elif args.forward:
        if args.job_list:
            print("Error: --job_list is only supported in --historic mode.")
            parser.print_help()
            exit(1)
            
        if args.dates:
            print("Error: --dates cannot be used with --forward mode.")
            parser.print_help()
            exit(1)

        if not args.pairing:
            print("Error: --pairing is required when using --forward mode. Options: 'all', 'sequential', 'coseismic'.")
            parser.print_help()
            exit(1)

        if args.aoi:
            print("Error: --aoi cannot be used with --forward mode.")
            parser.print_help()
            exit(1)
            
        if not args.do_processing and not args.send_email:
            print("Warning: Running --forward without --do_processing or --send_email. The script will only update tracking files.")

        main_forward(args.pairing, args.resolution, args.do_processing, args.send_email, mode=args.mode, optical_backend=args.optical_backend)