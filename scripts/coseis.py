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
import csv
from shapely import wkt
from shapely.geometry import mapping, shape, box, Point, Polygon, LineString, MultiLineString, MultiPolygon
from shapely.ops import unary_union
import subprocess
from datetime import datetime, timedelta, timezone
from collections import defaultdict
from itertools import combinations
import logging
import yagmail
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
root_dir = os.path.join(os.getcwd(), "data")  # Defaults to ./data; change is

# Global variables
OPTICAL_CLOUD_THRESHOLD = 20.0  # Maximum cloud cover percentage for optical data

PENDING_EARTHQUAKES_FILE = "pending_earthquakes.json"
PRIMARY_RECIPIENTS = ['cole.speed@jpl.nasa.gov', 'cole.speed@yahoo.com',
                      'mary.grace.p.bato@jpl.nasa.gov', 'mgbato@gmail.com',
                      'eric.j.fielding@jpl.nasa.gov', 'emre.havazli@jpl.nasa.gov',
                      'bryan.raimbault@jpl.nasa.gov', 'karen.an@jpl.nasa.gov',
                      'ines.fenni@jpl.nasa.gov', 'alexander.handwerger@jpl.nasa.gov',
                      'brett.a.buzzanga@jpl.nasa.gov', 'dmelgarm@uoregon.edu', 'msolares@uoregon.edu'
                      ]

SECONDARY_RECIPIENTS = ['cole.speed@jpl.nasa.gov', 'cole.speed@yahoo.com']

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
                if (magnitude >= 6.0) and (alert in alert_list) and (depth <= 40.0) and within_Coastline_buffer:
                    significant_earthquakes.append(earthquake)

    # Base significance on magnitude, depth, and distance from land for forward-looking data
    if mode =='forward':
        for earthquake in earthquakes:
            magnitude = earthquake.get('mag')
            depth = earthquake.get('coordinates', [])[2] if earthquake.get('coordinates') else None
            within_Coastline_buffer = withinCoastline(earthquake, coastline)
            if all(var is not None for var in (magnitude, depth)):
                if (magnitude >= 6.0) and (depth <= 40.0) and within_Coastline_buffer:
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

            # Convert the GeoJSON data to a GeoDataFrame for visualization
            frame_dataframe = gpd.GeoDataFrame.from_features(data['features'], crs='EPSG:4326')

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
                    return None  # Or raise, or return, depending on your flow

        except Exception as e:
            print(f"Unexpected error while processing ASF DAAC response: {e}")
            return None


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
    start_date = rupture_date - timedelta(days=90)  # 90 days before the earthquake
    start_date= start_date.replace(hour=0, minute=0, second=0)

    # Set the end date based on the processing mode
    if processing_mode == 'historic':
        end_date = rupture_date + timedelta(days=30)    # 30 days after the earthquake
        end_date = end_date.replace(hour=23, minute=59, second=59)
    elif processing_mode == 'forward':
        today = datetime.now()
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
                    date = isoparse(start_time).date().isoformat()
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


def find_optical_pairs(optical_scenes, rupture_time, title, aoi_polygon, job_list=True):
    """
    Generate Optical Pairs grouped by Relative Orbit.
    Prioritizes: 1. True Coverage Area %, 2. Cloud Cover, 3. Time.
    Returns: (jobs, scene_features)
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
            
            # Split by Platform (S2A/S2B/S2C)
            platforms_on_date = defaultdict(list)
            for s in scenes_on_date:
                platforms_on_date[s['platform']].append(s)
            
            for platform, scenes in platforms_on_date.items():
                unique_scenes = []
                seen_tiles = set()
                scenes.sort(key=lambda x: x['granule_id'], reverse=True)
                
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

        if not pre_candidates or not post_candidates:
            continue

        # Sort: Max Coverage > Low Cloud > Low Delta
        pre_candidates.sort(key=lambda x: (-x['coverage'], x['cc'], x['delta']))
        best_pre = pre_candidates[0]
        
        post_candidates.sort(key=lambda x: (-x['coverage'], x['cc'], x['delta']))
        best_post = post_candidates[0]
        
        print(f"  Orbit {orbit_id}: Pre={best_pre['date']} ({best_pre['coverage']:.1f}% Cov), Post={best_post['date']} ({best_post['coverage']:.1f}% Cov)")

        # Build GeoJSON Features
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
            # Job Name includes metadata
            job_name = f"{title}-S2-R{orbit_id}-{best_pre['date']}_{best_post['date']}"
            
            # Strip .SAFE extension from IDs
            primary_ids = [s['granule_id'].replace('.SAFE', '') for s in best_pre['scenes']]
            secondary_ids = [s['granule_id'].replace('.SAFE', '') for s in best_post['scenes']]
            
            job = {
                "name": job_name,
                "job_type": "AUTORIFT",
                "job_parameters": {
                    "reference": primary_ids,
                    "secondary": secondary_ids
                }
            }
            jobs.append(job)
            
    return jobs, scene_features

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
    rupture_date_dt = datetime.strptime(rupture_date, '%Y-%m-%d')
    
    # Reformatting for dictionary keys for later use    
    flight_direction = 'A' if flight_direction == 'ASCENDING' else 'D'
    path_number = f"{int(path_number):03}" 

    # Pair SLCs by date 
    slc_by_date = {}
    for slc in SLCs:
        date_obj = datetime.strptime(slc['date'], "%Y-%m-%d")
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


def filter_S1C(job):
    """Return True if the job does NOT contain any S1C granules in primary or secondary."""
    job_params = job.get("job_parameters", {})
    granules = job_params.get("granules", [])
    secondary = job_params.get("secondary_granules", [])

    # Reject if any granule in either list starts with S1C_
    return not any(g.startswith("S1C_") for g in granules + secondary)


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


def run_dockerized_topsApp(json_data):
    """
    Run dockerized topsApp InSAR processing workflow using the provided JSON data.
    Outputs are added to the root dir + an extension for each pair.
    :param json_data: JSON object containing the parameters for dockerized topsApp
    """
    # Extract the parameters from the JSON data
    reference_scenes = json_data['reference-scenes']
    secondary_scenes = json_data['secondary-scenes']
    frame_id = json_data['frame-id']
    estimate_ionosphere_delay = json_data['estimate-ionosphere-delay']
    esd_coherence_threshold = json_data['esd-coherence-threshold']
    compute_solid_earth_tide = json_data['compute-solid-earth-tide']
    goldstein_filter_power = json_data['goldstein-filter-power']
    output_resolution = json_data['output-resolution']
    unfiltered_coherence = json_data['unfiltered-coherence']
    dense_offsets = json_data['dense-offsets']
    process = 'coseis_sar'
    
    # Construct the command
    command = [
        "isce2_topsapp",
        "--reference-scenes", " ".join(reference_scenes),
        "--secondary-scenes", " ".join(secondary_scenes),
        "--frame-id", str(frame_id),
        "--estimate-ionosphere-delay", str(estimate_ionosphere_delay),
        "--esd-coherence-threshold", str(esd_coherence_threshold),
        "--compute-solid-earth-tide", str(compute_solid_earth_tide),
        "--goldstein-filter-power", str(goldstein_filter_power),
        "--output-resolution", str(output_resolution),
        "--unfiltered-coherence", str(unfiltered_coherence),
        "--dense-offsets", str(dense_offsets),
        "++process", process
    ]
    
    print('=========================================')
    print("Running the dockerized topsApp InSAR processing workflow...")
    print('=========================================')

    # Run the command
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        print("Command executed successfully!")
        print("Output:\n", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error occurred while running the command.")
        print("Error Output:\n", e.stderr)
    return


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


def send_email(subject, body, attachment=None, recipients=None):
    """
    Send an email with the earthquake information to a specified list of recipients.
    :param subject: Email subject
    :param body: Email body content
    :param attachment: Path to an optional file to attach
    :param recipients: List of email addresses to send the email to
    """
    if recipients is None:
        recipients = PRIMARY_RECIPIENTS
    
    GMAIL_USER = 'aria.hazards.jpl@gmail.com'
    GMAIL_PSWD = os.environ['GMAIL_APP_PSWD']
    yag = yagmail.SMTP(GMAIL_USER, GMAIL_PSWD)

    yag.send(
             bcc=recipients,
             subject=subject,
             contents=[body],
             attachments=[attachment] if attachment else None
             )
    return


def process_earthquake(eq, aoi, pairing_mode, job_list, resolution=90, mode='sar'):
    """
    Process earthquake event and generate the necessary SLC pairs for InSAR processing.
    :param eq: dictionary containing earthquake data
    :param aoi: Area of Interest (AOI) as a GeoJSON file
    :param pairing_mode: 'all', 'sequential', or 'coseismic' for specifying desired SLC pairing
    :param job_list: True if the JSON objects are for HYP3 job submission, False otherwise
    :param resolution: Output resolution for the topsApp processing, default is 90m
    :param mode: 'sar' for SAR processing, 'optical' for optical processing
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
        geojson.dump(aoi, f, indent=2)

    all_jobs = []

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

        s2_scenes = search_copernicus_public(aoi, start_search, end_search)
        
        # UPDATED: Unpack both return values
        s2_jobs, s2_features = find_optical_pairs(s2_scenes, rupture_time, title, aoi, job_list)
        
        if s2_jobs:
            print(f"Generated {len(s2_jobs)} Sentinel-2 jobs.")
            all_jobs.append(s2_jobs)
            
            # NEW: Save the footprints to GeoJSON
            if s2_features:
                fc = {"type": "FeatureCollection", "features": s2_features}
                scene_filename = f"{title}_selected_scenes.geojson"
                with open(scene_filename, 'w') as f:
                    json.dump(fc, f, indent=2)
                print(f"Saved selected scene footprints to {scene_filename}")
        else:
            print("No optical pairs found.")

    return all_jobs


def get_next_pass(AOI, timestamp_dir, satellite="sentinel-1"):
    """
    Get the next satellite pass over the given AOI.
    Uses the next_pass.py script to determine the next overpass.
    :param AOI: Shapely Polygon object representing the AOI
    :param satellite: Satellite name ('sentinel-1', 'sentinel-2', or 'landsat')
    :return: Next overpass time or None if an error occurs
    """
    import next_pass
    import plot_maps
    from datetime import date, datetime

    min_lon, min_lat, max_lon, max_lat = AOI.bounds
    bbox = [min_lat, max_lat, min_lon, max_lon]

    print("=========================================")
    print(f"Querying next pass for {satellite} over AOI...")
    print(f"BBOX: {min_lat}, {max_lat}, {min_lon}, {max_lon}")
    print("=========================================")
    
    args = SimpleNamespace(bbox=bbox, sat=satellite, event_date=date.today)

    result = next_pass.find_next_overpass(args)
    result_s1 = result["sentinel-1"] 
    result_s2 = result["sentinel-2"]
    result_l = result["landsat"]
    plot_maps.make_overpasses_map(result_s1, result_s2, result_l, args.bbox, timestamp_dir)

    # loop over results and display only missions that were requested
    for _, mission_result in result.items():
        if mission_result:
            s1_next_collect_info = mission_result.get("next_collect_info",
                                     "No collection info available.")
            return s1_next_collect_info
        else:
            print("No next pass data available.")
            return None


def load_pending_earthquakes() -> List[Dict[str, Any]]:
    """Loads the list of pending earthquakes from the JSON file."""
    if not os.path.exists(PENDING_EARTHQUAKES_FILE):
        return []
    with open(PENDING_EARTHQUAKES_FILE, "r") as f:
        return json.load(f)


def save_pending_earthquakes(earthquakes: List[Dict[str, Any]]):
    """Saves the list of pending earthquakes to the JSON file."""
    with open(PENDING_EARTHQUAKES_FILE, "w") as f:
        json.dump(earthquakes, f, indent=4)


def store_pending_earthquake(eq: dict, aoi: Polygon):
    """Stores a newly detected significant earthquake for later processing."""
    pending_list = load_pending_earthquakes()
    pending_list.append({
        "id": eq.get('id'),
        "title": to_snake_case(eq.get('title')),
        "time": eq.get('time'),
        "coords": eq.get('coordinates'),
        "aoi": mapping(aoi)  # Store AOI geometry as a GeoJSON-compatible dictionary
    })
    save_pending_earthquakes(pending_list)
    print(f"Stored pending earthquake '{eq.get('title')}' for future checks.")


def check_pending_slcs(pairing_mode: str, resolution: int):
    """
    Checks the pending earthquakes queue for available SLCs to form co-seismic pairs.
    If pairs are found, it generates a job list and sends an email.
    """
    pending_list = load_pending_earthquakes()
    if not pending_list:
        print("No pending earthquakes to check for SLCs.")
        return

    print('=========================================')
    print("Checking for available co-seismic SLCs for pending earthquakes...")
    print('=========================================')

    updated_pending_list = []
    completed_jobs_list = []

    for eq in pending_list:
        event_id = eq['id']
        title = eq['title']
        time = eq['time']
        
        # FFM check and AOI update 
        ffm_url = get_ffm_geojson_url(event_id)
        if ffm_url:
            print(f"FFM found for {title}. Updating AOI...")
            aoi = load_aoi_from_json(ffm_url)
            eq['aoi'] = mapping(aoi) # Update the AOI in the dictionary
        else:
            # If no FFM, use the stored AOI
            aoi_geom = eq['aoi']
            aoi = shape(aoi_geom)

        print(f"Processing '{title}' with AOI: {aoi}")
        
        path_frame_numbers, _ = get_path_and_frame_numbers(aoi, time)
        
        found_pairs_for_eq = False
        for (flight_direction, path_number), frame_numbers in path_frame_numbers.items():
            frame_numbers = list(set(fn[0] for fn in frame_numbers))
            SLCs = get_SLCs(flight_direction, path_number, frame_numbers, time, processing_mode='forward')
            
            # Find coseismic pairs only
            isce_jobs = find_reference_and_secondary_pairs(SLCs, time, flight_direction, path_number,
                                                           title, pairing_mode='coseismic', job_list=True, resolution=resolution)

            if isce_jobs:
                # Add found jobs to the master list
                if filter_S1C(isce_jobs[0]):
                    completed_jobs_list.extend(isce_jobs)
                    found_pairs_for_eq = True
        
        if found_pairs_for_eq:
            # Earthquake has a completed coseismic pair, so don't re-add to the pending list
            print(f"Co-seismic pair found for '{title}'. Job will be generated.")
        else:
            # No coseismic pair found yet, so keep it in the pending list (with updated AOI)
            print(f"No co-seismic pair found yet for '{title}'. Will check again later.")
            updated_pending_list.append(eq)

    save_pending_earthquakes(updated_pending_list)

    if completed_jobs_list:
        current_time = datetime.now(timezone.utc).strftime("%Y-%m-%d_%H-%M-%S_UTC")
        job_list_filename = f'jobs_list_{current_time}.json'
        with open(job_list_filename, 'w') as f:
            json.dump(completed_jobs_list, f, indent=4)
        print(f"Generated job list: {job_list_filename}")

        # Send second email
        send_email(
            subject=f"COSEISMIC JOB LIST READY",
            body=(
                f"A job list has been generated for one or more co-seismic events.\n\n"
                f"The job list file, '{job_list_filename}', is ready and stored on the server.\n\n"
                f"This is an automated message. Please do not reply."
            ),
            recipients=SECONDARY_RECIPIENTS
        )
        print('=========================================')
        print('Job list alert emailed to secondary recipients.')
        print('=========================================')


def main_forward(pairing_mode=None):
    """
    Runs the main query and processing workflow in forward processing mode.
    Used to produce co-seismic product for new earthquakes when new SLC data becomes available.
    :param pairing_mode: 'all', 'sequential', or 'coseismic' for specifying desired SLC pairing
    """

    print('=========================================')
    print("Running cronjob to check for new earthquakes...")
    print('=========================================')
    
    # Check for New Earthquakes
    geojson_data = check_for_new_data(USGS_api_hourly)

    start_date = datetime.now().strftime('%Y-%m-%d')
    current_time = datetime.now(timezone.utc).strftime("%Y-%m-%d at %H:%M:%S UTC")

    if geojson_data:
        # Parse GeoJSON and create variables for each feature's properties
        earthquakes = parse_geojson(geojson_data)
        eq_sig = check_significance(earthquakes, start_date, end_date=None, mode = 'forward')

        if eq_sig is not None:
            for eq in eq_sig:
                # Check for duplicate entry in the pending queue
                pending_list = load_pending_earthquakes()
                if any(p['id'] == eq.get('id') for p in pending_list):
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
                next_pass_info = get_next_pass(aoi, timestamp_dir, satellite="sentinel-1")
                
                # Rename sentinel-1 overpass and opera granule maps
                original_filename = timestamp_dir / "satellite_overpasses_map.html"
                s1_map_filename = timestamp_dir / f"{title_snake}_Sentinel-1_Next_Overpasses.html"

                if os.path.exists(original_filename):
                    os.rename(original_filename, s1_map_filename)
                    print(f"Renamed to {s1_map_filename}")
                else:
                    print(f"Expected file {original_filename} not found.")

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

                send_email(
                    subject=f"{message_dict['title']}",
                    body=(
                        f"{message_dict['title']} at {message_dict['time']} UTC\n"
                        f"------------------------------------------------------------------------------------------------------------------------\n"
                        f"Epicenter coordinates (lat, lon): ({message_dict['coordinates'][1]}, {message_dict['coordinates'][0]})\n"
                        f"Depth: {message_dict['depth']} km\n"
                        f"For more details, visit the USGS Earthquake Hazard Portal page for this event: {message_dict['url']}\n"
                        f"------------------------------------------------------------------------------------------------------------------------\n"
                        f"Next acquisition times and relative orbits for Sentinel-1 frames intersecting the earthquake AOI:\n"
                        f"{next_pass_info}\n"
                        f"------------------------------------------------------------------------------------------------------------------------\n"
                        f"View an interactive map of intersecting Sentinel-1 orbits by downloading the attachment and launching it your web browser.\n"
                        f"------------------------------------------------------------------------------------------------------------------------\n"
                        f"This is an automated message. Please do not reply. For product-specific inquiries contact Dr. Cole Speed (<a href=\"mailto:cole.speed@jpl.nasa.gov\">cole.speed@jpl.nasa.gov</a>) and Dr. Grace Bato (<a href=\"mailto:bato@jpl.nasa.gov\">bato@jpl.nasa.gov</a>)."
                    ),
                    attachment=s1_map_filename,
                    recipients=PRIMARY_RECIPIENTS
                )

                print('=========================================')
                print('Alert emailed to recipients.')
                print('=========================================')

                # Store the new earthquake for later checks
                store_pending_earthquake(eq, aoi)
        else:
            print(f"No new significant earthquakes found as of {current_time}.")

    # Check ASF DAAC for available SLCs for poending earthquakes
    check_pending_slcs(pairing_mode, args.resolution)


def main_historic(start_date, end_date = None, aoi = None, pairing_mode = None, job_list = False, resolution=90, mode='sar'):
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
    """
    if start_date and not end_date:
        print('=========================================')
        print(f"Running historic processing in single-date mode for date: {start_date}")
        print('=========================================')
        # Fetch GeoJSON data from the USGS Earthquake Hazard Portal for a single date
        geojson_data = get_historic_earthquake_data_single_date(USGS_api_alltime, str(start_date))

    elif start_date and end_date:
        print('=========================================')
        print(f"Running historic processing in date range mode for dates: {start_date} to {end_date}")
        print('=========================================')
        # Fetch GeoJSON data from the USGS Earthquake Hazard Portal for the date range provided
        geojson_data = get_historic_earthquake_data_date_range(USGS_api_alltime, str(start_date), str(end_date))

    if geojson_data:
        # Parse GeoJSON and create variables for each feature's properties
        earthquakes = parse_geojson(geojson_data)
        eq_sig = check_significance(earthquakes, start_date, end_date, mode='historic')

        if eq_sig is not None:
            jobs_dict = []
            for eq in eq_sig:
                try:
                    eq_jsons = process_earthquake(eq, aoi, pairing_mode, job_list, resolution, mode)
                except Exception as e:
                    print(f"Error processing {eq['title']}: {e}")
                    continue

                if job_list:    # produce the list of jobs needed for HYPE3 processing, but do not run any jobs
                    for i, eq_json in enumerate(eq_jsons):
                        for j, json_data in enumerate(eq_json):
                            if mode == 'sar':
                                if filter_S1C(json_data):  # apply filter to remove jobs contain Sentinel-1C granules
                                    jobs_dict.append(json_data)
                            else:
                                jobs_dict.append(json_data)

                else:   # run dockerizedTopsApp locally
                    
                    if mode == 'sar':
                        # Create directories for each group of SLCs based on the JSON data provided
                        dirnames, total = create_directories_from_json(eq_jsons, root_dir)
                        
                        for i, eq_json in enumerate(eq_jsons):
                            for j, json_data in enumerate(eq_json):
                                
                                # Navigate to working directory for job
                                working_dir = dirnames[i][j]
                                os.chdir(working_dir)

                                # Add json to the working directory 
                                with open('isce_params.json', 'w') as f:
                                    json.dump(json_data, f, indent=4)
                                
                                print(f'working directory: {os.getcwd()}')
                                print(f'Running dockerized topsApp for dates {json_data["secondary-date"]} to {json_data["reference-date"]}...')
                                try:
                                    run_dockerized_topsApp(json_data)
                                except:
                                    print('Error running dockerized topsApp')
                                    continue
                    else:
                        print("Optical processing mode selected. Running jobs locally is not implemented yet.")

            if job_list and jobs_dict:  # Write only once after all earthquakes are processed
                current_time = datetime.now(timezone.utc).strftime("%Y-%m-%d_%H-%M-%S_UTC")
                with open(f'jobs_list_{current_time}.json', 'w') as f:
                    json.dump(jobs_dict, f, indent=4)
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

    args = parser.parse_args()

    if args.historic and args.forward:
        print("Error: You cannot specify both --historic and --forward.")
        parser.print_help()
        exit(1)

    if args.historic:
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

        main_historic(start_date, end_date, aoi=aoi, pairing_mode=args.pairing, job_list=args.job_list, resolution=args.resolution, mode=args.mode)

    elif args.forward:
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

        main_forward(args.pairing)

    else:
        print("Error: Please specify either --historic or --forward to run the processing workflow.")
        parser.print_help()
        exit(1)