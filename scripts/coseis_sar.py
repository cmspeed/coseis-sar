import argparse
import asf_search as asf
import requests
import json
import geojson
from shapely.geometry import mapping, shape, Point, Polygon, LineString, MultiLineString, MultiPolygon
from shapely.ops import unary_union
import subprocess
from datetime import datetime, timedelta, timezone
from collections import defaultdict
from itertools import combinations
import logging
import re
import os

# Set logging level to WARNING to suppress DEBUG and INFO logs
logging.basicConfig(level=logging.WARNING)

# API Endpoints
USGS_api_hourly = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_hour.geojson"  # USGS Earthquake API - Hourly
USGS_api_30day = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_month.geojson"  # USGS Earthquake API - Monthly
USGS_api_alltime = "https://earthquake.usgs.gov/fdsnws/event/1/query" # USGS Earthquake API - All Time
coastline_api = "https://raw.githubusercontent.com/OSGeo/PROJ/refs/heads/master/docs/plot/data/coastline.geojson" # Coastline API
ASF_DAAC_API = "https://api.daac.asf.alaska.edu/services/search/param"
root_dir = "/u/trappist-r0/colespeed/work/coseis/earthquakes/"

def get_historic_earthquake_data_single_date(eq_api, input_date):
    """
    Fetch data from the USGS Earthquake Portal for a single date and returns it as a GeoJSON object.
    The data returned will depend on the parameters included with the API request.
    :param: eq_api - USGS API endpoint
    :param: input_date - date in the format 'YYYY-MM-DD'
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

        # Save to a GeoJSON file
        output_file = f"historic_earthquakes_{input_date}.geojson"
        with open(output_file, "w") as f:
            json.dump(earthquakes, f, indent=2)
        
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
    :param: eq_api - USGS API endpoint
    :param: start_date - start date in the format 'YYYY-MM-DD'
    :param: end_date - end date in the format 'YYYY-MM-DD'
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
            "minmagnitude": 6.0,
            "maxdepth": 30.0
        }

        # Fetch data from the USGS Earthquake API
        response = requests.get(eq_api, params=params)
        response.raise_for_status()  # Raise error if request fails
        
        # Parse the response as GeoJSON
        earthquakes = geojson.loads(response.text)

        # Save to a GeoJSON file
        output_file = f"historic_earthquakes_{start_date}_to_{end_date}.geojson"
        with open(output_file, "w") as f:
            json.dump(earthquakes, f, indent=2)
        
        return earthquakes

    except requests.RequestException as e:
        print(f"Error accessing primary API: {e}")
        return None
    except geojson.GeoJSONDecodeError as e:
        print(f"Error parsing GeoJSON data: {e}")
        return None

def check_for_new_data(eq_api):
    """
    Fetch data from the USGS Earthquake Portal and returns it as a GeoJSON object.
    The data returned is updated hourly and includes all earthquakes occuring during that period.
    :param: eq_api - USGS API endpoint
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
    Fetch coastline data from https://raw.githubusercontent.com/OSGeo/PROJ/refs/heads/master/docs/plot/data/coastline.geojson and return it as a GeoJSON object.
    The data returned will be in the form of a MultiPolygon. The data are also written to a file: "coastline_buffered.geojson".
    :param: coastline_api - API endpoint for the coastline data
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
    :param: geojson_data - GeoJSON object containing earthquake data
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
        
        # Append the dictionary to the list
        earthquakes.append(feature_dict)
    return earthquakes

def withinCoastline(earthquake, coastline):
    """
    Determine if earthquake epicenter is within 0.5 decimal degrees (~55 km) of the coastline.
    This is one filtering parameter to determine if an earthquake is "significant" within the scope of this project.
    :param: earthquake - dictionary containing earthquake data
    :param: coastline - shapely Polygon object representing the coastline
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

def check_significance(earthquakes):
    """
    Check the significance of each earthquake based on its 
    (1) magnitude, (2) USGS alert level, (3) depth, and (4) distance from land.
    Criteria: Magnitude >= 6.0, USGS alert level = ['yellow','orange','red], and Depth <= 30.0 km, and within 0.5 degrees (~55 km) of the coastline.
    :param: earthquakes - list of dictionaries containing earthquake data
    :return: List of dictionaries containing significant earthquakes
    """
    print('=========================================')
    print("Checking for significant earthquakes...")
    print('=========================================')

    significant_earthquakes = []
    alert_list = ['yellow', 'orange', 'red']
    coastline = get_coastline(coastline_api)

    for earthquake in earthquakes:
        magnitude = earthquake.get('mag')
        alert = earthquake.get('alert')
        depth = earthquake.get('coordinates', [])[2] if earthquake.get('coordinates') else None
        within_Coastline_buffer = withinCoastline(earthquake, coastline)
        if all(var is not None for var in (magnitude, alert, depth)):
            if (magnitude >= 6.0) and (alert in alert_list) and (depth <= 30.0) and within_Coastline_buffer:
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
        significant_earthquakes_to_geojson(significant_earthquakes)
        return significant_earthquakes
    else:
        return None

def significant_earthquakes_to_geojson(significant_earthquakes):
    """
    Write the significant earthquakes to a GeoJSON file, namely: "significant_earthquakes_full_record.geojson"
    :param: significant_earthquakes - list of dictionaries containing significant earthquake data
    """
    geojson_features = []

    # Loop through each earthquake and create a GeoJSON feature
    for eq in significant_earthquakes:
        # Create a GeoJSON feature for each earthquake
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Point",  # Using Point geometry for each earthquake
                "coordinates": eq["coordinates"][:2]  # Only take the first two values: longitude and latitude
            },
            "properties": {
                "magnitude": eq["mag"],
                "place": eq["place"],
                "time": eq["time"],
                "alert": eq["alert"],
                "url": eq["url"],
                "depth": eq["coordinates"][2]  # Adding the depth (third value in coordinates) as a property
            }
        }
        geojson_features.append(feature)

    # Create the final GeoJSON structure
    geojson_data = {
        "type": "FeatureCollection",
        "features": geojson_features
    }

    # Save the GeoJSON data to a file
    with open('significant_earthquakes_full_record.geojson', 'w') as f:
        geojson.dump(geojson_data, f)

    return
        
def make_aoi(coordinates):
    """
    Create an Area of Interest (AOI) polygon based on the given coordinates.
    The AOI is a square with a side length of 1 degree (~111 km) centered on the earthquake's epicenter.
    The AOI is written to a GeoJSON file, namely: "AOI.geojson".
    :param: coordinates - list containing the longitude and latitude of the earthquake's epicenter
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
        (X - half_side, Y - half_side),  # bottom-left
        (X + half_side, Y - half_side),  # bottom-right
        (X + half_side, Y + half_side),  # top-right
        (X - half_side, Y + half_side)   # top-left
    ]
    
    # Create the square polygon
    AOI = Polygon(square_coords)

    # Write to a geojson file
    with open('AOI.geojson', 'w') as f:
        geojson.dump(AOI, f, indent=2)

    print(f"Area of Interest (AOI) created: {AOI}")
    print('=========================================')
    return AOI

def convert_time(time):
    """
    Convert the given Unix timestamp in milliseconds to a UTC datetime object.
    :param time_ms: Unix timestamp in milliseconds (int or float)
    :return: Datetime object in UTC in this format: 'YYYY-MM-DDTHH:MM:SS'
    """
    timestamp_s = time / 1000 # Timestamp in seconds
    dt = datetime.fromtimestamp(timestamp_s, tz=timezone.utc) # Convert to datetime object in UTC
    dt = dt.replace(microsecond=0) # Remove microseconds
    return dt

def get_path_and_frame_numbers(AOI, time):
    """
    Query the ASF DAAC API for SLC data intersecting the Area of Interest (AOI) over the previous 24 days.
    This ensures all possible intersecting SLC fileIDs are returned for the given AOI.
    :param: AOI - Shapely Polygon object representing the Area of Interest
    :param: time - Unix timestamp representing the earthquake's origin time
    :return: Dictionary containing the path and frame numbers for each *unique* intersecting SLC.
    """
    # Establish the date range for the query
    rupture_date = convert_time(time)
    start_date = rupture_date - timedelta(days=24) # 24 days before the earthquake
    start_date = start_date.replace(hour=0, minute=0, second=0)
    end_date = rupture_date.replace(hour=11, minute=59, second=59) # the day of the earthquake

    # Format the datetime object into a string
    start_date = start_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    end_date = end_date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # Define the ASF query parameters
    params = {
        'intersectsWith': AOI.wkt,
        'dataset': 'SENTINEL-1',
        'processingLevel': 'SLC',
        'start': start_date,
        'end': end_date,
        'output': 'geojson'
    }

    print('Performing ASF DAAC API query to return path and frame numbers for SLCs intersecting AOI over the last 24 days...')

    try:
        # Fetch data from the ASF DAAC API
        response = requests.get(ASF_DAAC_API, params=params)
        response.raise_for_status()  # Raise error if request fails

        # Parse the response as GeoJSON
        data = geojson.loads(response.text)

        # Initialize an empty dictionary to store the path and frame numbers as sets
        path_frame_numbers = defaultdict(lambda: defaultdict(set))

        # Extract the path and frame numbers from the GeoJSON data
        for feature in data['features']:
            flight_direction = feature['properties']['flightDirection']
            path_number = feature['properties']['pathNumber']
            frame_number = feature['properties']['frameNumber']
            path_frame_numbers[flight_direction][path_number].add(frame_number)  # Use a set to avoid duplicates

        # Reformat the dictionary into the desired format with tuples as keys and lists as values
        reformatted = {
            (flight_direction, path_number): sorted(frame_numbers)  # Convert to sorted list for consistency
            for flight_direction, path_frame in path_frame_numbers.items()
            for path_number, frame_numbers in path_frame.items()
        }

        print('=========================================')
        print('Reformatted Path and Frame Numbers for SLCs:')
        print('=========================================')
        for key, value in reformatted.items():
            print(f"{key}: {value}")
        return reformatted

    except requests.RequestException as e:
        print(f"Error accessing ASF DAAC API: {e}")
        return None

def get_SLCs(flight_direction, path_number, frame_numbers, time):
    """
    Query the ASF DAAC API for SLC data based on the given path and frame numbers.
    The data are organized by flight direction, path number, and frame numbers.
    :param: flight_direction - 'ASCENDING' or 'DESCENDING'
    :param: path_number - Sentinel-1 path number
    :param: frame_numbers - list of Sentinel-1 frame numbers
    :param: time - Unix timestamp representing the earthquake's origin time
    :return: List of dictionaries containing SLC fileIDs and their respective dates
    """
    # Establish the date range for the query
    rupture_date = convert_time(time)
    start_date = rupture_date - timedelta(days=90)  # 90 days before the earthquake
    start_date= start_date.replace(hour=0, minute=0, second=0)
    end_date = rupture_date + timedelta(days=30)    # 30 days after the earthquake
    end_date = end_date.replace(hour=11, minute=59, second=59)

    # Format the datetime object into a string
    start_date = start_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    end_date = end_date.strftime('%Y-%m-%dT%H:%M:%SZ')

    # Define the query parameters
    params = {
        'flightDirection': flight_direction,
        'frame': ','.join(frame_numbers),
        'relativeOrbit': path_number,
        'dataset':'SENTINEL-1',
        'processingLevel': 'SLC',
        'start':start_date,
        'end':end_date,
        'output': 'geojson'
    }

    print('Performing ASF DAAC API query to return SLCs for the given path and frame numbers...')
    
    try:
        # Fetch data from the ASF DAAC API
        response = requests.get(ASF_DAAC_API, params=params)
        response.raise_for_status()  # Raise error if request fails

        # Parse the response as GeoJSON
        data = geojson.loads(response.text)

        # Extract the file IDs from the GeoJSON data
        SLCs = []
        for feature in data['features']:

            # Extract the date part from startTime
            start_time = feature['properties']['startTime']
            path = feature['properties']['pathNumber']
            frame = feature['properties']['frameNumber']
            date = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S.%fZ').date().isoformat()
            
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
    
    except requests.RequestException as e:
        print(f"Error accessing ASF DAAC API: {e}")
        return None

def generate_pairs(pairs, mode):
    """
    Generate pairs of SLCs based on the selected pairing mode.
    :param: pairs - List of SLC pairs sorted by date
    :param: mode - 'sequential' for consecutive pairs, 'all' for all possible pairs
    :return: List of SLC pairs based on the mode
    """
    if mode == 'sequential':
        return [(pairs[i], pairs[i + 1]) for i in range(len(pairs) - 1)]
    elif mode == 'all':
        all_pairs = list(combinations(pairs, 2))
        return all_pairs
    elif mode == 'coseismic':
        return []  # Co-seismic handled separately

def find_reference_and_secondary_pairs(SLCs, time, flight_direction, path_number, title, mode='sequential'):
    """
    Find the reference and secondary pairs of SLCs necessary to run dockerized topsApp, 
    and determine whether each pair is pre-seismic, co-seismic, or post-seismic based on the rupture date and SLC dates.
    :param: SLCs - List of dictionaries containing SLC fileIDs and their respective dates
    :param: time - Unix timestamp representing the earthquake's origin time
    :param: flight_direction - 'ASCENDING' or 'DESCENDING'
    :param: path_number - Sentinel-1 path number
    :param: title: - USGS title of the earthquake event, used for file organization
    :param: mode - 'all', 'sequential', or 'coseismic' for pairing mode
    :return: List of JSON objects containing the parameters for each pair of SLCs
    """
    # Get the rupture date in the format YYYY-MM-DD
    rupture_date = convert_time(time)
    rupture_date = rupture_date.strftime('%Y-%m-%d')
    rupture_date_dt = datetime.strptime(rupture_date, '%Y-%m-%d')
    
    # Reformatting for dictionary keys for later use    
    flight_direction = 'A' if flight_direction == 'ASCENDING' else 'D'
    path_number = f"{int(path_number):03}" 

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
    
    pre_seismic = [pair for pair in initial_pairs if pair[0] < rupture_date_dt]
    post_seismic = [pair for pair in initial_pairs if pair[0] > rupture_date_dt]
    co_seismic = []
    
    for i in range(len(initial_pairs) - 1):
        if initial_pairs[i][0] < rupture_date_dt < initial_pairs[i + 1][0]:
            co_seismic = [(initial_pairs[i], initial_pairs[i + 1])]
            break
    
    # Handle different modes
    if mode == 'coseismic':
        paired_results = {'co-seismic': co_seismic}
    else:
        paired_results = {
            'pre_seismic': generate_pairs(pre_seismic, mode),
            'co_seismic': co_seismic,
            'post_seismic': generate_pairs(post_seismic, mode)
        }
    
    isce_jsons = []
    for timing, pairs in paired_results.items():
        for (secondary, reference) in pairs:
            reference_date, reference_scenes = reference
            secondary_date, secondary_scenes = secondary
            reference_scenes_ids = [slc['fileID'] for slc in reference_scenes]
            secondary_scenes_ids = [slc['fileID'] for slc in secondary_scenes]
            json_output = make_json(title, timing, flight_direction, path_number, list(frame_numbers), 
                                    {'date': reference_date.strftime('%Y-%m-%d')}, 
                                    {'date': secondary_date.strftime('%Y-%m-%d')}, 
                                    reference_scenes_ids, secondary_scenes_ids)
            isce_jsons.append(json_output)
    
    return isce_jsons

def make_json(title, timing, flight_direction, path_number, frame_numbers, reference, secondary, reference_scenes, secondary_scenes):
    """Create a JSON object containing parameters for dockerized topsApp.
    Note: Not all params here are used in the final dockerized topsApp. Some are used for file organzation.
    Note: Several params are 'hardcoded', as these should not be modified for individual products.
    :param: title - USGS title of the earthquake event
    :param: timing - 'pre-seismic', 'co-seismic', or 'post-seismic'
    :param: flight_direction - 'A' or 'D' for ascending or descending
    :param: path_number - Sentinel-1 path number
    :param: frame_numbers - list of intersecting Sentinel-1 frame numbers
    :param: reference - dictionary containing the reference date
    :param: secondary - dictionary containing the secondary date
    :param: reference_scenes - list of reference SLC fileIDs
    :param: secondary_scenes - list of secondary SLC fileIDs
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

def create_directories_from_json(eq_jsons, root_dir):
    """
    Create directories for each group of SLCs based on the JSON data provided.
    The directories are created in the root directory specified and will have a name following this format: 'flight_directionpath_number_secondary_date_reference_date'.
    :param: eq_jsons - list of JSON objects containing the parameters for each pair of SLCs
    :param: root_dir - root directory where the directories will be created
    :return: List of directory names
    """
    dirnames = []
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
        
        dirnames.append(sub_dirnames)

    return dirnames

def run_dockerized_topsApp(json_data):
    """
    Run dockerized topsApp InSAR processing workflow using the provided JSON data.
    Outputs are added to the root dir + an extension for each pair.
    :param: json_data - JSON object containing the parameters for dockerized topsApp
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
    
    # Construct the command
    command = [
        "taskset", "-c", "65-128",
        "nohup", "isce2_topsapp",
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
    ]
    
    print('=========================================')
    print("Running the dockerized topsApp InSAR processing workflow...")
    print('=========================================')

    # Set the environment variable for XLA
    env = os.environ.copy()
    env["XLA_PYTHON_CLIENT_MEM_FRACTION"] = "0.20"

    # Write nohup.out in the current working directory
    nohup_out_file = os.path.join(os.getcwd(), "nohup.out")
    print(f"Writing nohup output to: {nohup_out_file}")

    # Run the command
    with open(nohup_out_file, "w") as nohup_out:
        try:
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True,
                env=env  # Pass the environment with XLA variable set
            )
            print("Command executed successfully!")
            print("Output:\n", result.stdout)
        except subprocess.CalledProcessError as e:
            print("Error occurred while running the command.")
            print("Error Output:\n", e.stderr)
    return

def to_snake_case(input_string):
    """
    Convert the given string to snake_case.
    :param: input_string - string to convert to snake_case
    :return: snake_case version of the input string
    """
    # Replace non-alphanumeric characters with spaces
    cleaned_string = re.sub(r'[^\w\s]', '', input_string)
    # Replace spaces with underscores and convert to lowercase
    snake_case_string = re.sub(r'\s+', '_', cleaned_string.strip()).lower()
    return snake_case_string

def main_forward():
    """
    Runs the main query and processing workflow for the forward processing mode.
    Produces co-seismic displacement product for new earthquakes.
    """
    # Fetch GeoJSON data from the USGS Earthquake Hazard Portal
    #geojson_data = check_for_new_data(USGS_api_hourly)
    geojson_data = check_for_new_data(USGS_api_30day)
    
    if geojson_data:
        # Parse GeoJSON and create variables for each feature's properties
        earthquakes = parse_geojson(geojson_data)
        eq_sig = check_significance(earthquakes)

        if eq_sig is not None:
            for eq in eq_sig:
                coords = eq.get('coordinates', [])
                aoi = make_aoi(coords)
            
                # Query the ASF DAAC API for SAR data within the AOI
                SLCs = query_asfDAAC(aoi, eq.get('time'))

            # Find SLC pairs for InSAR processing
            reference, secondary = find_SLC_pairs_for_infg(SLCs, aoi, eq.get('time'))

        else:
            print("No significant earthquakes found in the last hour.")
    
def main_historic(start_date, end_date = None, pairing_mode = None):
    """
    Runs the main query and processing workflow for the historic processing mode.
    Produces pre-seismic, co-seismic, and post-seismic displacement products for historic earthquakes.
    Pre-seismic and post-seismic data are generated for 90 days before and 30 days after the event.
    :param: start_date - the start date in YYYY-MM-DD format
    :param: end_date - the optional end date in YYYY-MM-DD format
    :param: pairing_mode - 'all', 'sequential', or 'coseismic' for SLC pairing
    """
    if str(start_date) and str(end_date) is None:
        # Fetch GeoJSON data from the USGS Earthquake Hazard Portal for a single date
        geojson_data = get_historic_earthquake_data_single_date(USGS_api_alltime, start_date)

    elif str(start_date) and str(end_date):
        # Fetch GeoJSON data from the USGS Earthquake Hazard Portal for the date range provided
        geojson_data = get_historic_earthquake_data_date_range(USGS_api_alltime, start_date, end_date)

    if geojson_data:
        # Parse GeoJSON and create variables for each feature's properties
        earthquakes = parse_geojson(geojson_data)
        eq_sig = check_significance(earthquakes)

        if eq_sig is not None:
            for eq in eq_sig:
                title = eq.get('title', '')
                title = to_snake_case(title)
                print(f"title: {title}")
                coords = eq.get('coordinates', [])
                aoi = make_aoi(coords)
                
                path_frame_numbers = get_path_and_frame_numbers(aoi, eq.get('time'))
                print('path_frame_numbers:', path_frame_numbers)

                eq_jsons = []
                for (flight_direction, path_number), frame_numbers in path_frame_numbers.items():
                    SLCs = get_SLCs(flight_direction, path_number, frame_numbers, eq.get('time'))
                    isce_jsons = find_reference_and_secondary_pairs(SLCs, eq.get('time'), flight_direction, path_number, title, pairing_mode)
                    eq_jsons.append(isce_jsons)

                dirnames = create_directories_from_json(eq_jsons, root_dir)
                
                for i, eq_json in enumerate(eq_jsons):
                    for j, json_data in enumerate(eq_json): 
                        working_dir = dirnames[i][j]
                        os.chdir(working_dir)
                        # Add json to the directory 
                        with open('isce_params.json', 'w') as f:
                            json.dump(json_data, f, indent=4)
                        print(f'working directory: {os.getcwd()}')
                        # print(f'Running dockerized topsApp for dates {json_data["secondary-date"]} to {json_data["reference-date"]}...')
                        # try:
                        #     run_dockerized_topsApp(json_data)
                        # except:
                        #     print('Error running dockerized topsApp')
                        #     continue

            return eq_jsons
                    
        else:
            print(f"No significant earthquakes found betweeen {start_date} and {end_date}.")

if __name__ == "__main__":
    """
    Run the main function based on the input arguments provided, in either historic or forward processing mode.
    Example usage for historic processing: python coseis.py --historic --start_date 2021-08-14 --end_date 2021-08-15 --pairing all
    Example usage for forward processing: python coseis.py --forward
    """
    parser = argparse.ArgumentParser(
        description="Run historic or forward processing based on input arguments."
    )

    parser.add_argument("--historic", action="store_true", help="Run historic processing.")
    parser.add_argument("--forward", action="store_true", help="Run forward processing.")
    parser.add_argument("--start_date", help="The start date in YYYY-MM-DD format (required for historic processing).")
    parser.add_argument("--end_date", help="The optional end date in YYYY-MM-DD format.")
    parser.add_argument("--pairing", choices=["all", "sequential","coseismic"], help="Specify the SLC pairing mode. Required for historic processing."
    )

    args = parser.parse_args()

    # Ensure only one mode is selected
    if args.historic and args.forward:
        print("Error: You cannot specify both --historic and --forward.")
        parser.print_help()
        exit(1)

    if args.historic:
        if not args.start_date:
            print("Error: --start_date is required when using --historic mode.")
            parser.print_help()
            exit(1)

        if not args.pairing:
            print("Error: --pairing is required when using --historic mode. Options: 'all', 'sequential', 'coseismic.")
            parser.print_help()
            exit(1)

        main_historic(args.start_date, args.end_date, pairing_mode=args.pairing)

    elif args.forward:
        if args.start_date or args.end_date or args.pairing:
            print("Error: --start_date, --end_date, and --pairing cannot be used with --forward mode.")
            parser.print_help()
            exit(1)

        main_forward()

    else:
        print("Error: Please specify either --historic or --forward to run the processing workflow.")
        parser.print_help()
        exit(1)