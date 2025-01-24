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
    Fetches data from the USGS Earthquake Portal for a single date and returns it as a GeoJSON object.
    The data returned will depend on the parameters included with the API request.
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
    Fetches data from the USGS Earthquake Portal over the date range and returns it as a GeoJSON object.
    The data returned will depend on the parameters included with the API request.
    """
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
    Fetches data from the USGS Earthquake Portal and returns it as a GeoJSON object.
    The data returned is updated hourly and includes all earthquakes occuring during that period.
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
    Fetches the coastline data from the specified API and returns it as a GeoJSON object.
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
    Parses the features of a GeoJSON object and creates a dictionary for each earthquake (feature),
    with property names as the keys and property values as the values.

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
    """
    # Extract the coordinates of the earthquake's epicenter
    coords = earthquake.get('coordinates', [])
    if not coords:
        return None

    # Create a Point object from the earthquake's coordinates
    epicenter = Point(coords[:2])
    
    # Determine if the epicenter is within the coastline
    within_coastline = coastline.contains(epicenter)
    
    return within_coastline

def check_significance(earthquakes):
    """
    Checks the significance of each earthquake based on its 
    (1) magnitude, (2) USGS alert level, (3) depth, and (4) distance from land.
    Criteria: Magnitude >= 6.0, USGS alert level = ['yellow','orange','red], and Depth <= 30.0 km
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
        within_Coastline = withinCoastline(earthquake, coastline)
        if all(var is not None for var in (magnitude, alert, depth)):
            if (magnitude >= 6.0) and (alert in alert_list) and (depth <= 30.0):
                significant_earthquakes.append(earthquake)

    # Write significant earthquakes to a GeoJSON file
    if len(significant_earthquakes) > 0:
        print('=========================================')
        print(f"Found {len(significant_earthquakes)} significant earthquakes.")
        print('=========================================')
        print(significant_earthquakes)
        print('=========================================')
        significant_earthquakes_to_geojson(significant_earthquakes)
        return significant_earthquakes
    else:
        return None

def significant_earthquakes_to_geojson(significant_earthquakes):
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
        
def make_aoi(coordinates):
    """
    Create an Area of Interest (AOI) polygon based on the given coordinates.
    """
    print('=========================================')
    print("Creating Area of Interest (AOI) polygon...")
    print('=========================================')

    # Extract the X and Y coordinates of the earthquake's epicenter
    X = coordinates[0]
    Y = coordinates[1]

    # Define the side length of the square AOI in decimal degrees
    #side_length = 0.5  # 1 degree is ~111 km at the equator
    side_length = 1.0

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
    Convert the given time string to a Unix timestamp.
    """
    # Convert milliseconds to seconds
    timestamp_s = time / 1000

    # Convert to datetime object in UTC timezone
    dt = datetime.fromtimestamp(timestamp_s, tz=timezone.utc)
    dt = dt.replace(microsecond=0)
    return dt

def get_path_and_frame_numbers(AOI, time):
    """
    Query the ASF DAAC API for S1 data intersecting the Area of Interest (AOI) over the previous 24 days.
    Return the *unique* path and frame numbers for each intersecting SLC.
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

    # Establish the date range for the query
    rupture_date = convert_time(time)
    start_date = rupture_date - timedelta(days=90)  # 3 months before the earthquake
    start_date= start_date.replace(hour=0, minute=0, second=0)
    end_date = rupture_date + timedelta(days=90)    # 3 months before the earthquake
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
            date = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%S.%fZ').date().isoformat()
            SLC = {
                'fileID': feature['properties']['fileID'],
                'date': date
            }

            SLCs.append(SLC)

        # Print the SLCs
        print('=========================================')
        print(f"Found {len(SLCs)} SLCs for the {flight_direction} path {path_number} and frame numbers {frame_numbers}.")
        # for SLC in SLCs:
        #     print(SLC)
        return SLCs
    
    except requests.RequestException as e:
        print(f"Error accessing ASF DAAC API: {e}")
        return None

def find_reference_and_secondary_pairs(SLCs, flight_direction, path_number, time, title):
    """
    Find the reference and secondary pairs of SLCs and determine whether each pair is 
    pre-seismic, co-seismic, or post-seismic based on the rupture date and SLC dates."""

    # Get the rupture date in the format YYYY-MM-DD
    rupture_date = convert_time(time)
    rupture_date = rupture_date.strftime('%Y-%m-%d')

    # Reformatting for dictionary keys for later use    
    flight_direction = 'A' if flight_direction == 'ASCENDING' else 'D'
    path_number = f"{int(path_number):03}" 
    
    # Sort SLCs by date in ascending order
    SLCs_sorted = sorted(SLCs, key=lambda x: x['date'])
    
    # Initialize the result dictionary
    #pairs_dict = {}
    
    isce_jsons = []

    # Create pairs
    for i in range(1, len(SLCs_sorted)):
        reference = SLCs_sorted[i]  # Newer
        secondary = SLCs_sorted[i - 1]  # Older
        
        # Skip if the dates are the same
        if reference['date'] == secondary['date']:
            continue
        
        # Determine the timing of pair based on the rupture_date
        if reference['date'] < rupture_date:
            timing = 'pre-seismic'
        elif reference['date'] > rupture_date and secondary['date'] < rupture_date:
            timing = 'co-seismic'
        else:
            timing = 'post-seismic'

        # Create unique key for each pair
        key = (flight_direction,
            path_number,
            secondary['date'],
            reference['date'],
            timing
            )
        
        # Collect fileIDs for reference/secondary
        reference_scenes = [
            entry['fileID'] for entry in SLCs if entry['date'] == reference['date']
        ]
        secondary_scenes = [
            entry['fileID'] for entry in SLCs if entry['date'] == secondary['date']
        ]
        
        # Store in the nested dictionary
        # pairs_dict[key] = {
        #     ('reference-scenes', 'secondary-scenes'): [reference_scenes, secondary_scenes]
        # }
    
        # Create the JSON for this pair
        json_output = make_json(title, flight_direction, path_number, {'date': reference['date']}, {'date': secondary['date']}, reference_scenes, secondary_scenes)
        
        # Append the JSON to the list
        isce_jsons.append(json_output)

    # Print the pairs in a readable format
    # print("\nReadable Output:\n")
    # pp = pprint.PrettyPrinter(indent=4)
    # pp.pprint(pairs_dict)
    
    return isce_jsons

def make_json(title, flight_direction, path_number, reference, secondary, reference_scenes, secondary_scenes):
    """Create a JSON object containing parameters to make a directory ."""
    isce_json = {
        "title": title,
        "flight-direction": flight_direction,
        "path-number": path_number,
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
    dirnames = []
    for isce_jsons in eq_jsons:
        for json_data in isce_jsons:
            title = json_data['title']
            flight_direction = json_data['flight-direction']
            path_number = json_data['path-number']
            secondary_date = json_data['secondary-date'].replace("-", "")
            reference_date = json_data['reference-date'].replace("-", "")
            
            # Build the full path
            base_path = os.path.join(root_dir, title, flight_direction + path_number)
            sub_path = f"{flight_direction}{path_number}_{secondary_date}_{reference_date}"
            full_path = os.path.join(base_path, sub_path)
            
            # Create directories, ensuring no overwriting
            os.makedirs(full_path, exist_ok=True)
            dirnames.append(full_path)
            print(f"Created: {full_path}")

    return dirnames

def run_dockerized_topsApp(json_data):
    """
    Run dockerized topsApp InSAR processing workflow using the provided JSON data.
    Outputs are added to the root dir + an extension for each pair.
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

    # Run the command
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

def query_asfDAAC(time):
    """
    Query the ASF DAAC API for S1 data within the Area of Interest (AOI).
    """

    # Establish the date range for the query
    rupture_date = convert_time(time)
    start_date = rupture_date - timedelta(days=12)  # 12 days before the earthquake
    start_date= start_date.replace(hour=0, minute=0, second=0)
    end_date = rupture_date + timedelta(days=12)    # 12 days after the earthquake
    end_date = end_date.replace(hour=11, minute=59, second=59)

    # Format the datetime object into a string
    start_date = start_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    end_date = end_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    
    # Define the query parameters
    params = {
        'intersectsWith': AOI.wkt,
        'dataset': 'SENTINEL-1',
        'processingLevel': 'SLC',
        'start':start_date,
        'end':end_date,
        'output': 'geojson'
    }

    print('Performing ASF DAAC API query...')
    print("Query Parameters:")
    print(json.dumps(params, indent=4))

    try:
        # Fetch data from the ASF DAAC API
        response = requests.get(ASF_DAAC_API, params=params)
        response.raise_for_status()  # Raise error if request fails

        # Parse the response as GeoJSON
        data = geojson.loads(response.text)


        # Filter the features based on the "flightDirection" in the "properties" field
        ascending_data = {'type': 'FeatureCollection', 'features': [feature for feature in data['features'] if feature['properties'].get('flightDirection') == 'ASCENDING']}
        descending_data = {'type': 'FeatureCollection', 'features': [feature for feature in data['features'] if feature['properties'].get('flightDirection') == 'DESCENDING']}
        
        # Write the filtered data to separate GeoJSON files
        with open('ASF_query_ascending.geojson', 'w') as f_ascending:
            geojson.dump(ascending_data, f_ascending, indent=2)

        with open('ASF_query_descending.geojson', 'w') as f_descending:
            geojson.dump(descending_data, f_descending, indent=2)

        # Extract the file IDs from the GeoJSON data
        result = []
        for feature in data['features']:
            SLC = {
                'fileID': feature['properties']['fileID'],
                'flightDirection': feature['properties']['flightDirection'],
                'pathNumber': feature['properties']['pathNumber'],
                'frameNumber': feature['properties']['frameNumber'],
                'startTime': feature['properties']['startTime'],
                'geometry': feature.geometry
            }
            result.append(SLC)
        print('=========================================')
        print(f"Found {len(result)} SLCs intersecting the AOI.")
        for SLC in result:
            print(f"SLC ID: {SLC['fileID']}")
        print('=========================================')
        return result
    
    except requests.RequestException as e:
        print(f"Error accessing ASF DAAC API: {e}")
        return None

def restructure_groups(ascending_groups, descending_groups):
    """
    Restructure (1) style groups into (2) style groups to make them compatible 
    with the original function.
    """
    # Flatten and merge the groups in ascending_groups
    if ascending_groups:
        ascending_group = []
        for group in ascending_groups:
            for item in group:
                ascending_group.append(item)
    else:
        ascending_group = []

    # Flatten and merge the groups in descending_groups
    if descending_groups:
        descending_group = []
        for group in descending_groups:
            for item in group:
                descending_group.append(item)
    else:
        descending_group = []
    
    return ascending_group, descending_group

def find_all_groups(groups, rupture_datetime):
    """
    Find all possible groups of SLCs for InSAR processing based on the flight direction and path number.
    """
    print('=========================================')
    print('Finding all possible ASCENDING AND DESCENDING tracks')
    print('=========================================')

    ascending_groups = []
    descending_groups = []

    for key, slcs in groups.items():
        flight_direction, _ = key  # Extract flight direction from the key

        if flight_direction == "A":
            ascending = (key, slcs)
            ascending_groups.append(ascending)

        if flight_direction == "D":
            descending = (key, slcs)
            descending_groups.append(descending)

    all_groups = [ascending_groups, descending_groups]

    print('=========================================')
    print('Make REFERENCE and SECONDARY LISTS...')
    print('=========================================')
    
    # Initialize dictionaries to store the "REFERENCE" and "SECONDARY" dictionaries
    split_groups = {}

    # Iterate through the best groups (best ascending and best descending)
    ascending_groups = []
    descending_groups = []
    for group in all_groups:
        for (flight_direction, path_number), slcs in group:
            # Sort SLCs by startTime to ensure correct order
            slcs.sort(key=lambda x: x['startTime'])
            
            # Initialize empty dictionaries for REFERENCE and SECONDARY
            reference_slcs = defaultdict(list)
            secondary_slcs = defaultdict(list)

            # Split the SLCs based on rupture_datetime
            for slc in slcs:
                slc['startTime'] = datetime.strptime(slc['startTime'], '%Y-%m-%dT%H:%M:%S.%fZ')
            
                # Extract the date part from the startTime
                slc_date = slc['startTime'].date()  # Get only the date (YYYY-MM-DD)
        
                if slc['startTime'] >= rupture_datetime:
                    reference_slcs[slc_date].append(slc['fileID'].split('-SLC')[0])
                else:
                    secondary_slcs[slc_date].append(slc['fileID'].split('-SLC')[0])
            
            # Reverse the secondary dictionary by flipping the order of the dates
            reversed_secondary = dict(reversed(list(secondary_slcs.items())))

            # Store the split SLCs in the dictionary
            split_groups[(flight_direction, path_number)] = {
                "REFERENCE": reference_slcs,
                "SECONDARY": reversed_secondary
            }

            # Store the REFERENCE and SECONDARY dictionaries in list based on flight direction
            if flight_direction == "A":
                ascending_group = [reference_slcs, reversed_secondary]
                ascending_groups.append(ascending_group)

            elif flight_direction == "D":
                descending_group = [reference_slcs, reversed_secondary]
                descending_groups.append(descending_group)

        print('=========================================')
        print('REFERENCE and SECONDARY SLCs:')
        print('=========================================')

    print(f'ascending_groups: {ascending_groups}')
    print(f'descending_groups: {descending_groups}')

    ascending_group, descending_group = restructure_groups(ascending_groups, descending_groups)
    print(f'ascending_group: {ascending_group}')
    print(f'descending_group: {descending_group}')
    
    return ascending_group, descending_group

def find_optimal_groups(AOI, groups, rupture_datetime):
    """
    Find the optimal groups of SLCs for InSAR processing based on the overlap with the Area of Interest (AOI).
    """
    print('=========================================')
    print('Finding optimal ASCENDING and DESCENDING tracks')
    print('=========================================')

    overlap_areas = {}  # To store the total overlap area for each group
    for key, slcs in groups.items():
        # Compute the union of all SLC geometries in the group
        slc_geometries = [shape(slc['geometry']) for slc in slcs]
        slc_union = unary_union(slc_geometries)
        # Compute overlap with the AOI
        total_overlap_area = 0
        overlap_area = slc_union.intersection(AOI).area
        total_overlap_area += overlap_area
        # Store the total overlap area for the group
        overlap_areas[key] = total_overlap_area
        # Print the group and its overlap area
        print(f"{key} with total overlap of the AOI of {total_overlap_area}")

    # Initialize variables to track the groups with the highest overlap
    best_ascending = None
    best_descending = None
    max_overlap_ascending = 0
    max_overlap_descending = 0

    # Iterate through groups and find the ones with the highest overlap
    for key, slcs in groups.items():
        flight_direction, _ = key  # Extract flight direction from the key
        
        # Get the total overlap area from the overlap_areas dictionary
        total_overlap_area = overlap_areas.get(key, 0)

        if flight_direction == "ASCENDING" and total_overlap_area > max_overlap_ascending:
            max_overlap_ascending = total_overlap_area
            best_ascending = (key, slcs)

        if flight_direction == "DESCENDING" and total_overlap_area > max_overlap_descending:
            max_overlap_descending = total_overlap_area
            best_descending = (key, slcs)

    best_groups = [best_ascending, best_descending]

    # print("Best ASCENDING Group:", best_ascending[0], "with total overlap area", max_overlap_ascending)
    # print("Best DESCENDING Group:", best_descending[0], "with total overlap area", max_overlap_descending)

    # Initialize dictionaries to store the "REFERENCE" and "SECONDARY" dictionaries
    split_groups = {}

    # Iterate through the best groups (best ascending and best descending)
    for (flight_direction, path_number), slcs in best_groups:
        # Sort SLCs by startTime to ensure correct order
        slcs.sort(key=lambda x: x['startTime'])
        
        # Initialize empty dictionaries for REFERENCE and SECONDARY
        reference_slcs = defaultdict(list)
        secondary_slcs = defaultdict(list)

        # Split the SLCs based on rupture_datetime
        for slc in slcs:
            slc['startTime'] = datetime.strptime(slc['startTime'], '%Y-%m-%dT%H:%M:%S.%fZ')
        
            # Extract the date part from the startTime
            slc_date = slc['startTime'].date()  # Get only the date (YYYY-MM-DD)
    
            if slc['startTime'] >= rupture_datetime:
                reference_slcs[slc_date].append(slc['fileID'].split('-SLC')[0])
            else:
                secondary_slcs[slc_date].append(slc['fileID'].split('-SLC')[0])
        
        # Reverse the secondary dictionary by flipping the order of the dates
        reversed_secondary = dict(reversed(list(secondary_slcs.items())))

        # Store the split SLCs in the dictionary
        split_groups[(flight_direction, path_number)] = {
            "REFERENCE": reference_slcs,
            "SECONDARY": reversed_secondary
        }

        # Store the REFERENCE and SECONDARY dictionaries in list based on flight direction
        if flight_direction == "ASCENDING":
            ascending_group = [reference_slcs, reversed_secondary]
        elif flight_direction == "DESCENDING":
            descending_group = [reference_slcs, reversed_secondary]

    print('=========================================')
    print('REFERENCE and SECONDARY SLCs:')
    print('=========================================')
    # Print the result
    # for (flight_direction, path_number), splits in split_groups.items():
    #     print(f"Flight Direction: {flight_direction}, Path Number: {path_number}")
    #     print(f"REFERENCE SLCs: {dict(splits['REFERENCE'])}")  # Convert defaultdict to dict for display
    #     print(f"SECONDARY SLCs: {dict(splits['SECONDARY'])}")  # Convert defaultdict to dict for display
    
    print(f'ascending_group: {ascending_group}')
    print(f'descending_group: {descending_group}')

    return ascending_group, descending_group

def find_SLC_pairs_for_infg(SLCs, AOI, event_datetime, optimal = False):
    """
    Find the SLCs that cover the Area of Interest (AOI) with appropriate acquisition times bounding the event to produce the InSAR product.
    """

    print('=========================================')
    print('Finding SLC pairs for InSAR processing...')
    print('=========================================')

    rupture_datetime = convert_time(event_datetime).strftime('%Y-%m-%dT%H:%M:%SZ')
    rupture_datetime = datetime.strptime(rupture_datetime, '%Y-%m-%dT%H:%M:%SZ')

    # Find the SLCs that intersect the AOI with the most overlap
    groups = defaultdict(list)

    for SLC in SLCs:
        # Extract flightDirection and pathNumber
        flight_direction = SLC.get('flightDirection', '').upper()
        path_number = SLC.get('pathNumber')

        # Ensure both keys are present
        if (flight_direction is not None) and (path_number is not None):
            formatted_path_number = f"{int(path_number):03d}"
            if flight_direction == 'ASCENDING':
                direction = 'A'
            elif flight_direction == 'DESCENDING':
                direction = 'D'
            groups[(direction, formatted_path_number)].append(SLC)

    print('=========================================')
    print('Groups of SLCs:', groups)
    print('=========================================')

    if optimal:
        ascending_group, descending_group = find_optimal_groups(AOI, groups, rupture_datetime)

    else:
        ascending_group, descending_group = find_all_groups(groups, rupture_datetime)

    return ascending_group, descending_group

def to_snake_case(input_string):
    # Replace non-alphanumeric characters with spaces
    cleaned_string = re.sub(r'[^\w\s]', '', input_string)
    # Replace spaces with underscores and convert to lowercase
    snake_case_string = re.sub(r'\s+', '_', cleaned_string.strip()).lower()
    return snake_case_string

def main_forward():

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
    
def main_historic(start_date, end_date = None):
    """
    Query the USGS Earthquake API for significant earthquakes between the given date (and end date, if provided)
    and return json files for ascending and descending groups of reference/secondary SLCs for InSAR processing.
    """

    if start_date and end_date is None:
        # Fetch GeoJSON data from the USGS Earthquake Hazard Portal for a single date
        geojson_data = get_historic_earthquake_data_single_date(USGS_api_alltime, start_date)

    elif start_date and end_date:
        # Fetch GeoJSON data from the USGS Earthquake Hazard Portal for the date range provided
        geojson_data = get_historic_earthquake_data_date_range(USGS_api_alltime, start_date, end_date)

    else:
        print("Must provide (at least) a start date for historic processing. Exiting...")
        return

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
                    isce_jsons = find_reference_and_secondary_pairs(SLCs, flight_direction, path_number, eq.get('time'), title)
                    eq_jsons.append(isce_jsons)

                dirnames = create_directories_from_json(eq_jsons, root_dir)

                for eq_json in eq_jsons:
                    for i, json_data in enumerate(eq_json):
                        os.chdir(dirnames[i])
                        print(f'working directory: {os.getcwd()}')
                        print(f'Running dockerized topsApp for dates {json_data["secondary-date"]} to {json_data["reference-date"]}...')
                        try:
                            run_dockerized_topsApp(json_data)
                        except:
                            print('Error running dockerized topsApp')
                            continue
            return eq_jsons
                    
        else:
            print(f"No significant earthquakes found betweeen {start_date} and {end_date}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run historic or forward processing based on input arguments."
    )
    parser.add_argument("start_date", nargs="?", help="The start date in YYYY-MM-DD format.")
    parser.add_argument("end_date", nargs="?", help="The optional end date in YYYY-MM-DD format.")

    args = parser.parse_args()

    if args.start_date:
        # Call main_historic with one or two arguments based on what's provided
        main_historic(args.start_date, args.end_date)
    else:
        main_forward()