import requests
import geojson
from shapely.geometry import Point, Polygon
from datetime import datetime, timedelta, timezone


# Configuration
#USGS_api = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_hour.geojson"  # USGS Earthquake API - Hourly
USGS_api = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_month.geojson"  # USGS Earthquake API - Monthly
coastline_api = "https://raw.githubusercontent.com/OSGeo/PROJ/refs/heads/master/docs/plot/data/coastline.geojson" # Coastline API
ASF_DAAC_API = "https://api.daac.asf.alaska.edu/services/search/param"

def check_for_new_data(eq_api):
    """
    Fetches data from the USGS Earthquake Portal and returns it as a GeoJSON object.
    The data returned is updated hourly and includes all earthquakes occuring during that period.
    """
    try:
        # Fetch data from the USGS Earthquake API
        response = requests.get(eq_api)
        response.raise_for_status()  # Raise error if request fails
        
        # Parse the response as GeoJSON
        data = geojson.loads(response.text)
        return data

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
        data = geojson.loads(response.text)
        return data
    
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

def check_significance(earthquakes):
    """
    Checks the significance of each earthquake based on its 
    (1) magnitude, (2) USGS alert level, (3) depth, and (4) distance from land.
    Criteria: Magnitude >= 6.0, USGS alert level = ['yellow','orange','red], and Depth <= 30.0 km
    """
    significant_earthquakes = []
    alert_list = ['yellow', 'orange', 'red']

    for earthquake in earthquakes:
        magnitude = earthquake.get('mag')
        alert = earthquake.get('alert')
        depth = earthquake.get('coordinates', [])[2] if earthquake.get('coordinates') else None
        if (magnitude >= 6.0) and (alert in alert_list) and (depth <= 30.0):
            significant_earthquakes.append(earthquake)

    return significant_earthquakes

def make_aoi(coordinates):
    """
    Create an Area of Interest (AOI) polygon based on the given coordinates.
    """
    # Extract the X and Y coordinates of the earthquake's epicenter
    X = coordinates[0]
    Y = coordinates[1]

    # Define the side length of the square AOI in decimal degrees
    side_length = 0.5  # 1 degree is ~111 km at the equator

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
    
def query_ASFDAAC(AOI, time):
    """
    Query the ASF DAAC API for SAR data within the Area of Interest (AOI).
    """
    # Define the ASF DAAC API endpoint
    ASF_DAAC_API = "https://api.daac.asf.alaska.edu/services/search/param"

    # Establish the date range for the query
    rupture_date = convert_time(time)
    start_date = rupture_date - timedelta(days=12)  # 12 days before the earthquake
    start_date= start_date.replace(hour=0, minute=0, second=0)
    end_date = rupture_date + timedelta(days=12)    # 12 days after the earthquake
    end_date = end_date.replace(hour=11, minute=59, second=59)

    # Format the datetime object into a string
    start_date = start_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    end_date = end_date.strftime('%Y-%m-%dT%H:%M:%SZ')
    
    print(f"Start Date: {start_date}")
    print(f"End Date: {end_date}")
    
    # Define the query parameters
    params = {
        'intersectsWith': AOI.wkt,
        'dataset': 'SENTINEL-1',
        'start':start_date,
        'end':end_date,
        'output': 'geojson'
    }

    try:
        # Fetch data from the ASF DAAC API
        response = requests.get(ASF_DAAC_API, params=params)
        response.raise_for_status()  # Raise error if request fails
        
        # Parse the response as GeoJSON
        data = geojson.loads(response.text)

        # Open file in write mode and dump the GeoJSON data into it
        with open('ASF_query.geojson', 'w') as f:
            geojson.dump(data, f, indent=2)

        return data
    
    except requests.RequestException as e:
        print(f"Error accessing ASF DAAC API: {e}")
        return None

if __name__ == "__main__":
    # Fetch GeoJSON data from the USGS Earthquake Hazard Portal
    geojson_data = check_for_new_data(USGS_api)
    
    if geojson_data:
        # Parse GeoJSON and create variables for each feature's properties
        earthquakes = parse_geojson(geojson_data)
        eq_sig = check_significance(earthquakes)

        if eq_sig:
            print("Significant Earthquakes:")
            for eq in eq_sig:
                print(eq)
        else:
            print("No significant earthquakes found.")

        for eq in eq_sig:
            coords = eq.get('coordinates', [])
            aoi = make_aoi(coords)
            print(f"Area of Interest (AOI) for Earthquake: {aoi}")
            
            # Query the ASF DAAC API for SAR data within the AOI
            ASF_geojson = query_ASFDAAC(aoi, eq.get('time'))

            # Parse the ASF GeoJSON data
            ASF_data = parse_geojson(ASF_geojson)
    # # Fetch GeoJSON data from the coastline API
    # coastline_data = get_coastline(coastline_api)


