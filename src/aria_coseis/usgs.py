"""USGS earthquake catalog access: event queries, finite-fault models, rake, and custom event lists."""

import json
import time
from datetime import datetime, timezone

import geojson
import requests


def get_historic_earthquake_data_single_date(eq_api, input_date):
    """
    Fetch data from the USGS Earthquake Portal for a single date and returns it as a GeoJSON object.
    The data returned will depend on the parameters included with the API request.
    :param eq_api: USGS API endpoint
    :param input_date: date in the format 'YYYY-MM-DD'
    :return: GeoJSON object containing earthquake data
    """
    print("=========================================")
    print(f"Fetching historic earthquake data from {input_date}...")
    print("=========================================")
    try:
        # Parameters for the API request
        params = {
            "format": "geojson",
            "starttime": input_date + "00:00:00",
            "endtime": input_date + "23:59:59",
            "minmagnitude": 6.0,
            "maxdepth": 40.0,
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

    print("=========================================")
    print(f"Fetching historic earthquake data from {start_date} to {end_date}...")
    print("=========================================")
    try:
        # Parameters for the API request
        params = {
            "format": "geojson",
            "starttime": start_date,
            "endtime": end_date,
            "minmagnitude": 6.0,
            "maxdepth": 40.0,
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
    detail_url = "https://earthquake.usgs.gov/fdsnws/event/1/query"
    params = {"eventid": event_id, "format": "geojson"}

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


def parse_geojson(geojson_data):
    """
    Parse the features of a GeoJSON object and create a dictionary for each earthquake (feature),
    with property names as the keys and property values as the values.
    :param geojson_data: GeoJSON object containing earthquake data
    :return: List of dictionaries containing earthquake data
    """
    earthquakes = []

    # Loop through each feature in the GeoJSON data
    for feature in geojson_data["features"]:
        # Extract the properties of the feature
        properties = feature["properties"]

        # Extract the geometry (coordinates) of the feature
        geometry = feature["geometry"]
        coordinates = geometry["coordinates"] if geometry and "coordinates" in geometry else None

        # Create a dictionary for the current feature with property names as keys
        feature_dict = {key: value for key, value in properties.items()}

        # Add geometry coordinates to the dictionary
        feature_dict["coordinates"] = coordinates

        # Add the USGS ID from the GeoJSON data
        feature_dict["id"] = feature["id"]

        # Append the dictionary to the list
        earthquakes.append(feature_dict)
    return earthquakes


def get_event_rake(event_id):
    """
    Fetches the rake angles for a specific event ID from USGS. Example : [-170.21, -34.16]
    :param event_id: The USGS event ID for the earthquake.
    :return: A list of rake angles from both nodal planes, or an empty list if not available.
    """
    detail_url = "https://earthquake.usgs.gov/fdsnws/event/1/query"
    params = {"eventid": event_id, "format": "geojson"}

    print(f"Fetching rake for event_id: {event_id}")
    try:
        response = requests.get(detail_url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        # Access products
        products = data.get("properties", {}).get("products", {})

        candidates = products.get("moment-tensor", []) + products.get("focal-mechanism", [])

        if not candidates:
            print(f"  No moment-tensor or focal-mechanism products found for {event_id}.")
            return []

        # Iterate through ALL candidates until we find one with rake data
        for product in candidates:
            props = product.get("properties", {})

            # Check if this product has the nodal plane info
            r1 = props.get("nodal-plane-1-rake")
            r2 = props.get("nodal-plane-2-rake")

            # If both exist, we found a valid product
            if r1 is not None and r2 is not None:
                try:
                    rakes = [float(r1), float(r2)]
                    print(f"  Found rakes in product {product.get('code')}: {rakes}")
                    return rakes
                except ValueError:
                    continue

        print(
            f"  Warning: products found, but no 'nodal-plane-X-rake' properties present for {event_id}."
        )
        return []

    except Exception as e:
        print(f"Warning: Could not fetch rake for {event_id}: {e}")
        return []


def parse_custom_eq_list(file_path):
    """
    Parses a custom JSON list of earthquakes and converts them into the standard
    dictionary format expected by the coseis.py processing pipeline.
    :param file_path: Path to the custom earthquake list JSON file
    :return: List of earthquake dictionaries with standardized keys
    """
    print("=========================================")
    print(f"Loading custom earthquake list from: {file_path}")
    print("=========================================")

    with open(file_path, "r") as f:
        raw_data = json.load(f)

    earthquakes = []
    for item in raw_data:
        title = item.get("title", "Unknown_Event")

        # Safely extract epicenter data
        epi = item.get("epicenter", {})
        lon = epi.get("longitude")
        lat = epi.get("latitude")
        depth = epi.get("depth_km")

        if lon is None or lat is None:
            print(f"Skipping '{title}' - Missing coordinates.")
            continue

        # Extract USGS Event ID from URL if present
        url = epi.get("usgs_event_url", "")
        eq_id = url.split("/")[-1] if url else f"custom_{int(time.time())}"

        # Convert time string to Unix timestamp in milliseconds
        time_str = item.get("time")
        try:
            # Handle "YYYY-MM-DD HH:MM:SS UTC"
            clean_time = time_str.replace(" UTC", "")
            dt = datetime.strptime(clean_time, "%Y-%m-%d %H:%M:%S").replace(tzinfo=timezone.utc)
            time_ms = int(dt.timestamp() * 1000)
        except Exception as e:
            print(f"Skipping '{title}' - Time parsing error: {e}")
            continue

        # Construct the standardized dictionary
        eq = {
            "title": title,
            "coordinates": [lon, lat, depth],
            "time": time_ms,
            "id": eq_id,
            "url": url,
        }
        earthquakes.append(eq)
        print(f"Loaded: {title} ({eq_id})")

    return earthquakes
