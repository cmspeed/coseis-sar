"""Sentinel-1 SLC search on ASF DAAC: intersecting tracks/frames and SLCs per track."""

from dateutil import parser as dateparser
from dateutil.parser import isoparse
import requests
import geojson
import geopandas as gpd
from shapely.geometry import shape
from datetime import datetime, timedelta, timezone
from collections import defaultdict
from time import sleep

from aria_coseis.config import ASF_DAAC_API
from aria_coseis.utils import convert_time


def get_path_and_frame_numbers(AOI, time):
    """
    Query the ASF DAAC API for SLC data intersecting the Area of Interest (AOI) over a +/- 90 day window.
    This ensures all possible intersecting tracks are returned for the given AOI, avoiding data gap omissions.
    :param AOI: Shapely Polygon object representing the Area of Interest
    :param time: Unix timestamp representing the earthquake's origin time
    :return: Dictionary containing the path and frame numbers for each *unique* intersecting SLC.
    """
    # Establish the date range for the query
    rupture_date = convert_time(time)

    # Widen search to +/- 90 days to capture ALL intersecting tracks
    start_date = rupture_date - timedelta(days=90)
    start_date = start_date.replace(hour=0, minute=0, second=0)

    end_date = rupture_date + timedelta(days=90)
    end_date = end_date.replace(hour=23, minute=59, second=59)

    # Format the datetime object into a string
    start_date = start_date.strftime("%Y-%m-%dT%H:%M:%SZ")
    end_date = end_date.strftime("%Y-%m-%dT%H:%M:%SZ")

    # Define the ASF query parameters
    params = {
        "intersectsWith": AOI.wkt,
        "dataset": "SENTINEL-1",
        "processingLevel": "SLC",
        "beamSwath": "IW",
        "start": start_date,
        "end": end_date,
        "output": "geojson",
    }

    print(
        "Performing ASF DAAC API query to return path and frame numbers for SLCs intersecting AOI over the preceding 24 days..."
    )

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
            if not data.get("features"):
                print("No SLC features found intersecting the AOI.")
                return {}, gpd.GeoDataFrame()

            # Convert the GeoJSON data to a GeoDataFrame for visualization
            # Wrapped in try/except to handle cases where features lack geometry
            try:
                frame_dataframe = gpd.GeoDataFrame.from_features(data["features"], crs="EPSG:4326")
            except Exception as e:
                print(f"Warning: Could not create GeoDataFrame (likely missing geometry): {e}")
                return {}, gpd.GeoDataFrame()

            # Initialize an empty dictionary to store the path and frame numbers as sets
            path_frame_numbers = defaultdict(lambda: defaultdict(set))

            # Extract the path and frame numbers from the GeoJSON data
            for feature in data["features"]:
                flight_direction = feature["properties"]["flightDirection"]
                start_time = (
                    dateparser.isoparse(feature["properties"]["startTime"]).strftime(
                        "%Y-%m-%d %H:%M:%S"
                    )
                    + " UTC"
                )
                path_number = feature["properties"]["pathNumber"]
                frame_number = feature["properties"]["frameNumber"]
                path_frame_numbers[flight_direction][path_number].add(
                    (frame_number, start_time)
                )  # Use a set to avoid duplicates

            # Reformat the dictionary into usable format with tuples as keys and lists as values
            reformatted = {
                (flight_direction, path_number): sorted(frame_numbers)
                for flight_direction, path_frame in path_frame_numbers.items()
                for path_number, frame_numbers in path_frame.items()
            }

            print("=========================================")
            print("Reformatted Path and Frame Numbers for SLCs:")
            print("=========================================")
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


def get_SLCs(flight_direction, path_number, aoi_wkt, time, processing_mode):
    """
    Query the ASF DAAC API for SLC data based on the given path and AOI.
    The data are organized by flight direction and path number.
    :param processing_mode: 'historic', 'forward'
    :param flight_direction: 'ASCENDING' or 'DESCENDING'
    :param path_number: Sentinel-1 path number
    :param aoi_wkt: WKT string representation of the Area of Interest
    :param time: Unix timestamp representing the earthquake's origin time
    :return: List of dictionaries containing SLC fileIDs, dates, and geometries
    """
    # Establish the date range for the query
    rupture_date = convert_time(time)

    if processing_mode == "historic":
        start_date = rupture_date - timedelta(days=90)  # 90 days before the earthquake
        start_date = start_date.replace(hour=0, minute=0, second=0)
        end_date = rupture_date + timedelta(days=30)  # 30 days after the earthquake
        end_date = end_date.replace(hour=23, minute=59, second=59)
    elif processing_mode == "forward":
        # In forward mode for tracking, we want data acquired AFTER the earthquake
        start_date = rupture_date
        today = datetime.now(timezone.utc)
        end_date = today

    # Format the datetime object into a string
    start_date = start_date.strftime("%Y-%m-%dT%H:%M:%SZ")
    end_date = end_date.strftime("%Y-%m-%dT%H:%M:%SZ")

    # Define the query parameters
    params = {
        "flightDirection": flight_direction,
        "relativeOrbit": path_number,
        "intersectsWith": aoi_wkt,
        "dataset": "SENTINEL-1",
        "processingLevel": "SLC",
        "beamSwath": "IW",
        "start": start_date,
        "end": end_date,
        "output": "geojson",
    }

    print(
        f"Performing ASF DAAC API query to return SLCs for {flight_direction} path {path_number} intersecting the AOI..."
    )

    MAX_RETRIES = 10
    WAIT_SECONDS = 30

    # Specify cutoff date for filtering S1C and S1D granules based on the acquisition date (DockerizedTopsApp requirement)
    S1C_CUTOFF = datetime(2025, 5, 19, tzinfo=timezone.utc)
    S1D_CUTOFF = datetime(2026, 6, 24, tzinfo=timezone.utc)

    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(ASF_DAAC_API, params=params, timeout=160)
            response.raise_for_status()
            data = geojson.loads(response.text)

            SLCs = []
            for feature in data["features"]:
                start_time = feature["properties"]["startTime"]
                path = feature["properties"]["pathNumber"]
                frame = feature["properties"]["frameNumber"]
                file_id = feature["properties"]["fileID"]

                try:
                    dt_obj = isoparse(start_time)
                    date = dt_obj.strftime("%Y-%m-%dT%H:%M:%SZ")
                except Exception:
                    print(f"Warning: Unexpected date format in startTime: {start_time}")
                    date = None
                    dt_obj = None

                # Apply S1C and S1D filters
                if file_id.startswith("S1C") and dt_obj and dt_obj < S1C_CUTOFF:
                    continue
                if file_id.startswith("S1D") and dt_obj and dt_obj < S1D_CUTOFF:
                    continue

                SLC = {
                    "fileID": file_id,
                    "date": date,
                    "pathNumber": path,
                    "frameNumber": frame,
                    "geometry": shape(feature["geometry"]),
                }
                SLCs.append(SLC)

            print("=========================================")
            print(
                f"Found {len(SLCs)} valid SLCs for the {flight_direction} path {path_number} intersecting the AOI."
            )
            print("=========================================")
            return SLCs

        except requests.exceptions.RequestException as e:
            print(f"Request error from ASF DAAC API: {e}")
            if attempt < MAX_RETRIES - 1:
                print(f"Retrying in {WAIT_SECONDS} seconds...")
                sleep(WAIT_SECONDS)
            else:
                print("All retry attempts failed.")
                return None
