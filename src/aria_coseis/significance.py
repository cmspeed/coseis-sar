"""Significance filtering: magnitude, depth, proximity to land, and the optical rake filter."""

import requests
import json
import geojson
import csv
from shapely.geometry import mapping, Point, Polygon, LineString, MultiPolygon
from shapely.ops import linemerge

from aria_coseis.config import coastline_api
from aria_coseis.usgs import get_event_rake
from aria_coseis.utils import convert_time


def get_coastline(coastline_api):
    """
    Fetch coastline data from OSGEO/PROJ Github repo and return it as a GeoJSON object.
    The data returned will be in the form of a MultiPolygon covering the landmass interiors
    and a 0.5 degree ocean buffer. The data are also written to "coastline_buffered.geojson".

    :param coastline_api: API endpoint for the coastline data
    :return: GeoJSON object containing coastline data
    """
    try:
        # Fetch data from the specified API
        response = requests.get(coastline_api)
        response.raise_for_status()

        # Parse the response as GeoJSON
        coastline_data = geojson.loads(response.text)

        # Extract the LineString features from the GeoJSON data
        features = []
        for feature in coastline_data["features"]:
            if feature["geometry"]["type"] == "LineString":
                features.append(LineString(feature["geometry"]["coordinates"]))
            else:
                print("Coastline data is not in LineString format.")

        # Merge contiguous line segments together to form complete coastlines
        merged_lines = linemerge(features)

        # Ensure merged_lines is iterable (handles cases where linemerge returns a single LineString)
        if isinstance(merged_lines, LineString):
            merged_lines = [merged_lines]
        else:
            merged_lines = merged_lines.geoms

        # Convert merged LineStrings to Polygons
        polygons = []
        for line in merged_lines:
            if len(line.coords) < 3:
                continue

            if line.is_ring:
                polygons.append(Polygon(line))
            else:
                # Close the LineString and create a Polygon
                closed_line = LineString(list(line.coords) + [line.coords[0]])
                polygons.append(Polygon(closed_line))

        # Combine all polygons into a MultiPolygon
        coastline_polys = MultiPolygon(polygons)

        # Buffer the entire MultiPolygon by 0.5 degrees
        coastline_buffered = coastline_polys.buffer(0.5)

        # Convert the MultiPolygon to GeoJSON format
        geojson_data = {
            "type": "FeatureCollection",
            "features": [
                {"type": "Feature", "geometry": mapping(coastline_buffered), "properties": {}}
            ],
        }

        # Save to a GeoJSON file
        output_file = "coastline_buffered.geojson"
        with open(output_file, "w") as f:
            json.dump(geojson_data, f, indent=2)

        return coastline_buffered

    except requests.RequestException as e:
        print(f"Error accessing coastline API: {e}")
        return None
    except json.JSONDecodeError as e:
        print(f"Error parsing GeoJSON data: {e}")
        return None


def withinCoastline(earthquake, coastline):
    """
    Determine if earthquake epicenter is within 0.5 decimal degrees (~55 km) of the coastline.
    This is one filtering parameter to determine if an earthquake is "significant" within the scope of this project.
    :param earthquake: dictionary containing earthquake data
    :param coastline: shapely Polygon object representing the coastline
    :return: True if the epicenter is within the coastline, False otherwise
    """
    # Extract the coordinates of the earthquake's epicenter
    coords = earthquake.get("coordinates", [])
    if not coords:
        return None

    # Create a Point object from the earthquake's coordinates
    epicenter = Point(coords[:2])

    # Buffer the coastline by 0.5 degrees
    coastline_buffer = coastline.buffer(0.5)

    # Determine if the epicenter is within the coastline
    within_coastline_buffer = coastline_buffer.contains(epicenter)
    return within_coastline_buffer


def check_significance(earthquakes, start_date, end_date=None, sensor="sar", mode="historic"):
    """
    Check the significance of each earthquake based on its
    (1) magnitude and (2) depth (historic: M>=6.0 and <=40 km; forward: M>=5.5 and <=15 km or M>=6.0 and <=40 km),
    (3) distance from land (within 0.5 degrees, ~55 km of the coastline),
    and (4) if sensor is optical, rake angle (must be strike-slip: ~0 or ~180 degrees).

    :param earthquakes: list of dictionaries containing earthquake data
    :param start_date: start date in the format 'YYYY-MM-DD'
    :param end_date: end date in the format 'YYYY-MM-DD' (optional)
    :param sensor: determines if rake filter will be applied ('sar' or 'optical')
    :param mode: 'historic' or 'forward' - determines the magnitude/depth criteria used for filtering.
    :return: List of dictionaries containing significant earthquakes
    """
    print("=========================================")
    print(f"Checking for significant earthquakes (Sensor: {sensor.upper()})...")
    print("=========================================")

    significant_earthquakes = []
    coastline = get_coastline(coastline_api)

    # Rake tolerance (degrees)
    RAKE_TOLERANCE = 45.0

    for earthquake in earthquakes:
        magnitude = earthquake.get("mag")
        depth = earthquake.get("coordinates", [])[2] if earthquake.get("coordinates") else None

        # Filters applicable to both sensors
        is_candidate = False

        if all(var is not None for var in (magnitude, depth)):
            within_Coastline_buffer = withinCoastline(earthquake, coastline)

            if mode == "historic":
                if (magnitude >= 6.0) and (depth <= 40.0) and within_Coastline_buffer:
                    is_candidate = True
            elif mode == "forward":
                # Catch M>=5.5 & <=15km OR M>=6.0 & <=40km
                is_shallow_moderate = (magnitude >= 5.5) and (depth <= 15.0)
                is_deeper_larger = (magnitude >= 6.0) and (depth <= 40.0)
                if within_Coastline_buffer and (is_shallow_moderate or is_deeper_larger):
                    is_candidate = True

        if not is_candidate:
            continue

        # SAR: no rake requirement
        if sensor not in ["sentinel-2", "landsat"]:
            significant_earthquakes.append(earthquake)
            continue

        # Optical: require strike-slip rake
        title = earthquake.get("title", "Unknown Event")
        print(f"  Fetching rake for candidate: {title}...")

        rakes = get_event_rake(earthquake.get("id"))
        earthquake["rakes"] = rakes

        if not rakes:
            print(f"    -> Skipped (Optical mode requires rake data, none found)")
            continue

        # Check if ANY available rake satisfies the condition
        is_strike_slip = False
        accepted_rake_val = None

        for r in rakes:
            if (abs(r) <= RAKE_TOLERANCE) or (abs(r) >= (180.0 - RAKE_TOLERANCE)):
                is_strike_slip = True
                accepted_rake_val = r
                break

        if is_strike_slip:
            print(f"    -> Accepted (Rake {accepted_rake_val}° fits strike-slip criteria)")
            significant_earthquakes.append(earthquake)
        else:
            print(f"    -> Skipped (Rakes {rakes} indicate dip-slip/oblique motion)")

    # Output / Logging
    if len(significant_earthquakes) > 0:
        print("=========================================")
        print(f"Found {len(significant_earthquakes)} significant earthquakes.")
        print("=========================================")
        for eq in significant_earthquakes:
            print(f"Name: {eq['title']}")
            print(f"Rupture Date/Time: {convert_time(eq['time'])} UTC")
            print(f"Magnitude: {eq['mag']}")
            print(f"Depth: {eq['coordinates'][2]} km")
            print(f"Alert Level: {eq['alert']}")
            if "rakes" in eq:
                print(f"Rakes: {eq['rakes']}")
            print("=========================================")

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
        dt = convert_time(eq["time"])

        # Create a GeoJSON feature for each earthquake
        feature = {
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": eq["coordinates"][:2]},
            "properties": {
                "place": eq["place"],
                "magnitude": eq["mag"],
                "date": dt.strftime("%Y-%m-%d"),
                "time_utc": dt.strftime("%H:%M:%S"),
                "longitude": eq["coordinates"][0],
                "latitude": eq["coordinates"][1],
                "depth_km": eq["coordinates"][2],
                "alert": eq["alert"],
                "url": eq["url"],
            },
        }
        geojson_features.append(feature)

    # Create the final GeoJSON structure
    geojson_data = {"type": "FeatureCollection", "features": geojson_features}

    # Output the data to GeoJSON and CSV files
    if start_date and end_date:
        with open(f"significant_earthquakes_{start_date}_to_{end_date}_M6.geojson", "w") as f:
            geojson.dump(geojson_data, f)

        with open(f"significant_earthquakes_{start_date}_to_{end_date}.csv", "w", newline="") as f:
            writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
            writer.writerow(
                [
                    "Place",
                    "Magnitude",
                    "Date",
                    "Time_utc",
                    "Longitude",
                    "Latitude",
                    "Depth_km",
                    "Alert",
                    "URL",
                ]
            )

            for eq in significant_earthquakes:
                dt = convert_time(eq["time"])
                writer.writerow(
                    [
                        eq["place"],
                        eq["mag"],
                        dt.strftime("%Y-%m-%d"),
                        dt.strftime("%H:%M:%S"),
                        eq["coordinates"][0],
                        eq["coordinates"][1],
                        eq["coordinates"][2],
                        eq["alert"],
                        eq["url"],
                    ]
                )
    else:
        with open(f"significant_earthquakes_{start_date}.geojson", "w") as f:
            geojson.dump(geojson_data, f)

        with open(f"significant_earthquakes_{start_date}.csv", "w", newline="") as f:
            writer = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
            writer.writerow(
                [
                    "Place",
                    "Magnitude",
                    "Date",
                    "Time_utc",
                    "Longitude",
                    "Latitude",
                    "Depth_km",
                    "Alert",
                    "URL",
                ]
            )

            for eq in significant_earthquakes:
                dt = convert_time(eq["time"])
                writer.writerow(
                    [
                        eq["place"],
                        eq["mag"],
                        dt.strftime("%Y-%m-%d"),
                        dt.strftime("%H:%M:%S"),
                        eq["coordinates"][0],
                        eq["coordinates"][1],
                        eq["coordinates"][2],
                        eq["alert"],
                        eq["url"],
                    ]
                )
    return
