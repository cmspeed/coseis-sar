"""Area-of-interest construction and loading."""

from __future__ import annotations

import json
from urllib.parse import urlparse

import requests
from shapely.geometry import Polygon, box, shape
from shapely.ops import unary_union


def make_aoi(coordinates: list[float]) -> Polygon:
    """
    Create an Area of Interest (AOI) polygon based on the given coordinates.
    The AOI is a square with a side length of 1 degree (~111 km) centered on the earthquake's epicenter.
    :param coordinates: list containing the longitude and latitude of the earthquake's epicenter
    :return: Shapely Polygon object representing the AOI
    """
    print("=========================================")
    print("Creating Area of Interest (AOI) polygon...")
    print("=========================================")

    # Extract the X and Y coordinates of the earthquake's epicenter
    X = coordinates[0]
    Y = coordinates[1]

    # Define the side length of the square AOI in decimal degrees
    side_length = 1.0  # 1 degree is ~111 km at the equator

    # Calculate half side length
    half_side = side_length / 2

    # Define the square's vertices relative to the center point
    square_coords = [
        (X - half_side, Y - half_side),  # LL
        (X + half_side, Y - half_side),  # LR
        (X + half_side, Y + half_side),  # UR
        (X - half_side, Y + half_side),  # UL
    ]

    # Create the square polygon
    AOI = Polygon(square_coords)

    print(f"Area of Interest (AOI) created: {AOI}")
    print("=========================================")
    return AOI


def load_aoi_from_json(aoi_path_or_url: str) -> Polygon:
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
            with open(aoi_path_or_url, "r") as f:
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
