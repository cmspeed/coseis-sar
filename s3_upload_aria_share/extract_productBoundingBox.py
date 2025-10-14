#!/usr/bin/env python3

import os
import json
import re
from osgeo import gdal, ogr

# ===== CONFIGURATION =====
LOCAL_ROOT = "S3_local_copy"               # where your .nc files live
OUTPUT_DIR = "bounding_boxes_geojson"     # folder for per-file GeoJSONs
MASTER_GEOJSON = "all_bounding_boxes.geojson"  # aggregated output
SUBDATASET = "productBoundingBox"         # Sublayer name for the vector geometry
# =========================

def extract_earthquake_info(nc_path):
    """
    Extracts the earthquake magnitude and the full earthquake name (directory).
    The information is in the GREAT-GRANDPARENT directory name.
    
    Returns: A tuple (magnitude_float, earthquake_name_str) or (None, str).
    """
    # 1. Get the immediate parent directory (S1-COSEIS_...)
    parent_dir = os.path.dirname(nc_path)
    
    # 2. Get the grandparent directory (the product directory, e.g., S1-COSEIS_...)
    grandparent_dir = os.path.dirname(parent_dir)
    
    # 3. Get the great-grandparent directory name (the earthquake event name, e.g., m_61_...)
    earthquake_dir_name = os.path.basename(grandparent_dir)
    
    # --- Extract Magnitude ---
    # Use a regular expression to find the number immediately after 'm_'
    match = re.search(r'm_(\d+)', earthquake_dir_name)
    
    earthquake_magnitude = None
    if match:
        magnitude_int = int(match.group(1))
        # Convert the integer (e.g., 60) to a float (e.g., 6.0)
        earthquake_magnitude = float(magnitude_int) / 10.0
    else:
        print(f"Warning: Could not find magnitude pattern 'm_XX' in directory name: {earthquake_dir_name}")
        
    # --- Return both pieces of information ---
    # The earthquake name is always the directory name.
    return earthquake_magnitude, earthquake_dir_name


def extract_bounding_box(nc_path):
    """
    Extract the productBoundingBox geometry and add the magnitude and
    earthquake name properties.
    """
    # 1. Define the OGR DataSource string for the specific subdataset/layer
    ogr_ds_name = f'NETCDF:"{nc_path}":{SUBDATASET}'

    # 2. Open the OGR DataSource
    ds = ogr.Open(ogr_ds_name)
    if ds is None:
        raise Exception(f"Failed to open OGR layer: {ogr_ds_name}")

    # 3. Get the layer
    layer = ds.GetLayer(0)
    if layer is None or layer.GetFeatureCount() == 0:
        ds = None
        raise Exception("Layer is empty or could not be retrieved.")

    # 4. Get the first (and only) feature
    feature = layer.GetNextFeature()
    if feature is None:
        ds = None
        raise Exception("Could not retrieve the first feature from the layer.")
        
    # 5. Get the geometry object
    geom = feature.GetGeometryRef()
    if geom is None:
        ds = None
        raise Exception("Feature has no geometry reference.")

    # 6. Convert the OGR Geometry to a GeoJSON geometry dictionary
    geom_json = json.loads(geom.ExportToJson()) 

    # Clean up
    feature = None
    ds = None 
    
    # 7. Extract the magnitude and earthquake name
    earthquake_magnitude, earthquake_name = extract_earthquake_info(nc_path)
    
    # 8. Create the GeoJSON Feature object, including the new fields
    properties = {
        "filename": os.path.basename(nc_path),
        "source_layer": SUBDATASET,
        # Add the new 'earthquake_name' field
        "earthquake_name": earthquake_name
    }
    
    # Only add 'magnitude' if it was successfully extracted
    if earthquake_magnitude is not None:
        properties["magnitude"] = earthquake_magnitude

    return {
        "type": "Feature",
        "geometry": geom_json,
        "properties": properties
    }

# The save_feature and main functions from your previous code remain valid:

def save_feature(feature, nc_path):
    """Save a single feature as its own GeoJSON, mirroring directory structure."""
    rel_path = os.path.relpath(nc_path, LOCAL_ROOT)
    filename_base = os.path.splitext(os.path.basename(nc_path))[0]
    # rel_dir is the path from LOCAL_ROOT to the directory containing the magnitude/earthquake name
    rel_dir = os.path.dirname(os.path.dirname(rel_path))

    geojson_path = os.path.join(
        OUTPUT_DIR,
        rel_dir,
        os.path.basename(os.path.dirname(rel_path)) + "_" + filename_base + ".geojson"
    )

    os.makedirs(os.path.dirname(geojson_path), exist_ok=True)

    with open(geojson_path, "w") as f:
        json.dump({"type": "FeatureCollection", "features": [feature]}, f, indent=2)

    print(f"Saved {geojson_path}")
    return feature

def main():
    all_features = []
    
    if not os.path.exists(LOCAL_ROOT):
        print(f"Error: The directory '{LOCAL_ROOT}' does not exist.")
        print("Please ensure your .nc files are in the specified LOCAL_ROOT directory.")
        return

    for root, _, files in os.walk(LOCAL_ROOT):
        for fn in files:
            if fn.endswith(".nc"):
                nc_path = os.path.join(root, fn)
                try:
                    feature = extract_bounding_box(nc_path)
                    feature = save_feature(feature, nc_path)
                    all_features.append(feature)
                except Exception as e:
                    print(f"Failed to process {nc_path}: {e}")

    # Save aggregated GeoJSON
    if all_features:
        with open(MASTER_GEOJSON, "w") as f:
            json.dump({"type": "FeatureCollection", "features": all_features}, f, indent=2)
        print(f"\nMaster GeoJSON written to {MASTER_GEOJSON} with {len(all_features)} features.")
    else:
        print("\nNo features extracted, master GeoJSON not created.")

    print("\nAll done.")

if __name__ == "__main__":
    main()