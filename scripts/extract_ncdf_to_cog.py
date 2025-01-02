import os
import numpy as np
from gdal2tiles import generate_tiles
import geojson
import matplotlib.pyplot as plt
from matplotlib import colors
import rasterio
from shapely import wkt
from shapely.geometry import mapping
import subprocess
import shutil


def extract_footprint(netCDF_path, output_dir):
    """
    Function to extract the productBoundingBox layer from a NetCDF file and save it as a GeoJSON.
    """
    # Remove the main 'tiles' directory, if it already exists, and create a new one
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    #Process the productBoundingBox as a GeoJSON
    bbox_layer = "productBoundingBox"
    bbox_subdataset = f"NETCDF:\"{netCDF_path}\":{bbox_layer}"
    print(bbox_subdataset)
    print(f"Processing bounding box layer: {bbox_layer}")

    with rasterio.open(bbox_subdataset) as bbox_src:
        wkt_string = bbox_src.tags()['NC_GLOBAL#product_geometry_wkt']

        if wkt_string:
            # Parse the WKT string into a shapely Polygon object
            bounding_box = wkt.loads(wkt_string)
        else:
            raise ValueError("Bounding box WKT string not found in the metadata.")

        # Create a GeoJSON FeatureCollection with the bounding box
        geojson_feature = geojson.Feature(geometry=mapping(bounding_box), properties={})
        geojson_output = geojson.FeatureCollection([geojson_feature])

        # Write GeoJSON to file
        geojson_output_path = f"{output_dir}/{bbox_layer}.geojson"
        with open(geojson_output_path, "w") as geojson_file:
            geojson.dump(geojson_output, geojson_file)

        print(f"Saved GeoJSON: {geojson_output_path}")

    return


def tile_raster(input_raster, output_dir):
    """
    Function to tile a raster using gdal2tiles.py
    """

    # Remove the output directory if it already exists and create a new one
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    # Define the options for gdal2tiles
    options = {'zoom': (1, 10),
        'nb_processes': 4,
        }
    
    # Generate the tiles
    generate_tiles(input_raster, output_dir,
     **options)


def colorize_netCDF_layer_tiles(netcdf_path, output_dir):
    """
    Function to produce single-band cloud optimized GeoTIFFs from a NetCDF sublayers and a tiled colorized version.
    """
    # Remove the main 'tiles' directory, if it already exists, and create a new one
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    with rasterio.open(netcdf_path) as dataset:
        subdatasets = dataset.subdatasets
        if not subdatasets:
            print(f"No subdatasets found in {netcdf_path}")
            return

    # Define the rasters to process
    rasters = ['amplitude', 'azimuthPixelOffsets', 'rangePixelOffsets', 'unfilteredCoherence', 'unwrappedPhase']
    connected_components = f"NETCDF:\"{netcdf_path}\":/science/grids/data/connectedComponents"

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    for raster in rasters:
        subdataset_name = f"NETCDF:\"{netcdf_path}\":/science/grids/data/{raster}"
        print(f"Processing layer: {raster}")

        # Determine the colormap based on the raster layer
        if raster in ['azimuthPixelOffsets', 'rangePixelOffsets', 'unwrappedPhase']:
            cmap = plt.get_cmap("RdYlBu")
            is_grayscale = False
        else:
            cmap = plt.get_cmap("Greys_r")
            is_grayscale = True

        # Make a color table text file for input into gdalem
        with rasterio.open(subdataset_name) as src:
            data = src.read(1) 
            nodata_value = src.nodata
            if nodata_value is not None:
                data = np.ma.masked_equal(data, nodata_value)

        # Mask the data with connected components to ensure proper color scaling
        with rasterio.open(connected_components) as src:
            conn = src.read(1)
            msk = np.array(conn)*-1
            msk[msk<0] = 1
            data = np.ma.masked_array(data, mask=msk == 0)

        # Calculate 2nd and 98th percentiles for visualization
        p2 = np.percentile(data.compressed(), 2)
        p98 = np.percentile(data.compressed(), 98)
        print(f"{subdataset_name}: p2 = {p2}, p98 = {p98}")

        # Create a normalized colormap from the percentile range
        norm = colors.Normalize(vmin=p2, vmax=p98)

        # Generate color table content
        color_table_lines = []

        for value in np.linspace(p2, p98, 256):
            if is_grayscale:
                intensity = int(255 * (value - p2) / (p98 - p2))
                # Replace 0 with 1 in grayscale intensity
                intensity = max(intensity, 1)
                color_table_lines.append(f"{value:.2f} {intensity} {intensity} {intensity}")
            else:
                rgba = cmap(norm(value))
                rgb = tuple(max(int(c * 255), 1) for c in rgba[:3])  # Normalize to 0-255 and replace 0 with 1
                color_table_lines.append(f"{value:.2f} {rgb[0]} {rgb[1]} {rgb[2]}")

        # Handle nodata values by adding "nodata" for the nodata range
        if nodata_value is not None:
            color_table_lines.insert(0, f"{nodata_value:.2f} nodata nodata nodata")

        # Construct output paths
        color_table_file = os.path.join(output_dir, f"{raster}_color_table.txt")
        single_band_output = os.path.join(output_dir, f"{raster}_data.tif")
        output_colorized = os.path.join(output_dir, f"{raster}_colorized.tif")
        output_colorized_ds = os.path.join(output_dir, f"{raster}_colorized_ds.tif")
        tiled_output_dir = os.path.join(output_dir, f"{raster}_tiles")

        # Save the color table to the output directory
        with open(color_table_file, "w") as f:
            f.write("\n".join(color_table_lines))
        print(f"Color table saved as {color_table_file}")

        # Create a single-band output GeoTIFF from the original data
        subprocess.run(
            [
                "gdal_translate", subdataset_name, single_band_output,
                "-a_nodata", str(nodata_value), # Set nodata value
                "-co", "COMPRESS=DEFLATE",      # Use DEFLATE compression
                "-co", "PREDICTOR=2",           # Use horizontal differencing predictor
                "-co", "TILED=YES"              # Enable tiling for better performance
            ],
            check=True
        )

        # Run gdaldem color-relief to generate the colorized raster with three bands (RGB)
        subprocess.run(
            [
                "gdaldem", "color-relief", subdataset_name, color_table_file, output_colorized,
                "-co", "COMPRESS=DEFLATE",       # Use DEFLATE compression
                "-co", "PREDICTOR=2",            # Use horizontal differencing predictor
                "-co", "TILED=YES"               # Enable tiling for better performance
            ],
            check=True
        )

        # Assign nodata value to the colorized version
        subprocess.run(
            [
                "gdal_edit.py", "-a_nodata", '0.0', output_colorized
            ],
            check=True
        )
        
        # Reduce resolution of the colorized version
        output_resolution_x = 0.00008 # in decimal degrees (approx. 90 m)
        output_resolution_y = 0.00008 # in decimal degrees (approx. 90 m)

        subprocess.run(
            [
                "gdal_translate", output_colorized, output_colorized_ds,
                "-tr", str(output_resolution_x), str(output_resolution_y),  # Set resolution
                "-a_nodata", str(nodata_value),  # Set nodata value
                "-co", "COMPRESS=DEFLATE",       # Use DEFLATE compression
                "-co", "PREDICTOR=2",            # Use horizontal differencing predictor
                "-co", "TILED=YES"               # Enable tiling for better performance
            ],
            check=True
        )

        # Tile the colorized version using gdal2tiles
        tile_raster(output_colorized_ds, tiled_output_dir)

        # Delete intermediate files (except the final output)
        os.remove(color_table_file)
        os.remove(output_colorized)
        os.remove(output_colorized_ds)

    return


def colorize_netCDF_layer_COG(netcdf_path, output_dir):
    """
    Function to produce 4-band cloud optimized GeoTIFFs from a NetCDF sublayers.
    """
    # Remove the main 'cogs' directory, if it already exists, and create a new one
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    with rasterio.open(netcdf_path) as dataset:
        subdatasets = dataset.subdatasets
        if not subdatasets:
            print(f"No subdatasets found in {netcdf_path}")
            return

    # Define the rasters to process
    rasters = ['amplitude', 'azimuthPixelOffsets', 'rangePixelOffsets', 'unfilteredCoherence', 'unwrappedPhase']
    connected_components = f"NETCDF:\"{netcdf_path}\":/science/grids/data/connectedComponents"

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    for raster in rasters:
        subdataset_name = f"NETCDF:\"{netcdf_path}\":/science/grids/data/{raster}"
        print(f"Processing layer: {raster}")

        # Determine the colormap based on the raster layer
        if raster in ['azimuthPixelOffsets', 'rangePixelOffsets', 'unwrappedPhase']:
            cmap = plt.get_cmap("RdYlBu")
            is_grayscale = False
        else:
            cmap = plt.get_cmap("Greys_r")
            is_grayscale = True

        # Make a color table text file for input into gdalem
        with rasterio.open(subdataset_name) as src:
            data = src.read(1) 
            nodata_value = src.nodata

            if nodata_value is not None:
                data = np.ma.masked_equal(data, nodata_value)

        # Mask the data with connected components to ensure proper color scaling
        with rasterio.open(connected_components) as src:
            conn = src.read(1)
            msk = np.array(conn)*-1
            msk[msk<0] = 1
            data = np.ma.masked_array(data, mask=msk == 0)

        # Calculate 2nd and 98th percentiles for visualization
        p2 = np.percentile(data.compressed(), 2)
        p98 = np.percentile(data.compressed(), 98)
        print(f"{subdataset_name}: p2 = {p2}, p98 = {p98}")

        # Create a normalized colormap from the percentile range
        norm = colors.Normalize(vmin=p2, vmax=p98)

        # Generate color table content
        color_table_lines = []

        for value in np.linspace(p2, p98, 256):
            if is_grayscale:
                intensity = int(255 * (value - p2) / (p98 - p2))
                # Replace 0 with 1 in grayscale intensity
                intensity = max(intensity, 1)
                color_table_lines.append(f"{value:.2f} {intensity} {intensity} {intensity}")
            else:
                rgba = cmap(norm(value))
                rgb = tuple(max(int(c * 255), 1) for c in rgba[:3])  # Normalize to 0-255 and replace 0 with 1
                color_table_lines.append(f"{value:.2f} {rgb[0]} {rgb[1]} {rgb[2]}")

        # Handle nodata values by adding "nodata" for the nodata range
        if nodata_value is not None:
            color_table_lines.insert(0, f"{nodata_value:.2f} nodata nodata nodata")

        # Construct output paths
        color_table_file = os.path.join(output_dir, f"{raster}_color_table.txt")
        output_colorized = os.path.join(output_dir, f"{raster}_colorized.tif")
        output_colorized_modified = os.path.join(output_dir, f"{raster}_colorized_modified.tif")
        colorized_float32 = os.path.join(output_dir, f"{raster}_colorized_float32.tif")
        vrt_file = os.path.join(output_dir, f"{raster}_temp.vrt")
        final_output = os.path.join(output_dir, f"{raster}.tif")

        # Save the color table to the output directory
        with open(color_table_file, "w") as f:
            f.write("\n".join(color_table_lines))
        print(f"Color table saved as {color_table_file}")

        # Run gdaldem color-relief to generate the colorized raster with three bands (RGB)
        subprocess.run(
            [
                "gdaldem", "color-relief", subdataset_name, color_table_file, output_colorized,
                "-co", "COMPRESS=DEFLATE",       # Use DEFLATE compression
                "-co", "PREDICTOR=2",            # Use horizontal differencing predictor
                "-co", "TILED=YES"               # Enable tiling for better performance
            ],
            check=True
        )

        # Convert the color bands (RGB) to float32
        subprocess.run(
            [
                "gdal_translate", output_colorized, colorized_float32,
                "-ot", "Float32",               # Convert to float32
                "-co", "COMPRESS=DEFLATE",      # Use DEFLATE compression
                "-co", "PREDICTOR=2",           # Use horizontal differencing predictor
                "-co", "TILED=YES"              # Enable tiling for better performance
            ],
            check=True
        )

        # Modify pixel values of 0 to nodata_value in any band, which is needed for merging
        subprocess.run(
            [
                "gdal_calc.py",
                "-A", f"{colorized_float32}",  # Input raster
                "--calc", f"where(A==0, {nodata_value}, A)",  # Apply the conditional operation
                "--NoDataValue", str(nodata_value),  # Set NoData value dynamically
                "--allBands=A",                      # Process all bands from the input raster
                "--outfile", output_colorized_modified,  # Specify the output raster
                "--co", "COMPRESS=DEFLATE",          # Use DEFLATE compression
                "--co", "PREDICTOR=2",               # Use horizontal differencing predictor
                "--co", "TILED=YES"                  # Enable tiling for better performance
            ],
            check=True
        )

        # Create individual 1-band VRTs for each of the colorized bands (RGB)
        rgb_vrts = []
        for band_idx in range(1, 4):  # Bands 1, 2, 3 for RGB
            vrt_band_file = os.path.join(output_dir, f"{raster}_colorized_band{band_idx}.vrt")
            subprocess.run(
                [
                    "gdal_translate",
                    "-of", "VRT", 
                    "-b", str(band_idx), 
                    output_colorized_modified, 
                    vrt_band_file
                ],
                check=True
            )
            rgb_vrts.append(vrt_band_file)

        # Now build the 4-band VRT, including the original subdataset (band 1) and the 3 RGB bands
        subprocess.run(
            [
                "gdalbuildvrt", "-separate", "-overwrite", 
                vrt_file, subdataset_name, *rgb_vrts
            ],
            check=True
        )

        # Create the final 4-band GeoTIFF as a Cloud Optimized GeoTIFF (COG)
        subprocess.run(
            [
                "gdal_translate", vrt_file, final_output,
                "-a_nodata", str(nodata_value),     # Set nodata value
                "-co", "COMPRESS=DEFLATE",          # DEFLATE compression
                "-co", "PREDICTOR=2",               # Horizontal differencing predictor
                "-co", "TILED=YES",                 # Enable tiling (required for COG)
                "-co", "BIGTIFF=YES",               # Ensures support for large files
            ],
            check=True
        )

        print(f'Final COG saved as {final_output}')
        # Delete intermediate files (except the final output)
        temp_files = [color_table_file, output_colorized, colorized_float32, output_colorized_modified, vrt_file, *rgb_vrts]

        # Remove all temporary files
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                print(f"Deleted temporary file: {temp_file}")
    
    return


def run(nc, outdir_cogs, outdir_tiles, outdir_footprint):
    print("=====================================")
    print("Making COGs")
    print("=====================================")
    colorize_netCDF_layer_COG(nc, outdir_cogs)
    print("=====================================")
    print("Making tiles")
    print("=====================================")
    colorize_netCDF_layer_tiles(nc, outdir_tiles)
    print("=====================================")
    print("Extracting footprint")
    print("=====================================")
    extract_footprint(nc, outdir_footprint)
    print("=====================================")
    print("Derivatives generated successfully")
    print("=====================================")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process netCDF files.")
    parser.add_argument('nc', help="Path to the netCDF file")
    parser.add_argument('outdir_cogs', help="Output directory for COGs")
    parser.add_argument('outdir_tiles', help="Output directory for tiles")
    parser.add_argument('outdir_footprint', help="Output directory for footprints")

    args = parser.parse_args()

    run(args.nc, args.outdir_cogs, args.outdir_tiles, args.outdir_footprint)