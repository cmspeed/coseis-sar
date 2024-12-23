import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import rasterio
import tempfile
import subprocess


def extract_footprint(netCDF_path, output_dir):
    #Process the productBoundingBox as a GeoJSON
    bbox_layer = "productBoundingBox"
    bbox_subdataset = f"NETCDF:\"{netcdf_path}\":{bbox_layer}"
    print(f"Processing bounding box layer: {bbox_layer}")

    with rasterio.open(bbox_subdataset) as bbox_src:
        # Read the bounding box layer
        bbox_data = bbox_src.read(1)  # Assuming the layer has one band
        transform = bbox_src.transform
        crs = bbox_src.crs

        # Extract geometries
        shapes_generator = shapes(bbox_data, transform=transform)
        geojson_features = [
            {
                "type": "Feature",
                "geometry": mapping(shape(geom)),
                "properties": {"value": value},
            }
            for geom, value in shapes_generator
        ]

        # Create GeoJSON object
        geojson_output = {
            "type": "FeatureCollection",
            "features": geojson_features,
        }

        # Write GeoJSON to file
        geojson_output_path = f"{output_dir}/{bbox_layer}.geojson"
        with open(geojson_output_path, "w") as geojson_file:
            geojson.dump(geojson_output, geojson_file)

        print(f"Saved GeoJSON: {geojson_output_path}")


def colorize_netCDF_layer(netcdf_path, output_dir):
    """
    Function to colorize a layer using a gdalem.
    """
    with rasterio.open(netcdf_path) as dataset:
        subdatasets = dataset.subdatasets

        if not subdatasets:
            print(f"No subdatasets found in {netcdf_path}")
            return

        rasters = ['azimuthPixelOffsets', 'rangePixelOffsets', 'unfilteredCoherence', 'amplitude', 'unwrappedPhase']

        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        for raster in rasters:
            subdataset_name = f"NETCDF:\"{netcdf_path}\":/science/grids/data/{raster}"
            print(f"Processing layer: {raster}")

            # Determine the colormap based on the raster layer
            if raster in ['azimuthPixelOffsets', 'rangePixelOffsets', 'unwrappedPhase']:
                cmap = plt.get_cmap("RdBu")
            else:
                cmap = plt.get_cmap("gray")

            # Make a color table text file for input into gdalem
            with rasterio.open(subdataset_name) as src:
                data = src.read(1) 
                nodata_value = src.nodata

                if nodata_value is not None:
                    data = np.ma.masked_equal(data, nodata_value)

                # Calculate 2nd and 98th percentiles for visualization
                p2 = np.percentile(data.compressed(), 2)
                p98 = np.percentile(data.compressed(), 98)
                print(f"{subdataset_name}: p2 = {p2}, p98 = {p98}")

                # Create a normalized colormap from the percentile range
                norm = colors.Normalize(vmin=p2, vmax=p98)

                # Generate color table content
                color_table_lines = []

                for value in np.linspace(p2, p98, 256):
                    rgba = cmap(norm(value))
                    rgb = tuple(int(c * 255) for c in rgba[:3])  # Normalize to 0-255
                    color_table_lines.append(f"{value:.2f} {rgb[0]} {rgb[1]} {rgb[2]}")

                # Add colors for below p2 and above p98
                color_table_lines.insert(0, f"{p2:.2f} 0 0 255")  # Below p2
                color_table_lines.append(f"{p98:.2f} 255 0 0")    # Above p98

            # Construct output paths
            color_table_file = os.path.join(output_dir, f"{raster}_color_table.txt")
            output_colorized = os.path.join(output_dir, f"{raster}_colorized_deflate.tif")
            final_output = os.path.join(output_dir, f"{raster}_final.tif")


            # Save the color table to the output directory
            with open(color_table_file, "w") as f:
                f.write("\n".join(color_table_lines))

            print(f"Color table saved as {color_table_file}")

            # Use the text file with gdaldem, adding compression options
            subprocess.run(
                [
                    "gdaldem", "color-relief", subdataset_name, color_table_file, output_colorized,
                    "-co", "COMPRESS=DEFLATE",  # Use DEFLATE compression
                    "-co", "PREDICTOR=2",      # Use horizontal differencing predictor
                    "-co", "TILED=YES"         # Enable tiling for better performance
                ],
                check=True
            )


            # Combine original data and colorized raster into a multi-band GeoTIFF
            subprocess.run(
                [
                    "gdalbuildvrt", "-separate", "-overwrite", "temp.vrt",
                    subdataset_name, output_colorized
                ],
                check=True
            )
            subprocess.run(
                [
                    "gdal_translate", "temp.vrt", final_output,
                    "-co", "COMPRESS=DEFLATE", "-co", "PREDICTOR=2", "-co", "TILED=YES"
                ],
                check=True
            )
    
    return


def main():
    nc = '/u/trappist-r0/colespeed/work/coseis/earthquakes/Nevada_12-09-2024/A64/insar_20241216_20241204/S1-GUNW_CUSTOM-A-R-064-tops-20241216_20241204-015210-00120W_00038N-PP-ba76-v3_0_1/S1-GUNW_CUSTOM-A-R-064-tops-20241216_20241204-015210-00120W_00038N-PP-ba76-v3_0_1.nc'
    outdir = '/u/trappist-r0/colespeed/work/coseis/earthquakes/Nevada_12-09-2024/A64/insar_20241216_20241204/S1-GUNW_CUSTOM-A-R-064-tops-20241216_20241204-015210-00120W_00038N-PP-ba76-v3_0_1/cogs'
    #extract_netcdf_to_cog(nc, outdir)
    #extract_netcdf_to_cog_metadata(nc, outdir)
    #extract_netcdf_to_cog_w_gdal(nc, outdir)
    colorize_netCDF_layer(nc, outdir)

if __name__ == "__main__":
    main()
