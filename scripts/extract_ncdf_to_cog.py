import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import rasterio
import tempfile
import subprocess


def extract_netcdf_to_cog_w_gdal(netcdf_path, output_dir):
    """
    Processes subdatasets in a netCDF file to generate GeoTIFFs retaining original values but 
    with visualization metadata (color stretch and colormap).

    Parameters:
    - netcdf_path (str): Path to the netCDF file.
    - output_dir (str): Directory to save the outputs.
    """
    with rasterio.open(netcdf_path) as dataset:
        subdatasets = dataset.subdatasets

    if not subdatasets:
        print(f"No subdatasets found in {netcdf_path}")
        return

    rasters = ['azimuthPixelOffsets', 'rangePixelOffsets', 'unfilteredCoherence', 'amplitude', 'unwrappedPhase']

    for raster in rasters:
        subdataset_name = f"NETCDF:\"{netcdf_path}\":/science/grids/data/{raster}"
        print(f"Processing layer: {raster}")

        if raster == "unwrappedPhase":
            ds = gdal.Open(subdataset_name)

            # Read the data from the first band
            data = ds.GetRasterBand(1).ReadAsArray()
            nodata_value = ds.GetRasterBand(1).GetNoDataValue()

            # Mask nodata values
            if nodata_value is not None:
                data = np.ma.masked_equal(data, nodata_value)

            # Calculate percentiles for visualization
            p2 = np.percentile(data.compressed(), 2)
            p98 = np.percentile(data.compressed(), 98)
            print(f"{subdataset_name}: p2 = {p2}, p98 = {p98}")

            # Create a colormap from the percentile range
            cmap = plt.get_cmap("RdBu")
            norm = colors.Normalize(vmin=p2, vmax=p98)

            # Create the color table as a list of RGB tuples for the percentile range
            color_table = [
                tuple(int(c * 255) for c in cmap(norm(value))[:3]) + (255,)  # RGBA format
                for value in np.linspace(p2, p98, 256)
            ]

            # Extend the color table for out-of-range values (below p2 and above p98)
            color_table[0] = tuple(int(c * 255) for c in cmap(0)[:3]) + (255,)  # Below p2
            color_table[255] = tuple(int(c * 255) for c in cmap(1)[:3]) + (255,)  # Above p98

            # Output GeoTIFF path
            output_path = f"{output_dir}/{raster}_colorized_deflated_COG_w_gdal.tif"

            # Create a new GeoTIFF to store the colorized raster
            driver = gdal.GetDriverByName("GTiff")
            out_ds = driver.Create(
                output_path,
                ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32,  # Store as float32
                options=["COMPRESS=DEFLATE", "TILED=YES"]
            )

            # Set the projection and geotransform from the input dataset
            out_ds.SetProjection(ds.GetProjection())
            out_ds.SetGeoTransform(ds.GetGeoTransform())

            # Write the original data to the output file
            out_band = out_ds.GetRasterBand(1)
            out_band.WriteArray(data.filled(nodata_value))

            # Set the NoData value (if applicable)
            if nodata_value is not None:
                out_band.SetNoDataValue(nodata_value)

            # Add the colormap as metadata (for visualization purposes)
            # Store the color map as a list of RGBA values in the metadata
            out_ds.SetMetadata({
                'SCALE_MIN': str(p2),
                'SCALE_MAX': str(p98),
                'COLOR_TABLE': "RdBu",
                'COLOR_MAP': str(color_table),  # Store the color table as a string
            })

            # Close the dataset to save it
            out_ds = None
            ds = None

            print(f"Colorized raster written to {output_path}")

def extract_netcdf_to_cog_metadata(netcdf_path, output_dir):
    """
    Processes subdatasets in a netCDF file to generate GeoTIFFs retaining original values but 
    with visualization metadata (color stretch and colormap).

    Parameters:
    - netcdf_path (str): Path to the netCDF file.
    - output_dir (str): Directory to save the outputs.
    """
    with rasterio.open(netcdf_path) as dataset:
        subdatasets = dataset.subdatasets

    if not subdatasets:
        print(f"No subdatasets found in {netcdf_path}")
        return

    rasters = ['azimuthPixelOffsets', 'rangePixelOffsets', 'unfilteredCoherence', 'amplitude', 'unwrappedPhase']

    for raster in rasters:
        subdataset_name = f"NETCDF:\"{netcdf_path}\":/science/grids/data/{raster}"
        print(f"Processing layer: {raster}")

        if raster == "unwrappedPhase":
            with rasterio.open(subdataset_name) as src:
                data = src.read(1)  # Read the first band
                transform = src.transform
                crs = src.crs
                nodata_value = src.nodata

                if nodata_value is not None:
                    data = np.ma.masked_equal(data, nodata_value)

                # Calculate 2nd and 98th percentiles for visualization
                p2 = np.percentile(data.compressed(), 2)
                p98 = np.percentile(data.compressed(), 98)
                print(f"{raster}: p2 = {p2}, p98 = {p98}")

                # Create a colormap from the percentile range
                cmap = plt.get_cmap("RdBu")
                norm = colors.Normalize(vmin=p2, vmax=p98)

                # Create a color table for values within the percentile range
                color_table = {
                    i: tuple(int(c * 255) for c in cmap(norm(value))[:3]) + (255,)
                    for i, value in enumerate(np.linspace(p2, p98, 256))
                }

                # Extend the color table to handle values outside the range
                color_table[0] = tuple(int(c * 255) for c in cmap(0)[:3]) + (255,)  # Below p2
                color_table[255] = tuple(int(c * 255) for c in cmap(1)[:3]) + (255,)  # Above p98

                # Output GeoTIFF path
                output_path = f"{output_dir}/{raster}_colorized_deflated_COG.tif"

                # Write the GeoTIFF as a COG
                with rasterio.open(
                    output_path,
                    "w",
                    driver="GTiff",
                    height=data.shape[0],
                    width=data.shape[1],
                    count=1,
                    dtype="float32",  # Preserve original data values
                    crs=crs,
                    transform=transform,
                    nodata=nodata_value,
                    tiled=True,  # Enable tiling (required for COG)
                    compress="DEFLATE",  # Use DEFLATE compression
                    zlevel=9,  # Maximum compression level
                    blockxsize=256,
                    blockysize=256,
                ) as dst:
                    # Write the original data as float32
                    dst.write(data.filled(nodata_value), 1)

                    # Assign a color table (visualization only)
                    dst.colorinterp = [ColorInterp.palette]
                    dst.write_colormap(1, {i: color_table[i] for i in range(256)})

                    # Add overviews for COG compliance
                    overviews = [2, 4, 8, 16]  # Overview levels (can be adjusted)
                    dst.build_overviews(overviews, resampling=rasterio.enums.Resampling.nearest)

                    # Add COG-specific tags
                    dst.update_tags(
                        OVR_BLOCKSIZE=256  # Overview block size (COG requirement)
                    )

def extract_netcdf_to_cog_rgb(netcdf_path, output_dir, subset_filter=None):
    """
    Extract all subdatasets in a netCDF file into COGs.

    Parameters:
    - netcdf_path (str): Path to the netCDF file.
    - output_dir (str): Directory to save the output COGs.
    - subset_filter (callable, optional): A function to filter subdatasets by name. If None, processes all subdatasets.
    """
    # Open the netCDF file to list subdatasets
    with rasterio.open(netcdf_path) as dataset:
        subdatasets = dataset.subdatasets

    if not subdatasets:
        print(f"No subdatasets found in {netcdf_path}")
        return

    
    for subdataset in subdatasets:
        print('===============================')
        print('subdataset:', subdataset)
        print('===============================')

    rasters = ['azimuthPixelOffsets', 'rangePixelOffsets', 'unfilteredCoherence', 'amplitude', 'unwrappedPhase']

    # Extract raster layers
    for raster in rasters:
        subdataset_name = f"NETCDF:\"{netcdf_path}\":/science/grids/data/{raster}"

        # Open the subdataset
        with rasterio.open(subdataset_name) as src:
            data = src.read(1)  # Read the first band
            transform = src.transform
            crs = src.crs
            nodata_value = src.nodata  # Retrieve the nodata value

            # Mask the nodata values
            if nodata_value is not None:
                data = np.ma.masked_equal(data, nodata_value)

            # Clip coherence values to 0-1
            if raster == "unfilteredCoherence":
                data = np.clip(data, 0, 1)
            
            # Calculate 2nd and 98th percentiles
            vmin = np.percentile(data.compressed(), 2)  # Use compressed() to exclude masked values
            vmax = np.percentile(data.compressed(), 98)
            print(f"{raster}: 2nd percentile = {vmin}, 98th percentile = {vmax}")

            # Normalize the data to 0-1 for colormap application
            data_normalized = np.clip((data - vmin) / (vmax - vmin), 0, 1)

            # Apply colormap based on layer name
            if raster in ['unwrappedPhase', 'azimuthPixelOffsets', 'rangePixelOffsets']:
                cmap = plt.get_cmap("RdBu")
            else:
                cmap = plt.get_cmap("gray")

            rgb_data = cmap(data_normalized)[..., :3]  # Extract RGB channels

            # Convert to 8-bit integer
            rgb_data = (rgb_data * 255).astype("uint8")

            # Generate a unique output filename
            output_path = f"{output_dir}/{raster}.tif"
            print(output_path)
            # # Write the COG
            with rasterio.open(
                output_path,
                "w",
                driver="COG",
                height=data.shape[0],
                width=data.shape[1],
                count=3,
                dtype="uint8",
                crs=crs,
                transform=transform,
                compress="deflate",
            ) as dst:
                for i in range(3):  # Write R, G, B bands
                    dst.write(rgb_data[..., i], i + 1)

            print(f"Saved COG: {output_path}")

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
            output_colorized = os.path.join(output_dir, f"{raster}_colorized.tif")

            # Save the color table to the output directory
            with open(color_table_file, "w") as f:
                f.write("\n".join(color_table_lines))

            print(f"Color table saved as {color_table_file}")

            # Use the text file with gdaldem
            subprocess.run(
                ["gdaldem", "color-relief", subdataset_name, color_table_file, output_colorized],
                check=True
            )

            print(f"Colorized raster saved as {output_colorized}")

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
