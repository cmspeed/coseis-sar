import os
import json
import glob
from pathlib import Path
from osgeo import gdal
import time
from hyp3_autorift import geometry, utils
from hyp3_autorift.vend.testGeogrid import GeogridOptical, runGeogrid
from hyp3_autorift.vend.testautoRIFT import generateAutoriftProduct
from hyp3_autorift.process import apply_fft_filter, apply_wallis_nodata_fill_filter

import numpy as np
if not hasattr(np.lib, 'pad'):
    np.lib.pad = np.pad

# Create dummy class to hold metadata
class Dummy(object):
    pass

def apply_custom_landsat_filtering(image_path, nodata=0):
    """
    Mirrors the exact native single-image workflow of hyp3-autorift:
    - Evaluates the input scene via the native FFT filter.
    - If the filter is skipped, it discards the array and returns the original path untouched.
    - If the filter is applied, it writes the result out as a Float32 GeoTIFF.
    """
    filename = os.path.basename(image_path)
    
    needs_fft = "Landsat4" in filename or "Landsat5" in filename
    needs_wallis = "Landsat7" in filename or "Landsat8" in filename or "Landsat9" in filename
    
    if not (needs_fft or needs_wallis):
        return image_path
        
    print(f"  -> Opening {filename} to evaluate filter...")
    
    ds = gdal.Open(image_path)
    band = ds.GetRasterBand(1)
    original_array = band.ReadAsArray()
    
    if needs_fft:
        print(f"  -> Evaluation via native FFT Filter (NoData: {nodata})")
        filter_output = apply_fft_filter(original_array.copy(), nodata)
    elif needs_wallis:
        print(f"  -> Evaluation via native Wallis NoData Fill Filter (NoData: {nodata})")
        filter_output = apply_wallis_nodata_fill_filter(original_array.copy(), nodata)
        
    if isinstance(filter_output, tuple):
        filtered_array = filter_output[0]
    else:
        filtered_array = filter_output

    # --- NATIVE ROUTING CHECK ---
    # If the filter aborted/skipped, its max value will be capped at its internal 
    # Z-score limit (3.0), while your raw image has much higher values. 
    # In this case, mirror single-image mode: discard the array and use the raw file.
    if needs_fft and filtered_array.max() <= 3.0 and original_array.max() > 3.0:
        print("  -> Native threshold not met. Bypassing filter and using raw input scene.")
        ds = None
        return image_path

    # If the filter DID trigger, create the folder and write a clean Float32 file
    base_dir = os.path.dirname(os.path.abspath(image_path))
    filtered_dir = os.path.join(base_dir, 'filtered')
    os.makedirs(filtered_dir, exist_ok=True)
    
    filtered_path = os.path.join(filtered_dir, filename.replace('.tif', '_filtered.tif'))
    
    driver = gdal.GetDriverByName('GTiff')
    # Native hyp3 filtered outputs are written as Float32
    out_ds = driver.Create(filtered_path, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
    
    out_ds.SetGeoTransform(ds.GetGeoTransform())
    out_ds.SetProjection(ds.GetProjection())
    
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(filtered_array)
    out_band.SetNoDataValue(nodata)
    
    out_ds.FlushCache()
    out_ds = None
    ds = None
    
    print(f"  -> Filter applied successfully. Saved to: {filtered_path}")
    return filtered_path

def assign_nodata_if_missing(filepath, nodata_val=0):
    """
    Checks if a NoData value is defined in the GeoTIFF header.
    If it is missing, it injects the specified nodata_val into the metadata.
    """
    ds = gdal.Open(filepath, gdal.GA_Update)
    if ds is None:
        print(f"  -> Error: Could not open {filepath} to update NoData.")
        return
        
    band = ds.GetRasterBand(1)
    current_nodata = band.GetNoDataValue()
    
    if current_nodata is None:
        print(f"  -> NoData tag missing in header. Injecting NoData={nodata_val}...")
        band.SetNoDataValue(nodata_val)
    
    # Flush the cache to ensure the header is written to disk
    ds.FlushCache()
    ds = None

def resample_to_grid(input_path, output_filename, resolution, bbox=None, is_mask=False):
    """Downloads, clips to scene bbox, and resamples an S3 raster to target resolution."""
    if not input_path:
        return None
    print(f"  -> Cropping & Resampling {os.path.basename(input_path)} to {resolution}m...")
    alg = gdal.GRA_NearestNeighbour if is_mask else gdal.GRA_Bilinear
    
    # Bounding box format for gdal.Warp: (minX, minY, maxX, maxY)
    # Adding a 5km (5000m) buffer to ensure complete coverage over image corners
    if bbox:
        minx, miny, maxx, maxy = bbox
        target_bounds = (minx - 5000, miny - 5000, maxx + 5000, maxy + 5000)
    else:
        target_bounds = None

    gdal.Warp(
        output_filename, 
        input_path, 
        xRes=resolution, 
        yRes=resolution, 
        outputBounds=target_bounds, 
        resampleAlg=alg
    )
    return output_filename

def convert_offset_to_meters(offset_tif_path, output_tif_path, pixel_size_m):
    """
    Reads offset.tif (in pixels), multiplies DX and DY by the native 
    sensor resolution (in meters), and writes offset_m.tif.
    """
    if not os.path.exists(offset_tif_path):
        return
        
    print(f"Converting pixel offsets to meters (Native Pixel Size: {pixel_size_m}m)...")
    ds = gdal.Open(offset_tif_path, gdal.GA_ReadOnly)
    
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(
        output_tif_path, 
        ds.RasterXSize, 
        ds.RasterYSize, 
        ds.RasterCount, 
        gdal.GDT_Float32
    )
    out_ds.SetGeoTransform(ds.GetGeoTransform())
    out_ds.SetProjection(ds.GetProjection())
    
    for i in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(1)
        arr = band.ReadAsArray()
        
        # Multiply DX (Band 1) and DY (Band 2) by pixel size in meters
        # Keep Band 3 (InterpMask) and Band 4 (ChipSize) unchanged
        if i in (1, 2):
            out_arr = arr * pixel_size_m
        else:
            out_arr = arr
            
        out_band = out_ds.GetRasterBand(i)
        out_band.WriteArray(out_arr)
        out_band.FlushCache()
        
    out_ds = None
    ds = None
    print(f"Saved metric offset raster: {output_tif_path}")

def process_event(manifest_path, enable_filtering=False):
    # Read the manifest generated by coseis.py
    with open(manifest_path, 'r') as f:
        manifest = json.load(f)
        
    if manifest.get("status") != "DOWNLOADED_READY_FOR_AUTORIFT":
        print(f"Skipping {manifest.get('event_title')} - Status is {manifest.get('status')}")
        return

    title = manifest["event_title"]
    ref_path = manifest["pre_composite_path"]
    sec_path = manifest["post_composite_path"]
    
    print(f"\n{'='*50}")
    print(f"Starting autoRIFT processing for: {title}")
    print(f"{'='*50}")
    
    # Change working directory to the event folder to store outputs
    event_dir = os.path.dirname(manifest_path)
    original_dir = os.getcwd()
    os.chdir(event_dir)
    
    try:
        local_ref = os.path.basename(ref_path)
        local_sec = os.path.basename(sec_path)
        
        # Assign nodata to 0, if missing from metadata
        print("Verifying NoData metadata headers...")
        assign_nodata_if_missing(local_ref, nodata_val=0)
        assign_nodata_if_missing(local_sec, nodata_val=0)
        
        if enable_filtering:
            print("Applying Landsat pre-filtering (Filter isolation mode)...")
            filtered_ref = apply_custom_landsat_filtering(local_ref, nodata=0)
            filtered_sec = apply_custom_landsat_filtering(local_sec, nodata=0)
        else:
            print("Bypassing Landsat pre-filtering (Parameter isolation mode)...")
            filtered_ref = local_ref
            filtered_sec = local_sec
        
        parameter_file = '/vsicurl/https://its-live-data.s3.amazonaws.com/autorift_parameters/v001/autorift_solidearth_0120m.shp'

        print("Calculating bounding box from local image...")
        info_json = gdal.Info(local_ref, format='json')
        
        # --- 1. DYNAMIC PHYSICAL-TO-PIXEL CONVERSION ---
        # Detect native pixel size (e.g., 15m for Landsat 8, 10m for Sentinel-2, 30m for Landsat 5)
        pixel_size = abs(info_json['geoTransform'][1])
        print(f"Detected native pixel size: {pixel_size}m")
        
        # Define universal physical parameters in meters
        TARGET_GRID_M = 90
        TARGET_CHIP_MIN_M = 360
        TARGET_CHIP_MAX_M = 720
        TARGET_SEARCH_LIMIT_M = 45
        
        # Translate to pixel-based values based on sensor resolution
        grid_px = int(TARGET_GRID_M / pixel_size)
        chip_min_px = int(TARGET_CHIP_MIN_M / pixel_size)
        chip_max_px = int(TARGET_CHIP_MAX_M / pixel_size)
        search_limit_px = int(TARGET_SEARCH_LIMIT_M / pixel_size)
        
        # Safeguard: Ensure chip sizes perfectly divide by grid spacing for autoRIFT's C++ core
        if chip_min_px % grid_px != 0:
            raise ValueError(f"Math Error: Min chip ({chip_min_px}px) is not a multiple of Grid size ({grid_px}px)")
        
        print(f"Dynamic arrays: {chip_min_px}-{chip_max_px}px chips, {search_limit_px}px search limit")
        # -----------------------------------------------

        coords = info_json['wgs84Extent']['coordinates'][0]
        lons = [c[0] for c in coords]
        lats = [c[1] for c in coords]
        poly = geometry.polygon_from_bbox(x_limits=(min(lats), max(lats)), y_limits=(min(lons), max(lons)))

        print("Fetching regional parameters from S3 Shapefile...")
        parameter_info = utils.find_jpl_parameter_info(poly, parameter_file)

        # --- 2. INJECT PARAMETERS & NULLIFY SHAPEFILE ---
        # Force 90m output grid spacing for Geogrid
        parameter_info['geogrid']['grid_spacing'] = TARGET_GRID_M
        parameter_info['geogrid']['chipSizeX0'] = TARGET_CHIP_MIN_M

        # Strip Geogrid shapefile paths using the correct kwargs keys
        parameter_info['geogrid']['csminx'] = None
        parameter_info['geogrid']['csminy'] = None
        parameter_info['geogrid']['csmaxx'] = None
        parameter_info['geogrid']['csmaxy'] = None
        parameter_info['geogrid']['srx'] = None 
        parameter_info['geogrid']['sry'] = None

        # Strip autoRIFT shapefile paths using the correct kwargs keys
        parameter_info['autorift']['chip_size_min'] = None
        parameter_info['autorift']['chip_size_max'] = None
        parameter_info['autorift']['search_range'] = None
        
        # Inject dynamically calculated pixels for autoRIFT arrays
        parameter_info['autorift']['ChipSizeMinX'] = chip_min_px
        parameter_info['autorift']['ChipSizeMinY'] = chip_min_px
        parameter_info['autorift']['ChipSizeMaxX'] = chip_max_px
        parameter_info['autorift']['ChipSizeMaxY'] = chip_max_px
        parameter_info['autorift']['SearchLimitX'] = search_limit_px
        parameter_info['autorift']['SearchLimitY'] = search_limit_px
        # ------------------------------------------------

        # Force target grid by local resampling of 120m DEM from parameter file
        # Extract native projected bounding box (in UTM meters) to limit download extent
        gt = info_json['geoTransform']
        raster_xsize = info_json['size'][0]
        raster_ysize = info_json['size'][1]
        
        minx = gt[0]
        maxy = gt[3]
        maxx = minx + (raster_xsize * gt[1])
        miny = maxy + (raster_ysize * gt[5])
        scene_bbox = (min(minx, maxx), min(miny, maxy), max(minx, maxx), max(miny, maxy))

        print(f"Cropping & Resampling shapefile maps to scene bounds ({TARGET_GRID_M}m grid)...")
        g_info = parameter_info['geogrid']
        
        # Resample only the necessary spatial patch from S3 and preserve local TIFs for inspection
        g_info['dem'] = resample_to_grid(g_info.get('dem'), 'local_dem.tif', TARGET_GRID_M, bbox=scene_bbox)
        g_info['dhdx'] = resample_to_grid(g_info.get('dhdx'), 'local_dhdx.tif', TARGET_GRID_M, bbox=scene_bbox)
        g_info['dhdy'] = resample_to_grid(g_info.get('dhdy'), 'local_dhdy.tif', TARGET_GRID_M, bbox=scene_bbox)
        g_info['vx'] = resample_to_grid(g_info.get('vx'), 'local_vx.tif', TARGET_GRID_M, bbox=scene_bbox)
        g_info['vy'] = resample_to_grid(g_info.get('vy'), 'local_vy.tif', TARGET_GRID_M, bbox=scene_bbox)
        g_info['ssm'] = resample_to_grid(g_info.get('ssm'), 'local_ssm.tif', TARGET_GRID_M, bbox=scene_bbox, is_mask=True)
        # -------------------------------------------------

        print("Manually co-registering to avoid filename parsing...")
        obj = GeogridOptical()
        
        # Pass filtered paths to geogrid
        x1a, y1a, xsize1, ysize1, x2a, y2a, xsize2, ysize2, trans = obj.coregister(filtered_ref, filtered_sec)

        # Inject dummy metadata
        info_m = Dummy()
        info_m.startingX = trans[0]
        info_m.startingY = trans[3]
        info_m.XSize = trans[1]
        info_m.YSize = trans[5]
        info_m.numberOfLines = ysize1
        info_m.numberOfSamples = xsize1
        
        # Update filename to the filtered ref
        info_m.filename = filtered_ref 
        info_m.time = "20190101"  # Dummy start date

        info_s = Dummy()
        info_s.time = "20200101"  # EXACTLY one year later to fix search limit scaling

        print("Running Geogrid...")
        geogrid_info = runGeogrid(info_m, info_s, epsg=parameter_info['epsg'], optical_flag=1, **parameter_info['geogrid'])

        print("Running autoRIFT...")
        gdal.AllRegister()

        try:
            # Pass the filtered paths to autorift
            generateAutoriftProduct(
                filtered_ref,
                filtered_sec,
                nc_sensor="DUMMY_SENSOR", 
                optical_flag=True,
                ncname=None,
                geogrid_run_info=geogrid_info,
                **parameter_info['autorift'],
                parameter_file=parameter_file.replace('/vsicurl/', '')
            )
        except Exception as e:
            if "netCDF packaging not supported" in str(e):
                print("Successfully completed autoRIFT (skipped NetCDF packaging).")
            else:
                raise e

        # Rename the hardcoded outputs to match the event
        print("Renaming output files for data provenance...")
        out_prefix = f"{title}_autorift"
        
        if os.path.exists("velocity.tif"):
            os.rename("velocity.tif", f"{out_prefix}_velocity.tif")
            
        if os.path.exists("offset.tif"):
            # Generate the offset in meters before renaming the raw pixel offset file
            metric_offset_path = f"{out_prefix}_offset_m.tif"
            convert_offset_to_meters("offset.tif", metric_offset_path, pixel_size)
            
            # Rename raw pixel offset file
            os.rename("offset.tif", f"{out_prefix}_offset.tif")

        # Update the manifest to show processing is complete
        manifest["status"] = "PROCESSED_AUTORIFT_COMPLETE"
        with open(os.path.basename(manifest_path), 'w') as f:
            json.dump(manifest, f, indent=4)
            
        print(f"Finished {title}. Outputs saved in {event_dir}")

    finally:
        # ALWAYS revert back to the root directory so the next loop iteration works
        os.chdir(original_dir)

def main():

    # Start timer
    script_start_time = time.time()
    
    # Gets the absolute path to the project root (one level up from /scripts)
    root_dir = Path(__file__).resolve().parent.parent 
    
    # Search for all manifest files in the GEE_Optical_Downloads directory
    search_pattern = root_dir / "data" / "GEE_Optical_Downloads" / "**" / "*_autorift_manifest.json"
    manifest_files = glob.glob(str(search_pattern), recursive=True)
    
    if not manifest_files:
        print("No manifest files found. Ensure coseis.py has downloaded the composites.")
        return

    RUN_WITH_FILTERS = False

    for manifest_path in manifest_files:
        process_event(manifest_path, enable_filtering=RUN_WITH_FILTERS)

    # Stop the timer and calculate elapsed time at the very end
    script_end_time = time.time()
    total_seconds = script_end_time - script_start_time

    # Convert and print runtime
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    print(f"\nProcessing complete.")
    print(f"Total execution time: {int(hours)}h {int(minutes)}m {seconds:.2f}s")

if __name__ == "__main__":
    main()