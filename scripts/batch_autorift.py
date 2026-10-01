import os
import json
import glob
import shutil
import argparse
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

def apply_optical_filtering(image_path, nodata=0):
    """
    Mirrors the exact native single-image workflow of hyp3-autorift:
    - Evaluates Landsat 4/5 via the native FFT filter to remove whisk-broom striping.
    - Evaluates all other optical sensors (Landsat 7/8/9, Sentinel-2) via the Wallis filter to normalize illumination.
    - Writes the result out as a Float32 GeoTIFF.
    """
    filename = os.path.basename(image_path)
    
    # Check if the image is from an older whisk-broom sensor
    needs_fft = "Landsat4" in filename or "Landsat5" in filename
        
    print(f"  -> Opening {filename} to apply optical pre-filtering...")
    
    ds = gdal.Open(image_path)
    band = ds.GetRasterBand(1)
    original_array = band.ReadAsArray()
    
    if needs_fft:
        print(f"  -> Applying native FFT Filter for whisk-broom striping (NoData: {nodata})")
        filter_output = apply_fft_filter(original_array.copy(), nodata)
    else:
        # Default to Wallis for everything else (Sentinel-2, Landsat 7/8/9)
        print(f"  -> Applying native Wallis Filter for illumination normalization (NoData: {nodata})")
        filter_output = apply_wallis_nodata_fill_filter(original_array.copy(), nodata)
        
    if isinstance(filter_output, tuple):
        filtered_array = filter_output[0]
    else:
        filtered_array = filter_output

    # --- NATIVE FFT ROUTING CHECK ---
    # If the FFT filter aborted/skipped, mirror single-image mode: discard the array and use the raw file.
    if needs_fft and filtered_array.max() <= 3.0 and original_array.max() > 3.0:
        print("  -> Native FFT threshold not met. Bypassing filter and using raw input scene.")
        ds = None
        return image_path

    # Create the folder and write a clean Float32 file
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
    # Use gdal.OpenEx with IGNORE_COG_LAYOUT_BREAK to allow NoData edits on COGs
    ds = gdal.OpenEx(filepath, gdal.OF_UPDATE, open_options=["IGNORE_COG_LAYOUT_BREAK=YES"])
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
        band = ds.GetRasterBand(i) # Read band i to ensure all 4 bands are distinct
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

def snap_to_global_grid(input_path, output_filename, target_grid_m=90):
    """
    Forces the input image to snap its bounding box to exact multiples of the target grid.
    This mathematically guarantees pixel alignment across all independent autoRIFT runs within the same UTM zone.
    """
    print(f"  -> Snapping {os.path.basename(input_path)} to global {target_grid_m}m grid...")
    ds = gdal.Open(input_path)
    gt = ds.GetGeoTransform()
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    
    # Native pixel size (e.g., 10m for Sentinel-2, 30m for Landsat)
    x_res = abs(gt[1])
    y_res = abs(gt[5])
    
    minx = gt[0]
    maxy = gt[3]
    maxx = minx + (cols * gt[1])
    miny = maxy + (rows * gt[5])
    
    # Snap bounds to perfect multiples of the target grid (90m)
    import math
    new_minx = math.floor(minx / target_grid_m) * target_grid_m
    new_miny = math.floor(miny / target_grid_m) * target_grid_m
    new_maxx = math.ceil(maxx / target_grid_m) * target_grid_m
    new_maxy = math.ceil(maxy / target_grid_m) * target_grid_m
    
    gdal.Warp(
        output_filename,
        input_path,
        outputBounds=(new_minx, new_miny, new_maxx, new_maxy),
        xRes=x_res,
        yRes=y_res,
        resampleAlg=gdal.GRA_Bilinear
    )
    ds = None
    return output_filename

def mosaic_and_average(input_files, output_file):
    """
    Mosaics multiple overlapping GeoTIFFs by calculating the mean of valid pixels.
    This prevents hard seamline artifacts at the boundaries of adjacent satellite paths.
    """
    if not input_files:
        return
        
    print(f"  -> Building global grid for {len(input_files)} files...")
    # 1. Build a temporary VRT just to calculate the global bounding box and metadata
    vrt_path = "temp_mosaic_bounds.vrt"
    vrt_ds = gdal.BuildVRT(vrt_path, input_files)
    
    cols = vrt_ds.RasterXSize
    rows = vrt_ds.RasterYSize
    bands = vrt_ds.RasterCount
    geo_transform = vrt_ds.GetGeoTransform()
    projection = vrt_ds.GetProjection()
    
    # 2. Prepare accumulators for the average calculation
    sum_array = np.zeros((bands, rows, cols), dtype=np.float32)
    count_array = np.zeros((bands, rows, cols), dtype=np.float32)
    
    # 3. Add each file to the accumulators
    for file in input_files:
        print(f"  -> Merging {os.path.basename(file)} into average...")
        ds = gdal.Open(file)
        gt = ds.GetGeoTransform()
        
        # Calculate pixel offsets within the global mosaic grid
        x_offset = int(round((gt[0] - geo_transform[0]) / geo_transform[1]))
        y_offset = int(round((gt[3] - geo_transform[3]) / geo_transform[5]))
        
        for b in range(1, bands + 1):
            band = ds.GetRasterBand(b)
            arr = band.ReadAsArray().astype(np.float32)
            nodata = band.GetNoDataValue()
            
            # Mask out nodata and NaNs
            if nodata is not None:
                valid_mask = (arr != nodata) & (~np.isnan(arr))
            else:
                valid_mask = (arr != 0) & (~np.isnan(arr))
            
            # Slice the target region and add to accumulators where valid
            target_sum = sum_array[b-1, y_offset:y_offset+ds.RasterYSize, x_offset:x_offset+ds.RasterXSize]
            target_count = count_array[b-1, y_offset:y_offset+ds.RasterYSize, x_offset:x_offset+ds.RasterXSize]
            
            target_sum[valid_mask] += arr[valid_mask]
            target_count[valid_mask] += 1
            
        ds = None
        
    # 4. Calculate average and handle division by zero
    print("  -> Calculating mean across overlapping pixels...")
    with np.errstate(divide='ignore', invalid='ignore'):
        avg_array = np.divide(sum_array, count_array)
        avg_array[count_array == 0] = 0  # Restore NoData where no valid pixels exist
        
    # 5. Write the final blended output
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_file, cols, rows, bands, gdal.GDT_Float32, options=["COMPRESS=LZW"])
    out_ds.SetGeoTransform(geo_transform)
    out_ds.SetProjection(projection)
    
    for b in range(1, bands + 1):
        out_band = out_ds.GetRasterBand(b)
        out_band.WriteArray(avg_array[b-1])
        out_band.SetNoDataValue(0) # standard autoRIFT NoData
        
    out_ds.FlushCache()
    out_ds = None
    vrt_ds = None
    if os.path.exists(vrt_path):
        os.remove(vrt_path)

def process_event(manifest_path, enable_filtering=False):
    # Read the manifest generated by coseis.py
    with open(manifest_path, 'r') as f:
        manifest = json.load(f)
        
    if manifest.get("status") != "DOWNLOADED_READY_FOR_AUTORIFT":
        print(f"Skipping {manifest.get('event_title')} - Status is {manifest.get('status')}")
        return

    manifest["status"] = "PROCESSING_AUTORIFT"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=4)

    title = manifest["event_title"]
    sensor_str = manifest.get("sensor", "optical")
    level_str = manifest.get("optical_level", "unknown")
    
    # Check for track_pairs (new format), path_exports (old format), or single composite
    pairs_to_process = []
    
    if "track_pairs" in manifest:
        for path_num, paths in manifest["track_pairs"].items():
            pairs_to_process.append({
                "id": f"Path{path_num}",
                "ref": paths.get("pre_image"),
                "sec": paths.get("post_image")
            })
    elif "path_exports" in manifest:
        for path_num, paths in manifest["path_exports"].items():
            pairs_to_process.append({
                "id": f"Path{path_num}",
                "ref": paths.get("pre_composite_path"),
                "sec": paths.get("post_composite_path")
            })
    else:
        # Fallback for single-pair manifests
        pairs_to_process.append({
            "id": "Composite",
            "ref": manifest.get("pre_composite_path") or manifest.get("pre_image"),
            "sec": manifest.get("post_composite_path") or manifest.get("post_image")
        })
        
    # Clean the list of any pairs that failed to pull actual file paths
    pairs_to_process = [p for p in pairs_to_process if p["ref"] and p["sec"]]
        
    if not pairs_to_process:
        print("Error: Could not find valid image paths in manifest.")
        return

    print(f"\n{'='*50}")
    print(f"Starting autoRIFT processing for: {title} ({len(pairs_to_process)} regions)")
    print(f"{'='*50}")
    
    # Change working directory to the event folder to store outputs
    event_dir = os.path.dirname(manifest_path)
    original_dir = os.getcwd()
    os.chdir(event_dir)
    
    metric_offset_files = []
    velocity_files = []
    
    try:
        parameter_file = '/vsicurl/https://its-live-data.s3.amazonaws.com/autorift_parameters/v001/autorift_solidearth_0120m.shp'
        
        # Loop through each individual path/composite pair
        for pair in pairs_to_process:
            pair_id = pair["id"]
            
            # Use absolute paths so we can safely change working directories
            ref_path = os.path.abspath(pair["ref"])
            sec_path = os.path.abspath(pair["sec"])
            
            print(f"\n--- Processing Pair: {pair_id} ---")
            
            # Create an isolated working directory for this track to prevent 
            # autoRIFT's hardcoded intermediate files from cross-contaminating
            pair_work_dir = os.path.join(event_dir, f"processing_{pair_id}_{sensor_str}_{level_str}")
            os.makedirs(pair_work_dir, exist_ok=True)
            
            # Move into the isolated directory
            os.chdir(pair_work_dir)
            
            # Define universal physical parameters in meters early so we can use it for snapping
            TARGET_GRID_M = 90
            TARGET_CHIP_MIN_M = 360
            TARGET_CHIP_MAX_M = 720
            TARGET_SEARCH_LIMIT_M = 45

            # 1. Snap raw inputs to a universal global grid
            print("Snapping composite boundaries to a universal global grid...")
            aligned_ref_path = f"aligned_ref_{pair_id}.tif"
            aligned_sec_path = f"aligned_sec_{pair_id}.tif"
            
            snap_to_global_grid(ref_path, aligned_ref_path, target_grid_m=TARGET_GRID_M)
            snap_to_global_grid(sec_path, aligned_sec_path, target_grid_m=TARGET_GRID_M)

            # Assign nodata to 0, if missing from metadata
            print("Verifying NoData metadata headers...")
            assign_nodata_if_missing(aligned_ref_path, nodata_val=0)
            assign_nodata_if_missing(aligned_sec_path, nodata_val=0)
            
            if enable_filtering:
                print("Applying optical pre-filtering...")
                filtered_ref = apply_optical_filtering(aligned_ref_path, nodata=0)
                filtered_sec = apply_optical_filtering(aligned_sec_path, nodata=0)
            else:
                print("Bypassing optical pre-filtering...")
                filtered_ref = aligned_ref_path
                filtered_sec = aligned_sec_path

            print("Calculating bounding box from local image...")
            info_json = gdal.Info(filtered_ref, format='json')
            
            # --- 2. DYNAMIC PHYSICAL-TO-PIXEL CONVERSION ---
            # Detect native pixel size (e.g., 15m for Landsat 8, 10m for Sentinel-2)
            pixel_size = abs(info_json['geoTransform'][1])
            print(f"Detected native pixel size: {pixel_size}m")
            
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
            
            g_info['dem'] = resample_to_grid(g_info.get('dem'), 'local_dem.tif', TARGET_GRID_M, bbox=scene_bbox)
            g_info['dhdx'] = resample_to_grid(g_info.get('dhdx'), 'local_dhdx.tif', TARGET_GRID_M, bbox=scene_bbox)
            g_info['dhdy'] = resample_to_grid(g_info.get('dhdy'), 'local_dhdy.tif', TARGET_GRID_M, bbox=scene_bbox)
            g_info['vx'] = resample_to_grid(g_info.get('vx'), 'local_vx.tif', TARGET_GRID_M, bbox=scene_bbox)
            g_info['vy'] = resample_to_grid(g_info.get('vy'), 'local_vy.tif', TARGET_GRID_M, bbox=scene_bbox)
            g_info['ssm'] = resample_to_grid(g_info.get('ssm'), 'local_ssm.tif', TARGET_GRID_M, bbox=scene_bbox, is_mask=True)

            print("Manually co-registering...")
            obj = GeogridOptical()
            x1a, y1a, xsize1, ysize1, x2a, y2a, xsize2, ysize2, trans = obj.coregister(filtered_ref, filtered_sec)

            info_m = Dummy()
            info_m.startingX = trans[0]
            info_m.startingY = trans[3]
            info_m.XSize = trans[1]
            info_m.YSize = trans[5]
            info_m.numberOfLines = ysize1
            info_m.numberOfSamples = xsize1
            info_m.filename = filtered_ref 
            info_m.time = "20190101" 

            info_s = Dummy()
            info_s.time = "20200101" 

            print("Running Geogrid...")
            geogrid_info = runGeogrid(info_m, info_s, epsg=parameter_info['epsg'], optical_flag=1, **parameter_info['geogrid'])

            print("Running autoRIFT...")
            gdal.AllRegister()

            try:
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

            print("Renaming and moving output files...")
            out_prefix = f"{title}_{pair_id}_{sensor_str}_{level_str}_autorift"
            
            if os.path.exists("velocity.tif"):
                vel_path = os.path.join(event_dir, f"{out_prefix}_velocity.tif")
                os.rename("velocity.tif", vel_path)
                velocity_files.append(vel_path)
                
            if os.path.exists("offset.tif"):
                metric_offset_path = os.path.join(event_dir, f"{out_prefix}_offset_m.tif")
                convert_offset_to_meters("offset.tif", metric_offset_path, pixel_size)
                
                os.rename("offset.tif", os.path.join(event_dir, f"{out_prefix}_offset.tif"))
                metric_offset_files.append(metric_offset_path)
                
            # Step back up to the main event directory for the next iteration
            os.chdir(event_dir)
                
        # ==========================================
        # Post-Processing: Mosaicing Step
        # ==========================================
        
        # Mosaic offsets using pixel averaging to avoid seamline artifacts
        if metric_offset_files:
            final_mosaic_path = f"{title}_final_mosaic_{sensor_str}_{level_str}_offset_m.tif"
            if len(metric_offset_files) > 1:
                print(f"\n--- Blending {len(metric_offset_files)} paths into seamless offset mosaic ---")
                mosaic_and_average(metric_offset_files, final_mosaic_path)
                print(f"Created seamless offset mosaic: {final_mosaic_path}")
            else:
                shutil.copy(metric_offset_files[0], final_mosaic_path)
                print(f"Only one path found. Copied to final output: {final_mosaic_path}")

        # Mosaic velocities using pixel averaging
        if velocity_files:
            final_vel_path = f"{title}_final_mosaic_{sensor_str}_{level_str}_velocity.tif"
            if len(velocity_files) > 1:
                print(f"--- Blending {len(velocity_files)} paths into seamless velocity mosaic ---")
                mosaic_and_average(velocity_files, final_vel_path)
                print(f"Created seamless velocity mosaic: {final_vel_path}")
            else:
                shutil.copy(velocity_files[0], final_vel_path)

        # Update the manifest to show processing is complete
        manifest["status"] = "PROCESSED_AUTORIFT_COMPLETE"
        with open(os.path.basename(manifest_path), 'w') as f:
            json.dump(manifest, f, indent=4)
            
        print(f"\nFinished {title}. Outputs saved in {event_dir}")

    finally:
        # ALWAYS revert back to the root directory so the next loop iteration works
        os.chdir(original_dir)

def main():
    parser = argparse.ArgumentParser(description="Run autoRIFT batch processing on GEE optical downloads.")
    parser.add_argument("--filter", action="store_true", help="Enable FFT/Wallis pre-filtering for optical imagery.")
    parser.add_argument("--data_dir", type=Path, default=Path(__file__).resolve().parent / "data",
                        help="coseis.py data directory containing GEE_Optical_Downloads/. "
                             "Default: scripts/data (where coseis.py writes when run from scripts/).")
    args = parser.parse_args()

    # Start timer
    script_start_time = time.time()

    # Search for all manifest files in the GEE_Optical_Downloads directory
    search_pattern = args.data_dir.resolve() / "GEE_Optical_Downloads" / "**" / "*_autorift_manifest.json"

    manifest_files = glob.glob(str(search_pattern), recursive=True)
    
    if not manifest_files:
        print("No manifest files found. Ensure coseis.py has downloaded the composites.")
        return

    # Pass the argument directly into the loop
    for manifest_path in manifest_files:
        process_event(manifest_path, enable_filtering=args.filter)

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