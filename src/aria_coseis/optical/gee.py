"""Google Earth Engine composites: export, GCS download, merge and nodata handling."""
import os
import subprocess
import time


def export_gee_sentinel2_composite(aoi_polygon, start_date, end_date, title, stage, gcs_bucket, collection_id, optical_level, crs_epsg='EPSG:4326'):
    """
    Generates a cloud-free median composite in GEE and exports to Google Cloud Storage.
    :param aoi_polygon: Shapely Polygon representing the Area of Interest
    :param start_date: Start date for the image collection filter (YYYY-MM-DD)
    :param end_date: End date for the image collection filter (YYYY-MM-DD)
    :param title: Title of the earthquake event (used for file naming)
    :param stage: 'pre-event' or 'post-event' to indicate the timing of the composite
    :param gcs_bucket: Name of the Google Cloud Storage bucket to export the composite
    :param collection_id: Sentinel-2 collection ID (e.g., 'COPERNICUS/S2_HARMONIZED')
    :param optical_level: Optical processing level (e.g., TOA, SR) for file naming
    :param crs_epsg: EPSG code for the coordinate reference system to use in the export (default is 'EPSG:4326')
    :return: The GEE export task object
    """
    import ee

    # Convert Shapely polygon to GEE Geometry
    bounds = aoi_polygon.bounds
    ee_roi = ee.Geometry.Rectangle([bounds[0], bounds[1], bounds[2], bounds[3]])

    # Query the parameterized Harmonized Collection
    s2_col = ee.ImageCollection(collection_id) \
        .filterBounds(ee_roi) \
        .filterDate(start_date, end_date)

    # Apply Cloud Score+ Masking
    cs_plus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED')
    s2_linked = s2_col.linkCollection(cs_plus, ['cs_cdf'])

    def mask_clouds(img):
        mask = img.select('cs_cdf').gte(0.65)
        return img.updateMask(mask)

    s2_masked = s2_linked.map(mask_clouds)

    # Fetch distinct orbits intersecting the AOI to the local Python environment
    distinct_orbits = ee.List(s2_masked.aggregate_array('SENSING_ORBIT_NUMBER')).distinct().getInfo()
    
    if not distinct_orbits:
        print(f"  No Sentinel-2 data found for {stage} between {start_date} and {end_date}.")
        return {}, []

    print(f"  Found {len(distinct_orbits)} distinct SENSING_ORBIT_NUMBER(s): {distinct_orbits}")
    
    orbit_exports = {}
    all_unique_dates = []

    # Iterate over each distinct orbit
    for orbit in distinct_orbits:
        orbit_col = s2_masked.filter(ee.Filter.eq('SENSING_ORBIT_NUMBER', orbit))
        
        def get_date(img):
            return ee.Feature(None, {'date': img.date().format("yyyy-MM-dd'T'HH:mm:ss")})
        
        raw_dates = orbit_col.map(get_date).aggregate_array('date').getInfo()
        orbit_dates = sorted(list(set(raw_dates)))
        all_unique_dates.extend(orbit_dates)
        
        # Calculate median for JUST this orbit
        composite = orbit_col.median().select('B8').toUint16().clip(ee_roi)

        # Clean up the start and end dates for file naming
        clean_start = start_date.replace('T', '_').replace(':', '') + 'UTC' if 'T' in start_date else start_date + '_000000UTC'
        clean_end = end_date.replace('T', '_').replace(':', '') + 'UTC' if 'T' in end_date else end_date + '_000000UTC'

        # Create Export Task
        file_name = f"{title}_S2_Path{orbit}_{optical_level.upper()}_B8_{stage}_{clean_start}_to_{clean_end}"
        prefix = f"COSEIS_Composites/{title}/{file_name}"
        
        task = ee.batch.Export.image.toCloudStorage(
            image=composite,
            description=file_name[:100],
            bucket=gcs_bucket,
            fileNamePrefix=prefix,
            region=ee_roi,
            scale=10, 
            crs=crs_epsg,
            maxPixels=1e13,
            formatOptions={'cloudOptimized': False}
        )
        
        task.start()
        print(f"  Started GCS Export Task for Orbit {orbit}: {file_name}")
        
        orbit_exports[str(orbit)] = {
            'task': task,
            'gcs_uri': f"gs://{gcs_bucket}/{prefix}.tif",
            'prefix': prefix,
            'dates_included': orbit_dates
        }

    return orbit_exports, sorted(list(set(all_unique_dates)))


def export_gee_landsat_composite(aoi_polygon, start_date, end_date, title, stage, gcs_bucket, collection_id, band_name, scale, optical_level, crs_epsg='EPSG:4326'):
    """Generates a cloud-free median TOA composite for the explicitly provided Landsat mission.
    :param aoi_polygon: Shapely Polygon representing the Area of Interest
    :param start_date: Start date for the image collection filter (YYYY-MM-DD)
    :param end_date: End date for the image collection filter (YYYY-MM-DD)
    :param title: Title of the earthquake event (used for file naming)
    :param stage: 'pre-event' or 'post-event' to indicate the timing of the composite
    :param gcs_bucket: Name of the Google Cloud Storage bucket to export the composite
    :param collection_id: Landsat collection ID (e.g., 'LANDSAT/LC08/C02/T1_L2')
    :param band_name: Name of the optical band to export (e.g., 'SR_B4' for Landsat 8 Red)
    :param scale: Scale in meters for the export (e.g., 30 for Landsat)
    :param optical_level: Optical processing level (e.g., TOA, SR, RAW) for file naming
    :param crs_epsg: EPSG code for the coordinate reference system to use in the export (default is 'EPSG:4326')
    :return: The GEE export task object and a list of unique acquisition dates
    """
    import ee

    bounds = aoi_polygon.bounds
    ee_roi = ee.Geometry.Rectangle([bounds[0], bounds[1], bounds[2], bounds[3]])

    # Determine Landsat Mission from collection_id for dynamic file naming
    if 'LT05' in collection_id:
        mission_name = 'Landsat5'
    elif 'LE07' in collection_id:
        mission_name = 'Landsat7'
    elif 'LC08' in collection_id:
        mission_name = 'Landsat8'
    elif 'LC09' in collection_id:
        mission_name = 'Landsat9'
    else:
        mission_name = 'Landsat'

    # Query the explicit collection passed into the function
    l_col = ee.ImageCollection(collection_id) \
        .filterBounds(ee_roi) \
        .filterDate(start_date, end_date)

    # Apply USGS QA_PIXEL bitmask for clouds and shadows
    def mask_clouds(img):
        qa = img.select('QA_PIXEL')
        cloud_shadow_bitmask = (1 << 4)
        clouds_bitmask = (1 << 3)
        dilated_cloud_bitmask = (1 << 1)
        
        mask = qa.bitwiseAnd(cloud_shadow_bitmask).eq(0) \
            .And(qa.bitwiseAnd(clouds_bitmask).eq(0)) \
            .And(qa.bitwiseAnd(dilated_cloud_bitmask).eq(0))
            
        return img.updateMask(mask).select([band_name], ['OPTICAL_BAND'])

    l_masked = l_col.map(mask_clouds)

    # Fetch distinct paths intersecting the AOI to the local Python environment
    # .getInfo() pulls the list from GEE servers to local execution
    distinct_paths = ee.List(l_masked.aggregate_array('WRS_PATH')).distinct().getInfo()
    
    if not distinct_paths:
        print(f"  No Landsat data found for {stage} between {start_date} and {end_date}.")
        return {}, []

    print(f"  Found {len(distinct_paths)} distinct WRS_PATH(s): {distinct_paths}")
    
    path_exports = {}
    all_unique_dates = []

    # Iterate over each distinct path to create independent composites and export tasks
    for path in distinct_paths:
        path_col = l_masked.filter(ee.Filter.eq('WRS_PATH', path))
        
        # Map over the specific path collection to get acquisition dates
        def get_date(img):
            return ee.Feature(None, {'date': img.date().format("yyyy-MM-dd'T'HH:mm:ss")})
        
        raw_dates = path_col.map(get_date).aggregate_array('date').getInfo()
        path_dates = sorted(list(set(raw_dates)))
        all_unique_dates.extend(path_dates)
        
        # Calculate median for JUST this path to avoid cross-track blending
        composite_raw = path_col.median().select('OPTICAL_BAND').clip(ee_roi)
        
        # Scale based on the collection type to ensure autoRIFT gets clean uint16
        if 'TOA' in collection_id:
            composite = composite_raw.multiply(10000).toUint16()
        elif '_L2' in collection_id:
            # GEE C02 L2 Surface Reflectance scaling: (pixel * 0.0000275) - 0.2
            composite = composite_raw.multiply(0.0000275).subtract(0.2).multiply(10000).max(0).toUint16()
        else:
            composite = composite_raw.unmask(0).toUint16()
        
        # Clean up the start and end dates for file naming
        clean_start = start_date.replace('T', '_').replace(':', '') + 'UTC' if 'T' in start_date else start_date + '_000000UTC'
        clean_end = end_date.replace('T', '_').replace(':', '') + 'UTC' if 'T' in end_date else end_date + '_000000UTC'

        # Include the specific WRS_PATH and optical_level in the file name
        file_name = f"{title}_{mission_name}_{optical_level.upper()}_Path{path}_{band_name}_{stage}_{clean_start}_to_{clean_end}"
        prefix = f"COSEIS_Composites/{title}/{file_name}"
        
        task = ee.batch.Export.image.toCloudStorage(
            image=composite,
            description=file_name[:100],
            bucket=gcs_bucket,
            fileNamePrefix=prefix,
            region=ee_roi,
            scale=scale, 
            crs=crs_epsg,
            maxPixels=1e13,
            formatOptions={'cloudOptimized': False}
        )
        
        task.start()
        print(f"  Started GCS Export Task for Path {path}: {file_name}")
        
        # Store the task and metadata keyed by path string to build the manifest later
        path_exports[str(path)] = {
            'task': task,
            'gcs_uri': f"gs://{gcs_bucket}/{prefix}.tif",
            'prefix': prefix,
            'dates_included': path_dates
        }

    return path_exports, sorted(list(set(all_unique_dates)))


def wait_for_gee_tasks(tasks, timeout_mins=60):
    """
    Polls GEE until all provided tasks are either COMPLETED or FAILED.
    Includes a timeout to prevent infinite hangs on stuck GEE backend tasks.
    :param tasks: List of GEE export task objects to monitor
    :param timeout_mins: Maximum minutes to wait before canceling tasks.
    """
    print(f"Waiting for Google Earth Engine exports to complete (Timeout: {timeout_mins} mins)...")
    start_time = time.time()
    timeout_seconds = timeout_mins * 60

    for task in tasks:
        while task.active():
            elapsed = time.time() - start_time
            if elapsed > timeout_seconds:
                print(f"  Task {task.id} TIMED OUT after {timeout_mins} minutes. Canceling task.")
                try:
                    task.cancel()
                except Exception as e:
                    print(f"  Could not cancel task {task.id}: {e}")
                break
                
            print(f"  Task {task.id} is {task.status()['state']}... waiting 30 seconds.")
            time.sleep(30)
        
        status = task.status()
        if status['state'] == 'COMPLETED':
            print(f"  Task {task.id} COMPLETED.")
        elif status['state'] in ['CANCELLED', 'CANCELED']:
            print(f"  Task {task.id} CANCELLED due to timeout.")
        else:
            print(f"  Task {task.id} FAILED: {status.get('error_message', 'Unknown error')}")


def download_from_gcs(bucket_name, prefix, local_dir):
    """Downloads a blob from the GCS bucket to user's local machine.
    :param bucket_name: Name of the GCS bucket
    :param prefix: Prefix of the blob in the GCS bucket (including path) to identify the file(s) to download
    :param local_dir: Local directory where the file(s) should be downloaded
    """
    from google.cloud import storage

    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blobs = bucket.list_blobs(prefix=prefix)

    os.makedirs(local_dir, exist_ok=True)
    downloaded_files = []

    for blob in blobs:
        # Extract the actual filename GEE generated (including any chip suffixes)
        filename = os.path.basename(blob.name)
        dest_path = os.path.join(local_dir, filename)
        
        blob.download_to_filename(dest_path)
        print(f"  Successfully downloaded chip: {dest_path}")
        downloaded_files.append(dest_path)
        
    return downloaded_files


def merge_and_compress_chips(chip_paths, output_path):
    """Uses GDAL to merge GEE image chips into a single compressed GeoTIFF.
    :param chip_paths: List of file paths to the individual image chips downloaded from GCS
    :param output_path: Desired file path for the final merged and compressed GeoTIFF
    """
    print(f"  Stitching and compressing {len(chip_paths)} file(s)...")
    
    # Create temporary output name to prevent GDAL from overwriting input
    tmp_output = output_path + ".tmp.tif"
    
    # gdalwarp merges the chips and applies DEFLATE compression
    cmd = [
        "gdalwarp",
        "-co", "COMPRESS=DEFLATE",
        "-co", "PREDICTOR=2",
        "-co", "TILED=YES",
        "-co", "NUM_THREADS=ALL_CPUS",
        *chip_paths,
        tmp_output
    ]
    
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Delete the individual chips
    for chip in chip_paths:
        try:
            os.remove(chip)
        except OSError:
            pass
            
    # Rename the clean, temporary file to your final desired filename
    os.rename(tmp_output, output_path)
        
    print(f"  Successfully created unified composite: {output_path}")
    return output_path


def assign_nodata(filepath, nodata_val=0):
    """
    Opens the specified GeoTIFF and explicitly writes the NoData 
    value into the metadata header.
    """
    from osgeo import gdal

    # Open the file in update mode, explicitly passing the open option to break COG layout
    ds = gdal.OpenEx(filepath, gdal.OF_UPDATE, open_options=["IGNORE_COG_LAYOUT_BREAK=YES"])
    
    if ds is not None:
        band = ds.GetRasterBand(1)
        if band.GetNoDataValue() is None:
            band.SetNoDataValue(nodata_val)
        # FlushCache ensures the metadata is written securely to disk
        ds.FlushCache()
        ds = None
    else:
        print(f"Warning: Could not open {filepath} to assign NoData value.")
