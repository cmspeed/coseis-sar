"""
autoRIFT on the optical data store: `python -m aria_coseis.optical.autorift`.

    python -m aria_coseis.optical.autorift --version exp-g90 --profile g90 --event <event> [...]

Reads the composites of one recipe from the store (aria_coseis.optical.store), processes each
relative orbit in work/<version>/<event>/, and writes products/<version>/<event>/ (see
aria_coseis.optical.runs). What is made is set by a processing profile
(aria_coseis.optical.autorift_profiles), optionally with overrides; the first run of a version
records its profile, recipe and code in VERSION.json, and later runs must match it.
"""

from __future__ import annotations

import argparse
import inspect
import json
import math
import os
import shutil
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np
from hyp3_autorift import geometry, utils
from hyp3_autorift.process import apply_fft_filter, apply_wallis_nodata_fill_filter
from hyp3_autorift.vend.testautoRIFT import generateAutoriftProduct
from hyp3_autorift.vend.testGeogrid import GeogridOptical, runGeogrid
from osgeo import gdal, osr

from aria_coseis.optical import runs, store
from aria_coseis.optical.autorift_profiles import (
    PROFILES,
    AutoriftProfile,
    aligned_bounds,
    clipped_chip_sizes,
    driver_options,
    pixel_geometry,
    resolve_profile,
)

if not hasattr(np.lib, "pad"):  # autoRIFT 2.0.0 calls np.lib.pad (gone in numpy 2)
    np.lib.pad = np.pad

ITS_LIVE_PARAMETER_FILE = (
    "/vsicurl/https://its-live-data.s3.amazonaws.com/autorift_parameters/v001/"
    "autorift_solidearth_0120m.shp"
)
# Geogrid needs acquisition dates; one year apart makes its "velocity" output numerically equal
# to the displacement in meters (the composites span months, not single dates).
DUMMY_REFERENCE_DATE = "20190101"
DUMMY_SECONDARY_DATE = "20200101"


class Dummy:
    pass


def check_driver() -> None:
    """Stop unless hyp3_autorift is the fork's driver (branch coseis), which these runs need."""
    if "CorrectOpticalDy0Sign" not in inspect.getsource(generateAutoriftProduct):
        raise SystemExit(
            "hyp3_autorift is not the COSEIS driver: use cmspeed/hyp3-autorift branch coseis "
            "(e.g. PYTHONPATH=<checkout>/src)"
        )


def apply_optical_filtering(
    image_path: str, nodata: float = 0, write_zero_mask: bool = False
) -> str:
    """
    Wallis-fill filter an image (FFT filter for Landsat 4/5), as hyp3-autorift does for Landsat 7;
    writes filtered/<name>_filtered.tif next to it and returns that path.

    With write_zero_mask, the Wallis filter's mask of large nodata areas is written next to the
    filtered image as <name>_zeroMask.tif, where hyp3-autorift's driver looks for it.
    """
    filename = os.path.basename(image_path)
    needs_fft = "Landsat4" in filename or "Landsat5" in filename
    print(f"  -> Opening {filename} to apply optical pre-filtering...")
    ds = gdal.Open(image_path)
    original_array = ds.GetRasterBand(1).ReadAsArray()
    if needs_fft:
        print(f"  -> Applying native FFT Filter for whisk-broom striping (NoData: {nodata})")
        filter_output = apply_fft_filter(original_array.copy(), nodata)
    else:
        print(
            f"  -> Applying native Wallis Filter for illumination normalization (NoData: {nodata})"
        )
        filter_output = apply_wallis_nodata_fill_filter(original_array.copy(), nodata)

    zero_mask = None
    if isinstance(filter_output, tuple):
        filtered_array, zero_mask = filter_output
    else:
        filtered_array = filter_output

    # If the FFT filter skipped, use the raw image, as hyp3-autorift does
    if needs_fft and filtered_array.max() <= 3.0 and original_array.max() > 3.0:
        print("  -> Native FFT threshold not met. Bypassing filter and using raw input scene.")
        return image_path

    filtered_dir = os.path.join(os.path.dirname(os.path.abspath(image_path)), "filtered")
    os.makedirs(filtered_dir, exist_ok=True)
    filtered_path = os.path.join(filtered_dir, filename.replace(".tif", "_filtered.tif"))
    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(filtered_path, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32)
    out_ds.SetGeoTransform(ds.GetGeoTransform())
    out_ds.SetProjection(ds.GetProjection())
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(filtered_array)
    out_band.SetNoDataValue(nodata)
    out_ds = None
    if write_zero_mask and zero_mask is not None:
        root, ext = os.path.splitext(filtered_path)
        mask_ds = driver.Create(
            f"{root}_zeroMask{ext}", ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte
        )
        mask_ds.SetGeoTransform(ds.GetGeoTransform())
        mask_ds.SetProjection(ds.GetProjection())
        mask_ds.GetRasterBand(1).WriteArray(zero_mask.astype(np.uint8))
        mask_ds = None
    print(f"  -> Filter applied successfully. Saved to: {filtered_path}")
    return filtered_path


def assign_nodata_if_missing(filepath: str, nodata_val: float = 0) -> None:
    """Declare nodata_val as the nodata value of a GeoTIFF that has none."""
    ds = gdal.OpenEx(filepath, gdal.OF_UPDATE, open_options=["IGNORE_COG_LAYOUT_BREAK=YES"])
    band = ds.GetRasterBand(1)
    if band.GetNoDataValue() is None:
        print(f"  -> NoData tag missing in header. Injecting NoData={nodata_val}...")
        band.SetNoDataValue(nodata_val)
    ds = None


def raster_epsg(path: str) -> int:
    """Return the EPSG code of a raster's coordinate system."""
    srs = osr.SpatialReference(wkt=gdal.Open(path).GetProjection())
    srs.AutoIdentifyEPSG()
    return int(srs.GetAuthorityCode(None))


def raster_has_data(path: str, nodata: float = 0) -> bool:
    """True if any pixel differs from nodata (checked on a 4x-decimated read)."""
    ds = gdal.Open(path)
    arr = ds.GetRasterBand(1).ReadAsArray(
        buf_xsize=max(1, ds.RasterXSize // 4), buf_ysize=max(1, ds.RasterYSize // 4)
    )
    return bool(np.any(arr != nodata))


def resample_to_grid(
    input_path: str | None,
    output_filename: str,
    resolution: float,
    bbox: tuple[float, float, float, float] | None = None,
    is_mask: bool = False,
    dst_epsg: int | None = None,
    pad_m: float = 5000,
) -> str | None:
    """
    Crop (scene box plus pad_m) and resample a parameter-tile raster to the target resolution.

    With dst_epsg, the raster is also reprojected to that coordinate system and bbox is read in
    it. The parameter tiles come in their own projection (a UTM zone, or polar stereographic at
    high latitudes), which can differ from the composite's (the epicenter's UTM zone).
    Vector rasters (vx, vy, dhdx, dhdy) are resampled but their components are not rotated.
    """
    if not input_path:
        return None
    print(f"  -> Cropping & Resampling {os.path.basename(input_path)} to {resolution}m...")
    alg = gdal.GRA_NearestNeighbour if is_mask else gdal.GRA_Bilinear
    target_bounds = None
    if bbox:
        minx, miny, maxx, maxy = bbox
        target_bounds = (minx - pad_m, miny - pad_m, maxx + pad_m, maxy + pad_m)
    gdal.Warp(
        output_filename,
        input_path,
        xRes=resolution,
        yRes=resolution,
        outputBounds=target_bounds,
        resampleAlg=alg,
        dstSRS=f"EPSG:{dst_epsg}" if dst_epsg else None,
    )
    return output_filename


def write_flat_geogrid_inputs(image_path: str, grid_m: float, out_dir: str) -> dict[str, str]:
    """
    Zero DEM, slopes and reference velocity on a grid of grid_m aligned to multiples of grid_m.

    For map-projected optical images in the same CRS, Geogrid uses only the DEM's *grid*: node
    positions come from it and the offset-to-velocity conversion does not depend on the slopes.
    The DEM gets nodata -32767 because Geogrid compares velocities with the DEM's nodata value,
    which GDAL reports as 0 when unset (a zero velocity would then count as missing).
    """
    ds = gdal.Open(image_path)
    gt = ds.GetGeoTransform()
    minx, maxy = gt[0], gt[3]
    maxx = minx + ds.RasterXSize * gt[1]
    miny = maxy + ds.RasterYSize * gt[5]
    x0, y0, x1, y1 = aligned_bounds(minx, miny, maxx, maxy, grid_m)
    ncol, nrow = int(round((x1 - x0) / grid_m)), int(round((y1 - y0) / grid_m))
    driver = gdal.GetDriverByName("GTiff")
    paths = {}
    for name in ("dem", "dhdx", "dhdy", "vx", "vy"):
        path = os.path.join(out_dir, f"flat_{name}.tif")
        out = driver.Create(path, ncol, nrow, 1, gdal.GDT_Float32)
        out.SetGeoTransform((x0, grid_m, 0.0, y1, 0.0, -grid_m))
        out.SetProjection(ds.GetProjection())
        band = out.GetRasterBand(1)
        if name == "dem":
            band.SetNoDataValue(-32767)
        band.Fill(0)
        out = None
        paths[name] = path
    print(f"  -> Flat Geogrid inputs: {ncol} x {nrow} nodes at {grid_m:g} m, origin ({x0}, {y1})")
    return paths


def convert_offset_to_meters(
    offset_tif_path: str, output_tif_path: str, pixel_size_m: float
) -> None:
    """
    offset.tif (pixels) to offset_m: bands dx_m (east +), dy_m (north +), interp_mask (1 = filled
    by autoRIFT's interpolation), chip_size_x_px (chip that succeeded, image pixels; 0 = none).
    Invalid nodes are NaN, declared as nodata.
    """
    ds = gdal.Open(offset_tif_path)
    out_ds = gdal.GetDriverByName("GTiff").Create(
        output_tif_path,
        ds.RasterXSize,
        ds.RasterYSize,
        ds.RasterCount,
        gdal.GDT_Float32,
        options=["COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"],
    )
    out_ds.SetGeoTransform(ds.GetGeoTransform())
    out_ds.SetProjection(ds.GetProjection())
    for i in range(1, ds.RasterCount + 1):
        arr = ds.GetRasterBand(i).ReadAsArray()
        out_band = out_ds.GetRasterBand(i)
        out_band.WriteArray(arr * pixel_size_m if i in (1, 2) else arr)
        if i <= len(runs.PATH_PRODUCT_BANDS):
            out_band.SetDescription(runs.PATH_PRODUCT_BANDS[i - 1])
        out_band.SetNoDataValue(float("nan"))
    out_ds = None


def snap_to_global_grid(input_path: str, output_filename: str, target_grid_m: float) -> str:
    """Pad an image's bounds outward to multiples of the grid (lossless for 10 m-aligned S2)."""
    print(f"  -> Snapping {os.path.basename(input_path)} to global {target_grid_m:g}m grid...")
    ds = gdal.Open(input_path)
    gt = ds.GetGeoTransform()
    minx, maxy = gt[0], gt[3]
    maxx = minx + ds.RasterXSize * gt[1]
    miny = maxy + ds.RasterYSize * gt[5]
    gdal.Warp(
        output_filename,
        input_path,
        outputBounds=(
            math.floor(minx / target_grid_m) * target_grid_m,
            math.floor(miny / target_grid_m) * target_grid_m,
            math.ceil(maxx / target_grid_m) * target_grid_m,
            math.ceil(maxy / target_grid_m) * target_grid_m,
        ),
        xRes=abs(gt[1]),
        yRes=abs(gt[5]),
        resampleAlg=gdal.GRA_Bilinear,
    )
    return output_filename


def its_live_geogrid_inputs(
    filtered_ref: str, profile: AutoriftProfile, scene_bbox: tuple, info_json: dict
) -> dict:
    """ITS_LIVE parameter-tile inputs, cropped and resampled to the profile grid."""
    coords = info_json["wgs84Extent"]["coordinates"][0]
    lons = [c[0] for c in coords]
    lats = [c[1] for c in coords]
    poly = geometry.polygon_from_bbox(
        x_limits=(min(lats), max(lats)), y_limits=(min(lons), max(lons))
    )
    print("Fetching regional parameters from S3 Shapefile...")
    parameter_info = utils.find_jpl_parameter_info(poly, ITS_LIVE_PARAMETER_FILE)

    # The parameter tile is chosen by the scene center and has its own projection (that UTM
    # zone, or polar stereographic at high latitudes); the composite is in the epicenter's UTM
    # zone. Bring the tile rasters into the image's projection so the crop below and Geogrid's
    # output grid use the image's coordinates.
    image_epsg = raster_epsg(filtered_ref)
    if parameter_info["epsg"] != image_epsg:
        print(
            f"Parameter tile {parameter_info['name']} is EPSG:{parameter_info['epsg']}; "
            f"reprojecting its rasters to the image's EPSG:{image_epsg}"
        )
        parameter_info["epsg"] = image_epsg

    g_info = parameter_info["geogrid"]
    for key in ("csminx", "csminy", "csmaxx", "csmaxy", "srx", "sry"):
        g_info[key] = None
    for key in ("chip_size_min", "chip_size_max", "search_range"):
        parameter_info["autorift"][key] = None

    # The grid is the resampled DEM's. v0 cropped the scene box plus 5 km (grid off the grid_m
    # lattice); aligned profiles put it on multiples of grid_m.
    if profile.align_grid:
        bounds, pad = aligned_bounds(*scene_bbox, profile.grid_m, pad_m=5000), 0
    else:
        bounds, pad = scene_bbox, 5000
    print(f"Cropping & Resampling shapefile maps to scene bounds ({profile.grid_m:g}m grid)...")
    for key, name, is_mask in (
        ("dem", "local_dem.tif", False),
        ("dhdx", "local_dhdx.tif", False),
        ("dhdy", "local_dhdy.tif", False),
        ("vx", "local_vx.tif", False),
        ("vy", "local_vy.tif", False),
        ("ssm", "local_ssm.tif", True),
    ):
        g_info[key] = resample_to_grid(
            g_info.get(key),
            name,
            profile.grid_m,
            bbox=bounds,
            is_mask=is_mask,
            dst_epsg=image_epsg,
            pad_m=pad,
        )
    return parameter_info


def flat_geogrid_inputs(filtered_ref: str, profile: AutoriftProfile, work_dir: str) -> dict:
    """Parameter-file-free Geogrid/autoRIFT inputs (see write_flat_geogrid_inputs)."""
    flat = write_flat_geogrid_inputs(filtered_ref, profile.grid_m, work_dir)
    return {
        "name": "flat",
        "epsg": raster_epsg(filtered_ref),
        "geogrid": {
            **{key: flat[key] for key in ("dem", "dhdx", "dhdy", "vx", "vy")},
            "srx": None,
            "sry": None,
            "csminx": None,
            "csminy": None,
            "csmaxx": None,
            "csmaxy": None,
            "ssm": None,
        },
        "autorift": {
            "grid_location": "window_location.tif",
            "init_offset": "window_offset.tif",
            "search_range": None,
            "chip_size_min": None,
            "chip_size_max": None,
            "offset2vx": "window_rdr_off2vel_x_vec.tif",
            "offset2vy": "window_rdr_off2vel_y_vec.tif",
            "stable_surface_mask": None,
            "scale_factor": "window_scale_factor.tif",
            "mpflag": 0,
        },
    }


def process_pair(
    pair: dict[str, Any],
    profile: AutoriftProfile,
    options: dict[str, Any],
    work_dir: Path,
    out_path: Path,
) -> dict[str, Any]:
    """
    Run Geogrid and autoRIFT on one path's PRE/POST composites in work_dir; write the path's
    offset_m product to out_path. Returns the path's record (outputs, grid, timings or 'skipped').
    """
    result: dict[str, Any] = {"path": pair["path"], "timings_s": {}}
    t_start = time.time()
    print(f"\n--- Path {pair['path']} ---")
    for image in (pair["ref"], pair["sec"]):
        if not raster_has_data(image):
            print(f"Skipping path {pair['path']}: {os.path.basename(image)} has no data")
            return {**result, "skipped": "empty composite"}

    if work_dir.exists():
        shutil.rmtree(work_dir)  # nothing from an earlier attempt may leak into this one
    work_dir.mkdir(parents=True)
    os.chdir(work_dir)

    aligned_ref = snap_to_global_grid(pair["ref"], "aligned_ref.tif", profile.grid_m)
    aligned_sec = snap_to_global_grid(pair["sec"], "aligned_sec.tif", profile.grid_m)
    assign_nodata_if_missing(aligned_ref, 0)
    assign_nodata_if_missing(aligned_sec, 0)

    t0 = time.time()
    if profile.prefilter in ("wallis_hps", "wallis"):
        write_mask = profile.prefilter == "wallis"
        filtered_ref = apply_optical_filtering(aligned_ref, 0, write_zero_mask=write_mask)
        filtered_sec = apply_optical_filtering(aligned_sec, 0, write_zero_mask=write_mask)
    else:
        filtered_ref, filtered_sec = aligned_ref, aligned_sec
    result["timings_s"]["prefilter"] = round(time.time() - t0, 1)

    info_json = gdal.Info(filtered_ref, format="json")
    gt = info_json["geoTransform"]
    pixel_size = abs(gt[1])
    geom = pixel_geometry(profile, pixel_size)
    if driver_options(profile, geom) != options:
        raise ValueError(
            f"{pixel_size:g} m pixels give other driver options than the version's (10 m S2)"
        )
    width, height = info_json["size"]
    minx, maxy = gt[0], gt[3]
    maxx, miny = minx + width * gt[1], maxy + height * gt[5]
    scene_bbox = (min(minx, maxx), min(miny, maxy), max(minx, maxx), max(miny, maxy))
    clipped = clipped_chip_sizes(geom, abs(maxx - minx), abs(maxy - miny))
    if clipped:
        print(f"  [Warning] image smaller than the coarse filter window for chips {clipped} px")
        result["warnings"] = [f"coarse filter window wider than the image for chips {clipped} px"]

    if profile.geogrid_inputs == "its_live":
        parameter_info = its_live_geogrid_inputs(filtered_ref, profile, scene_bbox, info_json)
    else:
        parameter_info = flat_geogrid_inputs(filtered_ref, profile, str(work_dir))
    parameter_info["geogrid"]["grid_spacing"] = profile.grid_m
    parameter_info["geogrid"]["chipSizeX0"] = profile.chip_min_m

    _, _, xsize1, ysize1, _, _, _, _, trans = GeogridOptical().coregister(
        filtered_ref, filtered_sec
    )
    info_m = Dummy()
    info_m.startingX, info_m.startingY = trans[0], trans[3]
    info_m.XSize, info_m.YSize = trans[1], trans[5]
    info_m.numberOfLines, info_m.numberOfSamples = ysize1, xsize1
    info_m.filename = filtered_ref
    info_m.time = DUMMY_REFERENCE_DATE
    info_s = Dummy()
    info_s.time = DUMMY_SECONDARY_DATE

    t0 = time.time()
    geogrid_info = runGeogrid(
        info_m, info_s, epsg=parameter_info["epsg"], optical_flag=1, **parameter_info["geogrid"]
    )
    result["timings_s"]["geogrid"] = round(time.time() - t0, 1)

    t0 = time.time()
    gdal.AllRegister()
    generateAutoriftProduct(
        filtered_ref,
        filtered_sec,
        nc_sensor=None,  # no netCDF packaging (needs the fork's driver, branch coseis)
        optical_flag=True,
        ncname=None,
        geogrid_run_info=geogrid_info,
        parameter_file=ITS_LIVE_PARAMETER_FILE.replace("/vsicurl/", ""),
        **parameter_info["autorift"],
        **options,
    )
    result["timings_s"]["autorift"] = round(time.time() - t0, 1)

    offset_ds = gdal.Open("offset.tif")
    dx = offset_ds.GetRasterBand(1).ReadAsArray()
    result["grid"] = {
        "nodes": [offset_ds.RasterYSize, offset_ds.RasterXSize],
        "geotransform": list(offset_ds.GetGeoTransform()),
        "valid_fraction": round(float(np.isfinite(dx).mean()), 4),
    }
    offset_ds = None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    convert_offset_to_meters("offset.tif", str(out_path), pixel_size)
    result["offset_m"] = f"paths/{out_path.name}"
    result["timings_s"]["total"] = round(time.time() - t_start, 1)
    return result


def process_event(
    root: Path,
    version_record: dict[str, Any],
    event: str,
    profile: AutoriftProfile,
    options: dict[str, Any],
    provenance: dict[str, Any],
    force: bool = False,
) -> str:
    """Process one event for a version; returns its status. Errors are recorded, not raised."""
    version = version_record["version"]
    recipe = version_record["definition"]["composite_recipe"]
    out_dir = root / "products" / version / event
    run_path = out_dir / "run_info.json"
    if run_path.exists() and not force:
        status = json.loads(run_path.read_text()).get("status")
        if status in ("COMPLETE", "NO_DATA"):
            print(f"Skipping {event}: {version} is {status} (use --force to redo)")
            return status
    if out_dir.exists():
        shutil.rmtree(out_dir)  # a rerun replaces the event's products entirely
    out_dir.mkdir(parents=True)

    events_path = root / "events.json"
    events = json.loads(events_path.read_text()) if events_path.exists() else {}
    event_meta = events.get(event, {})
    pairs, skipped = runs.composite_pairs(root, event, recipe)
    run: dict[str, Any] = {
        "event": event,
        "version": version,
        "status": "PROCESSING",
        "started_utc": runs.utc_now(),
        "event_metadata": event_meta,
        "inputs": [{k: v for k, v in p.items() if k not in ("ref", "sec")} for p in pairs],
        "skipped": skipped,
        "paths": [],
        "code": provenance,
    }

    def save() -> None:
        run_path.write_text(json.dumps(run, indent=1) + "\n")

    save()
    print(f"\n{'=' * 60}\n{event}: {len(pairs)} path(s), version {version}\n{'=' * 60}")
    cwd = os.getcwd()
    try:
        for pair in pairs:
            work_dir = root / "work" / version / event / f"P{pair['path']}"
            out_path = out_dir / "paths" / runs.path_product_name(event, version, pair["path"])
            run["paths"].append(process_pair(pair, profile, options, work_dir, out_path))
            os.chdir(cwd)
            save()
        files = [str(out_dir / p["offset_m"]) for p in run["paths"] if p.get("offset_m")]
        if files:
            run["deliverable"] = runs.deliverable_name(event, version)
            metadata = runs.deliverable_metadata(version_record, event, event_meta, run)
            run["deliverable_grid"] = runs.write_deliverable(
                files, str(out_dir / run["deliverable"]), metadata
            )
            run["status"] = "COMPLETE"
        else:
            run["status"] = "NO_DATA"
    except Exception as e:  # noqa: BLE001 - record the failure and go on with the next event
        run["status"] = "FAILED"
        run["error"] = f"{type(e).__name__}: {e}"
        run["traceback"] = traceback.format_exc()
        print(f"FAILED {event}: {run['error']}")
    finally:
        os.chdir(cwd)
    run["finished_utc"] = runs.utc_now()
    save()
    aria_commit = version_record["definition"]["code"].get("aria_coseis", {}).get("commit", "")
    store.update_catalog(
        root / "products" / version / "catalog.csv",
        runs.catalog_row(event, event_meta, run, aria_commit),
    )
    return run["status"]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="aria-coseis-autorift",
        description="Run autoRIFT on composites from the optical data store into a product version.",
    )
    parser.add_argument("--version", required=True, help="Product version, e.g. exp-g90, v1-90m")
    parser.add_argument("--description", default="", help="What the version is for (first run)")
    parser.add_argument("--store", help="Store root (default: COSEIS_OPTICAL_DIR)")
    parser.add_argument("--recipe", default="med90", help="Composite recipe (default med90)")
    parser.add_argument(
        "--event", action="append", help="Event folder name (repeatable); default: all events"
    )
    parser.add_argument("--profile", choices=sorted(PROFILES), default="g90")
    parser.add_argument("--grid_m", type=float, help="Output posting, m")
    parser.add_argument("--chip_min_m", type=float, help="Smallest chip, m")
    parser.add_argument("--chip_max_m", type=float, help="Largest chip, m (chip_min x 2^k)")
    parser.add_argument("--search_px", type=int, help="Search half-width, image pixels (min 6)")
    parser.add_argument("--oversample", type=int, help="Correlation oversampling for every chip")
    parser.add_argument("--prefilter", choices=["wallis_hps", "wallis", "hps"])
    parser.add_argument("--geogrid_inputs", choices=["its_live", "flat"])
    parser.add_argument("--reference_velocity", choices=["null", "its_live"])
    parser.add_argument("--force", action="store_true", help="Redo events already processed")
    args = parser.parse_args(argv)

    root = store.require_store(args.store).resolve()
    check_driver()
    null_ref = None if args.reference_velocity is None else args.reference_velocity == "null"
    profile = resolve_profile(
        args.profile,
        grid_m=args.grid_m,
        chip_min_m=args.chip_min_m,
        chip_max_m=args.chip_max_m,
        search_px=args.search_px,
        oversample=args.oversample,
        prefilter=args.prefilter,
        geogrid_inputs=args.geogrid_inputs,
        null_reference_velocity=null_ref,
    )
    geom = pixel_geometry(profile, 10.0)
    for warning in geom.warnings:
        print(f"[Warning] {warning}")
    options = driver_options(profile, geom)
    provenance = store.code_provenance()
    definition = runs.version_definition(profile, args.recipe, geom, options, provenance)
    record = store.start_version(root, args.version, args.description, definition)
    env_listing = root / "products" / args.version / "environment.txt"
    if not env_listing.exists():
        store.write_environment_listing(env_listing)
    print(
        f"Version {args.version}: profile {profile.name}, recipe {args.recipe}, chips "
        f"{geom.chip_sizes_px} px, grid {geom.grid_px} px, search {geom.effective_search_px} px"
    )

    events = args.event or sorted(
        p.parent.name for p in (root / "composites").glob("*/composites.json")
    )
    started = time.time()
    statuses = {
        event: process_event(root, record, event, profile, options, provenance, args.force)
        for event in events
    }
    summary = {s: list(statuses.values()).count(s) for s in sorted(set(statuses.values()))}
    print(f"\nDone in {(time.time() - started) / 60:.1f} min: {summary}")


if __name__ == "__main__":
    main()
