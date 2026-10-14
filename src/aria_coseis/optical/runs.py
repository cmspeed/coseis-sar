"""
Optical autoRIFT runs on the data store: inputs, version definitions, records and deliverables.

A run of version V on event E reads the composites of one recipe from composites/E/, processes each
relative orbit ("path") in work/V/E/P<path>/ and writes to products/V/E/:

    E_V_offsets.tif                  the deliverable: mosaic of all paths (see write_deliverable)
    paths/E_V_P<path>_offset_m.tif   per path: dx_m, dy_m, interp_mask, chip_size_x_px
    run_info.json                    inputs (with checksums), per-path results, timings, status

The autoRIFT processing itself is in aria_coseis.optical.autorift. This module needs neither
hyp3_autorift nor GDAL (only write_deliverable imports GDAL), so it can be tested anywhere.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from aria_coseis.optical.autorift_profiles import AutoriftProfile, PixelGeometry

DELIVERABLE_BANDS = ("dx_m", "dy_m")
PATH_PRODUCT_BANDS = ("dx_m", "dy_m", "interp_mask", "chip_size_x_px")
SIGN_CONVENTION = (
    "dx positive east, dy positive north (map axes of the file's UTM projection): motion of the "
    "ground from the PRE composite to the POST composite"
)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def deliverable_name(event: str, version: str) -> str:
    return f"{event}_{version}_offsets.tif"


def path_product_name(event: str, version: str, path: int) -> str:
    return f"{event}_{version}_P{path}_offset_m.tif"


def composite_pairs(
    root: Path, event: str, recipe: str
) -> tuple[list[dict[str, Any]], dict[str, str]]:
    """
    PRE/POST composite pairs of one recipe for an event, one per relative orbit.

    Returns (pairs sorted by path, {path: reason} for paths that cannot be paired). Each pair has
    the path, both files' store paths, names, checksums and time windows.
    """
    index_path = root / "composites" / event / "composites.json"
    if not index_path.exists():
        return [], {}
    index = json.loads(index_path.read_text())
    by_path: dict[int, dict[str, list[tuple[str, dict[str, Any]]]]] = {}
    for name, entry in index.items():
        if name.startswith("_") or entry.get("recipe") != recipe:
            continue
        by_path.setdefault(int(entry["path"]), {"PRE": [], "POST": []})[entry["stage"]].append(
            (name, entry)
        )
    pairs, skipped = [], {}
    for path in sorted(by_path):
        stages = by_path[path]
        if len(stages["PRE"]) != 1 or len(stages["POST"]) != 1:
            skipped[str(path)] = (
                f"{len(stages['PRE'])} PRE and {len(stages['POST'])} POST composites"
            )
            continue
        (pre_name, pre), (post_name, post) = stages["PRE"][0], stages["POST"][0]
        pairs.append(
            {
                "path": path,
                "ref": str(root / "composites" / event / pre_name),
                "sec": str(root / "composites" / event / post_name),
                "ref_name": pre_name,
                "sec_name": post_name,
                "ref_sha256": pre["sha256"],
                "sec_sha256": post["sha256"],
                "pre_window": [pre["start"], pre["end"]],
                "post_window": [post["start"], post["end"]],
            }
        )
    return pairs, skipped


def code_signature(provenance: dict[str, Any]) -> dict[str, Any]:
    """The parts of store.code_provenance() that identify the code (not where it is installed)."""
    out = {}
    for module in ("aria_coseis", "hyp3_autorift", "autoRIFT"):
        info = provenance.get(module, {})
        if info.get("repo"):
            out[module] = {k: info[k] for k in ("commit", "dirty", "diff_sha256") if k in info}
        else:
            keys = ("installed", "distribution_version", "conda_package")
            out[module] = {k: info.get(k) for k in keys}
    return out


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def version_definition(
    profile: AutoriftProfile,
    recipe: str,
    geometry: PixelGeometry,
    options: dict[str, Any],
    provenance: dict[str, Any],
) -> dict[str, Any]:
    """What defines a product version; every run of the version must have the same definition."""
    return {
        "composite_recipe": recipe,
        "profile": _jsonable(profile.parameters()),
        "sentinel2_pixel_geometry": {
            "grid_px": geometry.grid_px,
            "chip_sizes_px": list(geometry.chip_sizes_px),
            "chip_to_grid_ratio": geometry.ratio,
            "effective_search_px": geometry.effective_search_px,
        },
        "driver_options": _jsonable(options),
        "code": code_signature(provenance),
        "deliverable": {
            "file": "<event>/<event>_<version>_offsets.tif",
            "bands": list(DELIVERABLE_BANDS),
            "units": "m",
            "nodata": "NaN",
            "sign_convention": SIGN_CONVENTION,
            "mosaic": "mean of the paths with a value at each node",
        },
    }


def catalog_row(
    event: str, event_meta: dict[str, Any], run: dict[str, Any], aria_commit: str
) -> dict[str, Any]:
    """One event's row in products/<version>/catalog.csv."""
    used = [str(p["path"]) for p in run.get("paths", []) if p.get("offset_m")]
    skipped = {**run.get("skipped", {})}
    for p in run.get("paths", []):
        if p.get("skipped"):
            skipped[str(p["path"])] = p["skipped"]
    return {
        "event": event,
        "event_id": event_meta.get("event_id") or "",
        "magnitude": event_meta.get("magnitude", ""),
        "status": run.get("status", ""),
        "paths_used": " ".join(sorted(used, key=int)),
        "paths_skipped": "; ".join(
            f"{p}: {r}" for p, r in sorted(skipped.items(), key=lambda kv: int(kv[0]))
        ),
        "processed_utc": run.get("finished_utc", ""),
        "aria_code": aria_commit[:7],
        "mosaic": f"{event}/{run['deliverable']}" if run.get("deliverable") else "",
    }


def deliverable_metadata(
    version_record: dict[str, Any], event: str, event_meta: dict[str, Any], run: dict[str, Any]
) -> dict[str, str]:
    """GDAL metadata items written into the deliverable GeoTIFF."""
    definition = version_record["definition"]
    done = {r["path"] for r in run["paths"] if r.get("offset_m")}
    used = [p for p in run["inputs"] if p["path"] in done]
    code = definition["code"]
    return {
        "COSEIS_PRODUCT": "Sentinel-2 optical surface offsets (autoRIFT feature tracking)",
        "COSEIS_VERSION": version_record["version"],
        "EVENT_FOLDER": event,
        "USGS_EVENT_ID": str(event_meta.get("event_id") or ""),
        "USGS_TITLE": str(event_meta.get("title") or ""),
        "EVENT_TIME_UTC": str(event_meta.get("time_utc") or ""),
        "EPICENTER_LON_LAT_DEPTH_KM": json.dumps(
            [event_meta.get("longitude"), event_meta.get("latitude"), event_meta.get("depth_km")]
        ),
        "UNITS": "m",
        "SIGN_CONVENTION": SIGN_CONVENTION,
        "COMPOSITE_RECIPE": definition["composite_recipe"],
        "PATHS": " ".join(str(p["path"]) for p in used),
        "INPUT_COMPOSITES": json.dumps(
            [
                {"name": p[f"{k}_name"], "sha256": p[f"{k}_sha256"]}
                for p in used
                for k in ("ref", "sec")
            ]
        ),
        "PROFILE": json.dumps(definition["profile"], sort_keys=True),
        "ARIA_COSEIS": json.dumps(code.get("aria_coseis"), sort_keys=True),
        "HYP3_AUTORIFT": json.dumps(code.get("hyp3_autorift"), sort_keys=True),
        "AUTORIFT": json.dumps(code.get("autoRIFT"), sort_keys=True),
        "CREATED_UTC": utc_now(),
    }


def write_deliverable(
    path_files: list[str], out_path: str, metadata: dict[str, str]
) -> dict[str, Any]:
    """
    Mosaic per-path offset_m files into the deliverable GeoTIFF: band 1 dx_m, band 2 dy_m (mean of
    the paths with a value at each node), NaN nodata, band names and the given metadata. All paths
    must share one grid (same projection and pixel size, origins on the same lattice).
    """
    import numpy as np
    from osgeo import gdal

    gdal.UseExceptions()
    datasets = [gdal.Open(p) for p in path_files]
    proj = datasets[0].GetProjection()
    gt0 = datasets[0].GetGeoTransform()
    res = gt0[1]
    xmin = min(ds.GetGeoTransform()[0] for ds in datasets)
    ymax = max(ds.GetGeoTransform()[3] for ds in datasets)
    xmax = max(ds.GetGeoTransform()[0] + ds.RasterXSize * res for ds in datasets)
    ymin = min(ds.GetGeoTransform()[3] - ds.RasterYSize * res for ds in datasets)
    cols, rows = int(round((xmax - xmin) / res)), int(round((ymax - ymin) / res))
    sums = np.zeros((2, rows, cols), np.float64)
    count = np.zeros((rows, cols), np.int16)
    for path, ds in zip(path_files, datasets):
        gt = ds.GetGeoTransform()
        if ds.GetProjection() != proj or gt[1] != res or gt[5] != -res:
            raise ValueError(f"{path} is not on the grid of {path_files[0]}")
        c0, r0 = (gt[0] - xmin) / res, (ymax - gt[3]) / res
        if abs(c0 - round(c0)) > 1e-6 or abs(r0 - round(r0)) > 1e-6:
            raise ValueError(f"{path} is offset from the grid of {path_files[0]}")
        window = (
            slice(int(round(r0)), int(round(r0)) + ds.RasterYSize),
            slice(int(round(c0)), int(round(c0)) + ds.RasterXSize),
        )
        dx = ds.GetRasterBand(1).ReadAsArray()
        dy = ds.GetRasterBand(2).ReadAsArray()
        valid = np.isfinite(dx) & np.isfinite(dy)
        sums[0][window][valid] += dx[valid]
        sums[1][window][valid] += dy[valid]
        count[window][valid] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = (sums / count).astype(np.float32)
    mean[:, count == 0] = np.nan

    driver = gdal.GetDriverByName("GTiff")
    out = driver.Create(
        out_path,
        cols,
        rows,
        len(DELIVERABLE_BANDS),
        gdal.GDT_Float32,
        options=["COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"],
    )
    out.SetGeoTransform((xmin, res, 0.0, ymax, 0.0, -res))
    out.SetProjection(proj)
    out.SetMetadata(metadata)
    for i, name in enumerate(DELIVERABLE_BANDS):
        band = out.GetRasterBand(i + 1)
        band.WriteArray(mean[i])
        band.SetDescription(name)
        band.SetNoDataValue(float("nan"))
        band.SetUnitType("m")
    out = None
    valid = count > 0
    return {"nodes": [rows, cols], "valid_nodes": int(valid.sum()), "paths": len(path_files)}
