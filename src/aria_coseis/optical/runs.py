"""
Optical autoRIFT runs on the data store: inputs, version definitions, records and deliverables.

A run of version V on event E reads the composites of one recipe from composites/E/, processes each
relative orbit ("path") in work/V/E/P<path>/ and writes to products/V/E/:

    E_V_offsets.tif                  the deliverable: mosaic of all paths (see write_deliverable)
    paths/E_V_P<path>_offset_m.tif   per path: dx_m, dy_m, interp_mask, chip_size_x_px
    run_info.json                    inputs (with checksums), per-path results, timings, status

The autoRIFT processing itself is in aria_coseis.optical.autorift; aria_coseis.optical.remosaic
makes a version from another version's per-path files with a different mosaic. This module needs
neither hyp3_autorift nor GDAL (only write_deliverable imports GDAL), so it can be tested anywhere.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from aria_coseis.optical import mosaic
from aria_coseis.optical.autorift_profiles import AutoriftProfile, PixelGeometry

DELIVERABLE_BANDS = {"mean": ("dx_m", "dy_m"), "aligned_weighted": ("dx_m", "dy_m", "sigma_h_m")}
DEFAULT_MOSAIC = "aligned_weighted"
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


def deliverable_spec(method: str, far_field_km: float = 0.0) -> dict[str, Any]:
    """
    The deliverable part of a version definition for a mosaic method, with the far-field
    correction if far_field_km > 0 (versions made before it existed have no such entry).
    """
    spec: dict[str, Any] = {
        "file": "<event>/<event>_<version>_offsets.tif",
        "bands": list(DELIVERABLE_BANDS[method]),
        "units": "m",
        "nodata": "NaN",
        "sign_convention": SIGN_CONVENTION,
    }
    if method == "mean":  # as recorded by the versions made before other methods existed
        spec["mosaic"] = "mean of the paths with a value at each node"
        return spec
    spec["mosaic"] = {
        "method": method,
        "description": (
            "constant per-path offsets (dx and dy separately) fitted by least squares to the "
            "median differences on clean overlap nodes and removed, relative to the path with "
            "the most valid nodes; then the mean of the paths weighted by 1 / local variance"
        ),
        "local_window_nodes": mosaic.WINDOW,
        "sigma_floor_m": mosaic.SIGMA_FLOOR_M,
        "min_overlap_nodes": mosaic.MIN_OVERLAP_NODES,
        "clean_overlap": (
            "both paths valid, not interpolated, each path's smallest chip, local scatter at most "
            "its 75th percentile; differences trimmed at 3 NMAD"
        ),
    }
    spec["band_descriptions"] = {
        "sigma_h_m": (
            "estimated 1-sigma horizontal uncertainty sqrt(sigma_dx^2 + sigma_dy^2), each from the "
            "local scatter of the contributing paths; does not include biases shared by all paths "
            "or smooth composite (snow, sun) errors"
        )
    }
    if far_field_km > 0:
        spec["far_field_correction"] = {
            "distance_km": far_field_km,
            "statistic": "median dx and dy of the nodes beyond the distance from the ruptures, subtracted",
            "ruptures": (
                "catalogued events with the epicentre in the product and the origin within its "
                "composite windows: USGS finite-fault sub-faults with >= 10 % of the peak slip, or a "
                "circle of half the Wells & Coppersmith (1994) surface-rupture length"
            ),
        }
    return spec


def version_definition(
    profile: AutoriftProfile,
    recipe: str,
    geometry: PixelGeometry,
    options: dict[str, Any],
    provenance: dict[str, Any],
    method: str = DEFAULT_MOSAIC,
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
        "deliverable": deliverable_spec(method),
    }


def derived_definition(
    source: dict[str, Any], method: str, provenance: dict[str, Any], far_field_km: float = 0.0
) -> dict[str, Any]:
    """
    Definition of a version made from the per-path files of the source version (its record)
    with another mosaic: the source definition with a new deliverable, the source version and
    the code that makes the mosaic. "code" stays that of the source, which made the paths.
    """
    if source.get("version") == "v0" or "definition" not in source:
        raise ValueError(f"cannot remosaic {source.get('version')}: not made by a store run")
    if "derived_from" in source["definition"]:
        original = source["definition"]["derived_from"]["version"]
        raise ValueError(f"{source['version']} is remosaicked itself; remosaic {original}")
    return {
        **source["definition"],
        "deliverable": deliverable_spec(method, far_field_km),
        "derived_from": {
            "version": source["version"],
            "uses": "the source version's per-path products, unchanged (hard links)",
        },
        "mosaic_code": {"aria_coseis": code_signature(provenance)["aria_coseis"]},
    }


def ensemble_definition(
    sources: list[dict[str, Any]], method: str, provenance: dict[str, Any], far_field_km: float
) -> dict[str, Any]:
    """
    Definition of a version mosaicked from the per-path products of several versions (their
    records), e.g. the same event from two composite recipes: every (version, path) is a member of
    the mosaic. The members must share the profile, so that they measure the same thing.
    """
    for src in sources:
        if (
            src.get("version") == "v0"
            or "definition" not in src
            or "derived_from" in src["definition"]
        ):
            raise ValueError(f"cannot use {src.get('version')}: not made by a store run")
    profiles = {json.dumps(s["definition"]["profile"], sort_keys=True) for s in sources}
    if len(profiles) != 1:
        raise ValueError("ensemble members must share one profile")
    if method != "aligned_weighted":
        raise ValueError("an ensemble needs the aligned_weighted mosaic (members are aligned)")
    spec = deliverable_spec(method, far_field_km)
    spec["mosaic"] = {
        **spec["mosaic"],
        "members": "every (version, path) of the source versions; aligned and weighted as paths are",
    }
    return {
        "ensemble_of": {s["version"]: s["definition"] for s in sources},
        "deliverable": spec,
        "derived_from": {
            "versions": [s["version"] for s in sources],
            "uses": "the source versions' per-path products, unchanged (hard links)",
        },
        "mosaic_code": {"aria_coseis": code_signature(provenance)["aria_coseis"]},
    }


def source_versions(definition: dict[str, Any]) -> list[str]:
    """The versions a derived version is made from."""
    d = definition["derived_from"]
    return list(d["versions"]) if "versions" in d else [d["version"]]


def mosaic_method(definition: dict[str, Any]) -> str:
    """The mosaic method of a version (versions made before methods existed use the mean)."""
    spec = definition["deliverable"]["mosaic"]
    return spec["method"] if isinstance(spec, dict) else "mean"


def catalog_row(
    event: str, event_meta: dict[str, Any], run: dict[str, Any], aria_commit: str
) -> dict[str, Any]:
    """One event's row in products/<version>/catalog.csv."""
    used = sorted({str(p["path"]) for p in run.get("paths", []) if p.get("offset_m")}, key=int)
    skipped = {**run.get("skipped", {})}
    for p in run.get("paths", []):
        if p.get("skipped"):
            skipped[str(p["path"])] = p["skipped"]
    return {
        "event": event,
        "event_id": event_meta.get("event_id") or "",
        "magnitude": event_meta.get("magnitude", ""),
        "status": run.get("status", ""),
        "paths_used": " ".join(used),
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
    done = {(r.get("version"), r["path"]) for r in run["paths"] if r.get("offset_m")}
    used = [p for p in run["inputs"] if (p.get("version"), p["path"]) in done]
    members = definition.get("ensemble_of")
    if members:  # an ensemble: the members' recipes and code, the shared profile
        first = next(iter(members.values()))
        recipe = " + ".join(m["composite_recipe"] for m in members.values())
        code = {
            mod: {v: m["code"].get(mod) for v, m in members.items()}
            for mod in ("aria_coseis", "hyp3_autorift", "autoRIFT")
        }
        profile = first["profile"]
    else:
        recipe, code, profile = (
            definition["composite_recipe"],
            definition["code"],
            definition["profile"],
        )
    meta = {
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
        "COMPOSITE_RECIPE": recipe,
        "PATHS": " ".join(sorted({str(p["path"]) for p in used}, key=int)),
        "INPUT_COMPOSITES": json.dumps(
            [
                {"name": p[f"{k}_name"], "sha256": p[f"{k}_sha256"]}
                for p in used
                for k in ("ref", "sec")
            ]
        ),
        "PROFILE": json.dumps(profile, sort_keys=True),
        "ARIA_COSEIS": json.dumps(code.get("aria_coseis"), sort_keys=True),
        "HYP3_AUTORIFT": json.dumps(code.get("hyp3_autorift"), sort_keys=True),
        "AUTORIFT": json.dumps(code.get("autoRIFT"), sort_keys=True),
        "CREATED_UTC": utc_now(),
    }
    if "derived_from" in definition:
        meta["DERIVED_FROM_VERSION"] = " ".join(source_versions(definition))
        meta["MOSAIC_CODE"] = json.dumps(definition["mosaic_code"]["aria_coseis"], sort_keys=True)
    return meta


def write_deliverable(
    path_files: dict[Any, str], out_path: str, metadata: dict[str, str], method: str = "mean"
) -> dict[str, Any]:
    """
    Mosaic per-path offset_m files ({path or member label: file}) into the deliverable GeoTIFF with
    aria_coseis.optical.mosaic: bands DELIVERABLE_BANDS[method], NaN nodata, band names and the
    given metadata, plus the per-path offsets removed. All paths must share one grid (same
    projection and pixel size, origins on the same lattice). A node is valid in a path only where
    both dx and dy are.
    """
    import numpy as np
    from osgeo import gdal

    gdal.UseExceptions()
    paths = sorted(path_files)
    files = [path_files[p] for p in paths]
    datasets = [gdal.Open(f) for f in files]
    proj = datasets[0].GetProjection()
    res = datasets[0].GetGeoTransform()[1]
    xmin = min(ds.GetGeoTransform()[0] for ds in datasets)
    ymax = max(ds.GetGeoTransform()[3] for ds in datasets)
    xmax = max(ds.GetGeoTransform()[0] + ds.RasterXSize * res for ds in datasets)
    ymin = min(ds.GetGeoTransform()[3] - ds.RasterYSize * res for ds in datasets)
    cols, rows = int(round((xmax - xmin) / res)), int(round((ymax - ymin) / res))
    stacks = np.full((4, len(files), rows, cols), np.nan, np.float32)  # dx, dy, interp, chip
    for k, (name, ds) in enumerate(zip(files, datasets)):
        gt = ds.GetGeoTransform()
        if ds.GetProjection() != proj or gt[1] != res or gt[5] != -res:
            raise ValueError(f"{name} is not on the grid of {files[0]}")
        c0, r0 = (gt[0] - xmin) / res, (ymax - gt[3]) / res
        if abs(c0 - round(c0)) > 1e-6 or abs(r0 - round(r0)) > 1e-6:
            raise ValueError(f"{name} is offset from the grid of {files[0]}")
        window = (
            slice(int(round(r0)), int(round(r0)) + ds.RasterYSize),
            slice(int(round(c0)), int(round(c0)) + ds.RasterXSize),
        )
        for b in range(4):
            stacks[b, k][window] = ds.GetRasterBand(b + 1).ReadAsArray()
    invalid = ~(np.isfinite(stacks[0]) & np.isfinite(stacks[1]))
    stacks[0][invalid] = np.nan
    stacks[1][invalid] = np.nan
    merged = mosaic.merge(*stacks, method=method)

    bands = DELIVERABLE_BANDS[method]
    arrays = {"dx_m": merged["dx"], "dy_m": merged["dy"], "sigma_h_m": merged.get("sigma_h")}
    metadata = {**metadata, "MOSAIC_METHOD": method}
    alignment = None
    if "alignment" in merged:
        alignment = {
            comp: {
                "reference_path": paths[a["reference_path_index"]],
                "offsets_removed_m": {str(p): o for p, o in zip(paths, a["offsets_m"])},
                "overlaps": [
                    {**o, "paths": [paths[i] for i in o["paths"]]} for o in a["overlaps"]
                ],
            }
            for comp, a in merged["alignment"].items()
        }
        metadata["PATH_OFFSETS_REMOVED_M"] = json.dumps(
            {comp: a["offsets_removed_m"] for comp, a in alignment.items()}, sort_keys=True
        )
        metadata["REFERENCE_PATH"] = str(alignment["dx"]["reference_path"])
    out = gdal.GetDriverByName("GTiff").Create(
        out_path,
        cols,
        rows,
        len(bands),
        gdal.GDT_Float32,
        options=["COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"],
    )
    out.SetGeoTransform((xmin, res, 0.0, ymax, 0.0, -res))
    out.SetProjection(proj)
    out.SetMetadata(metadata)
    for i, name in enumerate(bands):
        band = out.GetRasterBand(i + 1)
        band.WriteArray(arrays[name])
        band.SetDescription(name)
        band.SetNoDataValue(float("nan"))
        band.SetUnitType("m")
    out = None
    stats: dict[str, Any] = {
        "nodes": [rows, cols],
        "valid_nodes": int(np.isfinite(merged["dx"]).sum()),
        "paths": len(files),
    }
    if alignment is not None:
        stats["alignment"] = alignment
    return stats


def write_event_deliverable(
    version_record: dict[str, Any],
    event: str,
    event_meta: dict[str, Any],
    run: dict[str, Any],
    out_dir: Path,
    catalog: dict[str, Any] | None = None,
) -> bool:
    """
    Write the event's deliverable from the run's per-path products with the version's mosaic
    method, recording it in the run, then apply the version's far-field correction (which needs
    the catalogue of events, the store's events.json). Returns False when no path has a product.
    """
    files = {
        p.get("member", p["path"]): str(out_dir / p["offset_m"])
        for p in run["paths"]
        if p.get("offset_m")
    }
    if not files:
        return False
    run["deliverable"] = deliverable_name(event, version_record["version"])
    metadata = deliverable_metadata(version_record, event, event_meta, run)
    out_path = out_dir / run["deliverable"]
    run["deliverable_grid"] = write_deliverable(
        files, str(out_path), metadata, mosaic_method(version_record["definition"])
    )
    far = version_record["definition"]["deliverable"].get("far_field_correction")
    if far:
        if catalog is None:
            raise ValueError("the far-field correction needs the event catalogue")
        run["far_field"] = apply_far_field(out_path, run, catalog, far["distance_km"], event)
    return True


def apply_far_field(
    path: Path, run: dict[str, Any], catalog: dict[str, Any], distance_km: float, event: str
) -> dict[str, Any]:
    """
    Subtract the far-field median (aria_coseis.optical.farfield) from the deliverable's dx and dy
    in place, record it in the file's metadata, and save the ruptures used next to it.
    """
    import numpy as np
    from osgeo import gdal, osr
    from pyproj import Transformer

    from aria_coseis.optical import farfield

    gdal.UseExceptions()
    ds = gdal.Open(str(path), gdal.GA_Update)
    gt = ds.GetGeoTransform()
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    epsg = int(srs.GetAuthorityCode(None))
    x = gt[0] + (np.arange(ds.RasterXSize) + 0.5) * gt[1]
    y = gt[3] + (np.arange(ds.RasterYSize) + 0.5) * gt[5]
    to_ll = Transformer.from_crs(epsg, 4326, always_xy=True)
    corners = [to_ll.transform(cx, cy) for cx in x[[0, -1]] for cy in y[[0, -1]]]
    bounds = (
        min(c[0] for c in corners),
        min(c[1] for c in corners),
        max(c[0] for c in corners),
        max(c[1] for c in corners),
    )
    t0 = min(p["pre_window"][0] for p in run["inputs"])
    t1 = max(p["post_window"][1] for p in run["inputs"])
    quakes = farfield.events_in(catalog, bounds, t0, t1)
    if not any(q["folder"] == event for q in quakes) and event in catalog:
        quakes += farfield.events_in({event: catalog[event]}, (-180, -90, 180, 90), "0", "9")
    geom, records = farfield.rupture_geometry(quakes, epsg)
    dx = ds.GetRasterBand(1).ReadAsArray().astype(np.float64)
    dy = ds.GetRasterBand(2).ReadAsArray().astype(np.float64)
    result = farfield.far_field_shift(dx, dy, farfield.distance_to(geom, x, y), distance_km * 1000)
    rupture_file = path.with_name(path.name.replace("_offsets.tif", "_ruptures.geojson"))
    rupture_file.write_text(json.dumps(farfield.to_geojson(geom, epsg, records), indent=1) + "\n")
    record = {
        "distance_km": distance_km,
        **result,
        "ruptures": records,
        "rupture_file": rupture_file.name,
    }
    if result["shift_m"] is not None:
        for band, shift in zip((1, 2), result["shift_m"]):
            a = ds.GetRasterBand(band).ReadAsArray()
            ds.GetRasterBand(band).WriteArray((a - shift).astype(np.float32))
        ds.SetMetadataItem("FAR_FIELD_SHIFT_REMOVED_M", json.dumps(result["shift_m"]))
    else:
        ds.SetMetadataItem("FAR_FIELD_SHIFT_REMOVED_M", "none (too few far-field nodes)")
    ds.SetMetadataItem("FAR_FIELD_DISTANCE_KM", str(distance_km))
    ds.SetMetadataItem(
        "RUPTURES", json.dumps([{k: r.get(k) for k in ("event_id", "magnitude")} for r in records])
    )
    ds = None
    return record
