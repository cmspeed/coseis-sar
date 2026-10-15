"""
Composites of a new recipe for events already in the optical data store:
`python -m aria_coseis.optical.composites`.

    python -m aria_coseis.optical.composites plan --recipe seas90 --testbed
    python -m aria_coseis.optical.composites build --recipe seas90 --event <event> [--event ...]

Historic mode (aria_coseis.pipeline) finds events in a date range and makes their med90
composites. This module instead takes events by name from the store (events.json: USGS id, origin
time, epicentre) and makes the composites of another recipe (store.COMPOSITE_RECIPES) on exactly
the grid of the event's reference composites (med90), for the same paths. Products made from two
recipes can then be compared node by node.

For each event:
1. The AOI is made as the pipeline makes it: the USGS finite-fault model buffered by 0.15 degrees
   (as a rectangle), or a 1 x 1 degree box around the epicentre. Its UTM bounds, snapped outward
   to 10 m, must equal the reference composites' grid; otherwise the event is refused (e.g. if
   USGS has replaced the finite-fault model since).
2. The recipe gives each stage's window, or the recipe whose composite it reuses.
3. Per path and stage: the Sentinel-2 L1C scenes of that relative orbit in the window, cloud-masked
   with Cloud Score+ (and snow-masked if the recipe says so), median of B8, clipped to the AOI, and
   exported to Cloud Storage on the reference grid. A reused stage is not exported: the other
   recipe's composite is registered again under this recipe (hard link).
4. Exports are downloaded, chips merged, nodata 0 tagged, the grid checked, and the file is
   registered in the store with its window and scene list.

`plan` only reads (USGS and Earth Engine queries). `build` can be rerun: registered composites are
skipped, and the ids of started exports are kept in work/composites/<recipe>/<event>/tasks.json,
so a rerun waits for them instead of starting new ones.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from aria_coseis.aoi import load_aoi_from_json, make_aoi
from aria_coseis.optical import runs, store
from aria_coseis.usgs import get_ffm_geojson_url

EE_PROJECT = "coseis-1"
CLOUD_SCORE_COLLECTION = "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED"
AOI_BUFFER_DEG = 0.15  # as aria_coseis.pipeline for optical
GRID_M = 10
STAGES = ("PRE", "POST")
TESTBED_EVENTS = (
    "m_71_ridgecrest_earthquake_sequence",
    "m_78_pazarcik_earthquake_kahramanmaras_earthquake_sequence",
    "m_75_72_km_n_of_palu_indonesia",
    "m_77_2025_mandalay_burma_myanmar_earthquake",
)


def window_bound(origin: datetime, spec: list[Any]) -> datetime:
    """A window bound from a recipe: ["time", d] = origin + d days; ["date", d] = that day 00:00."""
    kind, days = spec
    moment = origin + timedelta(days=days)
    if kind == "date":
        return moment.replace(hour=0, minute=0, second=0, microsecond=0)
    if kind == "time":
        return moment
    raise ValueError(f"unknown window bound {spec!r}")


def stage_window(recipe: str, stage: str, origin: datetime) -> tuple[datetime, datetime]:
    """[start, end) of an exported stage of a recipe for an event at origin (UTC)."""
    start, end = store.COMPOSITE_RECIPES[recipe]["window_days"][stage]
    return window_bound(origin, start), window_bound(origin, end)


def snapped_grid(xmin: float, ymin: float, xmax: float, ymax: float) -> dict[str, Any]:
    """The 10 m grid of an export covering these UTM bounds (snapped outward, as Earth Engine)."""
    x0, y0 = math.floor(xmin / GRID_M) * GRID_M, math.ceil(ymax / GRID_M) * GRID_M
    x1, y1 = math.ceil(xmax / GRID_M) * GRID_M, math.floor(ymin / GRID_M) * GRID_M
    return {"x0": x0, "y0": y0, "cols": (x1 - x0) // GRID_M, "rows": (y0 - y1) // GRID_M}


def raster_grid(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Grid of a north-up 10 m GeoTIFF: origin, size and EPSG."""
    from osgeo import gdal

    gdal.UseExceptions()
    ds = gdal.Open(str(path))
    gt = ds.GetGeoTransform()
    if gt[1] != GRID_M or gt[5] != -GRID_M or gt[2] or gt[4]:
        raise ValueError(f"{path} is not a north-up {GRID_M} m grid")
    return {
        "x0": gt[0],
        "y0": gt[3],
        "cols": ds.RasterXSize,
        "rows": ds.RasterYSize,
        "epsg": int(ds.GetSpatialRef().GetAuthorityCode(None)),
    }


def utm_epsg(lon: float, lat: float) -> int:
    """UTM zone of a point, as the pipeline chooses the composites' CRS."""
    return (32600 if lat >= 0 else 32700) + math.floor((lon + 180) / 6) + 1


def event_aoi(event_meta: dict[str, Any]) -> tuple[Any, str]:
    """The optical AOI of an event, as aria_coseis.pipeline makes it, and where it came from."""
    url = get_ffm_geojson_url(event_meta["event_id"])
    if url:
        return load_aoi_from_json(url).buffer(AOI_BUFFER_DEG).envelope, url
    box = make_aoi([event_meta["longitude"], event_meta["latitude"]])
    return box, "1 x 1 degree box around the epicentre"


def origin_time(event_meta: dict[str, Any]) -> datetime:
    return datetime.strptime(event_meta["time_utc"][:19], "%Y-%m-%dT%H:%M:%S").replace(
        tzinfo=timezone.utc
    )


def planned_composites(
    root: Path,
    event: str,
    recipe: str,
    origin: datetime,
    reference: str = "med90",
    strict: bool = True,
) -> list[dict[str, Any]]:
    """
    The composites of a recipe for an event: one per path of the reference recipe and stage,
    with its file name, window and what to do ("registered", "reuse" or "export"). A reused
    composite that does not exist yet raises, or with strict=False is planned as "reuse later"
    (with the window its recipe gives).
    """
    definition = store.COMPOSITE_RECIPES[recipe]
    if "window_days" not in definition:
        raise ValueError(f"{recipe} is not made by this module (no window_days)")
    pairs, _ = runs.composite_pairs(root, event, reference)
    if not pairs:
        raise ValueError(f"{event} has no {reference} composite pairs to take the paths from")
    index_path = root / "composites" / event / "composites.json"
    index = json.loads(index_path.read_text())
    reused = definition.get("stage_from_recipe", {})
    plan = []
    for pair in pairs:
        for stage in STAGES:
            item: dict[str, Any] = {"path": pair["path"], "stage": stage}
            if stage in reused:
                source_recipe = reused[stage]
                matches = [
                    (name, e)
                    for name, e in index.items()
                    if not name.startswith("_")
                    and e.get("recipe") == source_recipe
                    and e.get("path") == pair["path"]
                    and e.get("stage") == stage
                ]
                if not matches and not strict:
                    start, end = stage_window(source_recipe, stage, origin)
                    source_name = store.composite_name(
                        "S2", pair["path"], "TOA", "B8", stage, source_recipe, start, end
                    )
                    item.update(action="reuse later", source=source_name)
                elif len(matches) != 1:
                    raise ValueError(
                        f"{event} P{pair['path']}: {recipe} {stage} reuses the {source_recipe} "
                        f"composite, and {len(matches)} exist (make {source_recipe} first)"
                    )
                else:
                    source_name, source_entry = matches[0]
                    start, end = (
                        datetime.strptime(source_entry[k], "%Y-%m-%dT%H:%M:%SZ").replace(
                            tzinfo=timezone.utc
                        )
                        for k in ("start", "end")
                    )
                    item.update(action="reuse", source=source_name)
            else:
                start, end = stage_window(recipe, stage, origin)
                item["action"] = "export"
            item.update(start=start, end=end)
            item["name"] = store.composite_name(
                "S2", pair["path"], "TOA", "B8", stage, recipe, start, end
            )
            if item["name"] in index:
                item["action"] = "registered"
            plan.append(item)
    return plan


def scene_collection(
    ee: Any, roi: Any, path: int, start: datetime, end: datetime, recipe: str
) -> Any:
    """The masked scenes of one relative orbit in [start, end), as the recipe says."""
    definition = store.COMPOSITE_RECIPES[recipe]
    snow = definition.get("snow_mask")
    stamp = "%Y-%m-%dT%H:%M:%S"
    col = (
        ee.ImageCollection(definition["collection"])
        .filterBounds(roi)
        .filterDate(start.strftime(stamp), end.strftime(stamp))
        .filter(ee.Filter.eq("SENSING_ORBIT_NUMBER", path))
        .linkCollection(ee.ImageCollection(CLOUD_SCORE_COLLECTION), ["cs_cdf"])
    )

    def mask(img: Any) -> Any:
        keep = img.select("cs_cdf").gte(definition["cloud_score_min"])
        if snow:
            is_snow = (
                img.normalizedDifference(["B3", "B11"])
                .gte(snow["ndsi_min"])
                .And(img.select("B8").gt(snow["b8_min"]))
                .And(img.select("B3").gt(snow["b3_min"]))
            )
            keep = keep.And(is_snow.Not())
        return img.updateMask(keep)

    return col.map(mask)


def ee_geometries(ee: Any, aoi: Any, grid: dict[str, Any]) -> tuple[Any, Any]:
    """(the AOI rectangle as the pipeline passes it to Earth Engine, the grid's UTM rectangle)."""
    roi = ee.Geometry.Rectangle(list(aoi.bounds))
    x1 = grid["x0"] + grid["cols"] * GRID_M
    y1 = grid["y0"] - grid["rows"] * GRID_M
    region = ee.Geometry.Rectangle(
        [grid["x0"], y1, x1, grid["y0"]], proj=f"EPSG:{grid['epsg']}", geodesic=False
    )
    return roi, region


def check_grid(ee: Any, root: Path, event: str, aoi: Any, reference: str) -> dict[str, Any]:
    """The reference composites' grid; raises unless the AOI reproduces it."""
    pairs, _ = runs.composite_pairs(root, event, reference)
    grids = {json.dumps(raster_grid(p[k]), sort_keys=True) for p in pairs for k in ("ref", "sec")}
    if len(grids) != 1:
        raise ValueError(f"{event}: the {reference} composites are not on one grid")
    grid = json.loads(grids.pop())
    proj = ee.Projection(f"EPSG:{grid['epsg']}")
    roi = ee.Geometry.Rectangle(list(aoi.bounds))
    ring = roi.transform(proj, 0.1).bounds(0.1, proj).coordinates().getInfo()[0]
    xs, ys = [p[0] for p in ring], [p[1] for p in ring]
    predicted = snapped_grid(min(xs), min(ys), max(xs), max(ys))
    if any(predicted[k] != grid[k] for k in predicted):
        raise ValueError(
            f"{event}: the AOI gives the grid {predicted}, not the {reference} grid {grid} "
            "(has the USGS finite-fault model changed?)"
        )
    return grid


def scenes(collection: Any) -> list[str]:
    """Ids of the scenes in a collection (all of them, including those masked out entirely)."""
    return sorted(collection.aggregate_array("system:index").getInfo())


def setup_event(ee: Any, root: Path, event: str, reference: str) -> dict[str, Any]:
    """Event metadata, origin, AOI and the checked grid."""
    events = json.loads((root / "events.json").read_text())
    meta = events.get(event)
    if not meta or not meta.get("event_id"):
        raise ValueError(f"{event} has no USGS record in events.json")
    aoi, aoi_source = event_aoi(meta)
    grid = check_grid(ee, root, event, aoi, reference)
    if grid["epsg"] != utm_epsg(meta["longitude"], meta["latitude"]):
        raise ValueError(f"{event}: the {reference} grid is not in the epicentre's UTM zone")
    return {
        "meta": meta,
        "origin": origin_time(meta),
        "aoi": aoi,
        "aoi_source": aoi_source,
        "grid": grid,
    }


def plan_event(ee: Any, root: Path, event: str, recipe: str, reference: str) -> None:
    """Print what build would do for one event, with the scenes in each window."""
    setup = setup_event(ee, root, event, reference)
    roi, _ = ee_geometries(ee, setup["aoi"], setup["grid"])
    g = setup["grid"]
    print(f"\n{event}  (origin {setup['origin']:%Y-%m-%d %H:%M:%S} UTC)")
    print(f"  AOI: {setup['aoi_source']}")
    print(
        f"  grid: EPSG {g['epsg']}, origin ({g['x0']}, {g['y0']}), {g['cols']} x {g['rows']} "
        f"px = the {reference} grid"
    )
    for item in planned_composites(root, event, recipe, setup["origin"], reference, strict=False):
        what = item["action"]
        if what == "export":
            n = len(
                scenes(scene_collection(ee, roi, item["path"], item["start"], item["end"], recipe))
            )
            what = f"export ({n} scenes in the window)"
        elif what in ("reuse", "reuse later"):
            what = f"{what} {item['source']}"
        print(f"  P{item['path']:<4} {item['stage']:<4} {item['name']}: {what}")


def _read_tasks(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text()) if path.exists() else {}


def _write_tasks(path: Path, tasks: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(tasks, indent=1, sort_keys=True) + "\n")


def _record(item: dict[str, Any], recipe: str, **extra: Any) -> dict[str, Any]:
    return {
        "sensor": "S2",
        "path": item["path"],
        "level": "TOA",
        "band": "B8",
        "stage": item["stage"],
        "recipe": recipe,
        "start": item["start"],
        "end": item["end"],
        **extra,
    }


def register_reused(root: Path, event: str, item: dict[str, Any], recipe: str) -> Path:
    """Register another recipe's composite, unchanged, as this recipe's (a hard link)."""
    source = root / "composites" / event / item["source"]
    return store.add_composite(
        root, event, source, _record(item, recipe, same_file_as=item["source"])
    )


def crop_to_grid(path: Path, got: dict[str, Any], grid: dict[str, Any]) -> bool:
    """
    Crop a raster on the same lattice that contains the grid to exactly the grid (in place).
    Returns False if it cannot (another projection or lattice, or the grid is not inside).
    """
    from osgeo import gdal

    if got["epsg"] != grid["epsg"]:
        return False
    dc, dr = (grid["x0"] - got["x0"]) / GRID_M, (got["y0"] - grid["y0"]) / GRID_M
    if dc != int(dc) or dr != int(dr) or dc < 0 or dr < 0:
        return False
    if dc + grid["cols"] > got["cols"] or dr + grid["rows"] > got["rows"]:
        return False
    tmp = path.with_name(path.name + ".crop.tif")
    gdal.Translate(
        str(tmp),
        str(path),
        srcWin=[int(dc), int(dr), grid["cols"], grid["rows"]],
        creationOptions=["COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES"],
    )
    os.replace(tmp, path)
    return True


def valid_fraction(path: Path) -> float:
    """Fraction of pixels with data (not 0)."""
    import numpy as np
    from osgeo import gdal

    ds = gdal.Open(str(path))
    band = ds.GetRasterBand(1)
    valid = 0
    for row in range(0, ds.RasterYSize, 2048):
        rows = min(2048, ds.RasterYSize - row)
        valid += int(np.count_nonzero(band.ReadAsArray(0, row, ds.RasterXSize, rows)))
    return round(valid / (ds.RasterXSize * ds.RasterYSize), 4)


def start_event(
    ee: Any, root: Path, event: str, recipe: str, reference: str, bucket: str
) -> dict[str, Any]:
    """
    Register the reused stages of one event and start its exports (or find those already
    started). Returns the event's job for finish_exports.
    """
    setup = setup_event(ee, root, event, reference)
    roi, region = ee_geometries(ee, setup["aoi"], setup["grid"])
    g = setup["grid"]
    staging = root / "work" / "composites" / recipe / event
    tasks_path = staging / "tasks.json"
    tasks = _read_tasks(tasks_path)
    store.write_recipes(root)
    job: dict[str, Any] = {
        "event": event,
        "grid": g,
        "staging": staging,
        "tasks": tasks,
        "exports": [],
        "outcome": {},
    }
    for item in planned_composites(root, event, recipe, setup["origin"], reference):
        name = item["name"]
        if item["action"] == "registered":
            job["outcome"][name] = "already registered"
            continue
        if item["action"] == "reuse":
            register_reused(root, event, item, recipe)
            job["outcome"][name] = f"registered (same file as {item['source']})"
            continue
        job["exports"].append(item)
        known = tasks.get(name)
        if known:
            state = ee.data.getTaskStatus(known["task_id"])[0]["state"]
            if state in ("READY", "RUNNING", "COMPLETED"):
                print(f"  {name}: export already {state}")
                continue
            print(f"  {name}: the earlier export is {state}; starting a new one")
        col = scene_collection(ee, roi, item["path"], item["start"], item["end"], recipe)
        image = col.select("B8").median().toUint16().clip(roi)
        prefix = f"COSEIS_Store/{recipe}/{event}/{name[:-4]}"
        task = ee.batch.Export.image.toCloudStorage(
            image=image,
            description=f"{name[:-4]}_{event}"[:100],
            bucket=bucket,
            fileNamePrefix=prefix,
            region=region,
            crs=f"EPSG:{g['epsg']}",
            crsTransform=[GRID_M, 0, g["x0"], 0, -GRID_M, g["y0"]],
            maxPixels=1e13,
            formatOptions={"cloudOptimized": False},
        )
        task.start()
        tasks[name] = {
            "task_id": task.id,
            "prefix": prefix,
            "scenes": scenes(col),
            "aoi_source": setup["aoi_source"],
            "started_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        }
        _write_tasks(tasks_path, tasks)
        print(f"  started export {name} ({len(tasks[name]['scenes'])} scenes)")
    return job


def register_export(
    root: Path, job: dict[str, Any], item: dict[str, Any], recipe: str, reference: str, bucket: str
) -> str:
    """Download a completed export, merge, tag nodata, check the grid and register it."""
    from aria_coseis.optical.gee import assign_nodata, download_from_gcs, merge_and_compress_chips

    name, task, grid = item["name"], job["tasks"][item["name"]], job["grid"]
    files = download_from_gcs(bucket, task["prefix"], str(job["staging"] / "download"))
    if not files:
        return f"export COMPLETED but no files under {task['prefix']}"
    staged = job["staging"] / name
    if len(files) > 1:
        merge_and_compress_chips(files, str(staged))
    else:
        os.replace(files[0], staged)
    assign_nodata(str(staged), nodata_val=0)
    got = raster_grid(staged)
    cropped = got != grid and crop_to_grid(staged, got, grid)
    if raster_grid(staged) != grid:
        return f"NOT registered: grid {got} differs from {grid}"
    record = _record(
        item,
        recipe,
        scenes=task["scenes"],
        aoi=task["aoi_source"],
        grid_of=reference,
        valid_fraction=valid_fraction(staged),
        gee_task=task["task_id"],
        cropped_from=got if cropped else None,
    )
    store.add_composite(
        root, job["event"], staged, record, source=f"gs://{bucket}/{task['prefix']}"
    )
    staged.unlink()  # the store keeps the file (hard link)
    return "registered"


def finish_exports(
    ee: Any,
    root: Path,
    jobs: list[dict[str, Any]],
    recipe: str,
    reference: str,
    bucket: str,
    timeout_h: float,
) -> None:
    """Wait for the jobs' exports, registering each as it completes; records the outcomes."""
    waiting = [(job, item) for job in jobs for item in job["exports"]]
    deadline = time.time() + timeout_h * 3600
    while waiting:
        still = []
        for job, item in waiting:
            status = ee.data.getTaskStatus(job["tasks"][item["name"]]["task_id"])[0]
            if status["state"] in ("READY", "RUNNING"):
                still.append((job, item, status["state"]))
                continue
            if status["state"] == "COMPLETED":
                what = register_export(root, job, item, recipe, reference, bucket)
            else:
                what = f"export {status['state']}: {status.get('error_message', '')}"
            job["outcome"][item["name"]] = what
            print(f"  {job['event']} {item['name']}: {what}")
        waiting = [(job, item) for job, item, _ in still]
        if not waiting:
            break
        if time.time() > deadline:
            for job, item in waiting:
                job["outcome"][item["name"]] = "still running at the timeout (rerun build)"
            break
        running = sum(state == "RUNNING" for _, _, state in still)
        print(f"  {len(waiting)} export(s) not done ({running} running); next check in 2 min")
        time.sleep(120)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="aria-coseis-composites",
        description="Make composites of a recipe for events in the optical data store, on the "
        "grid of their med90 composites.",
    )
    parser.add_argument(
        "command", choices=["plan", "build"], help="plan only reads; build exports"
    )
    choices = sorted(k for k, v in store.COMPOSITE_RECIPES.items() if "window_days" in v)
    parser.add_argument("--recipe", required=True, choices=choices)
    which = parser.add_mutually_exclusive_group(required=True)
    which.add_argument("--event", action="append", help="Event folder name (repeatable)")
    which.add_argument(
        "--testbed", action="store_true", help=f"The {len(TESTBED_EVENTS)} test events"
    )
    parser.add_argument("--reference", default="med90", help="Recipe giving the grid and paths")
    parser.add_argument("--store", help="Store root (default: COSEIS_OPTICAL_DIR)")
    parser.add_argument("--timeout_h", type=float, default=12, help="Stop waiting after (hours)")
    args = parser.parse_args(argv)

    import ee

    root = store.require_store(args.store).resolve()
    ee.Initialize(project=EE_PROJECT)
    events = list(TESTBED_EVENTS) if args.testbed else args.event
    print(f"Recipe {args.recipe}: {store.COMPOSITE_RECIPES[args.recipe]['description']}")
    if args.command == "plan":
        for event in events:
            plan_event(ee, root, event, args.recipe, args.reference)
        return
    bucket = os.environ.get("COSEIS_GCS_BUCKET")
    if not bucket:
        parser.error("COSEIS_GCS_BUCKET is not set")
    jobs, failed = [], {}
    for event in events:
        print(f"\n{event}")
        try:
            jobs.append(start_event(ee, root, event, args.recipe, args.reference, bucket))
        except Exception as e:  # noqa: BLE001 - report the event and go on with the others
            failed[event] = f"{type(e).__name__}: {e}"
            print(f"  FAILED: {failed[event]}")
    print("\nWaiting for the exports")
    finish_exports(ee, root, jobs, args.recipe, args.reference, bucket, args.timeout_h)
    for job in jobs:
        print(f"\n{job['event']}")
        for name, what in job["outcome"].items():
            print(f"  {name}: {what}")
    for event, error in failed.items():
        print(f"\n{event}\n  FAILED: {error}")


if __name__ == "__main__":
    main()
