"""
Composites from the Copernicus Data Space Ecosystem instead of Earth Engine:
`python -m aria_coseis.optical.cdse`.

    python -m aria_coseis.optical.cdse plan --recipe med90snowcdse --event <event>
    python -m aria_coseis.optical.cdse build --recipe med90snowcdse --event <event> [--event ...]

The *cdse recipes (store.COMPOSITE_RECIPES) are the Earth Engine recipes of the same name made from
the Data Space, on the same grid and paths (the event's med90 composites, as composites.py), with
Collection-1 scenes only (processing baseline >= 05.00) and the Level-2A scene classification as the
cloud mask. For each event, path and window:

1. The OData catalogue gives the path's Level-1C products intersecting the AOI rectangle in the
   window. Baselines below the recipe's minimum are left out; per tile and acquisition date the
   highest baseline (then the latest processing) is kept, and matched with the Level-2A product of
   the same acquisition. Scenes without a Level-2A product, or offline, are left out and listed.
2. B03, B08 and B11 of the Level-1C product and the Level-2A SCL (and, for the *cdseprob recipes,
   the Level-2A cloud probability MSK_CLDPRB_20m) are downloaded as single files (not the whole
   product) into work/cdse/<event>/raw/.
3. Per acquisition date, each tile is warped onto the grid (nearest neighbour), harmonized (DN -
   1000), masked (the recipe's SCL classes and cloud probability, dilated, and its snow rule) and the date's tiles
   mosaicked into one layer, work/cdse/<event>/P<path>/dates_<masks>/<date>.tif (0 = no clear
   value).
4. The composite is the per-pixel median of the window's date layers (uint16, nodata 0), clipped to
   the AOI rectangle; the number of dates used per pixel is written to work/cdse/<event>/counts/. It
   is registered in the store with its products.

Date layers serve all recipes whose windows contain the date. build can be rerun: downloads, date
layers and registered composites are kept. Credentials: ~/.netrc, machine
identity.dataspace.copernicus.eu (or dataspace.copernicus.eu). Runs in the autoRIFT environment
(GDAL with JPEG 2000, scipy, pyproj); not Earth Engine.
"""

from __future__ import annotations

import argparse
import json
import netrc
import os
import re
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

from aria_coseis.optical import composites, runs, store

CATALOGUE = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
DOWNLOAD = "https://download.dataspace.copernicus.eu/odata/v1/Products"
TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
)
NETRC_MACHINES = ("identity.dataspace.copernicus.eu", "dataspace.copernicus.eu")
DOWNLOAD_THREADS = 4  # the Data Space allows four concurrent downloads per user
L1C_BANDS = ("B03", "B08", "B11")
HARMONIZATION_OFFSET = 1000  # RADIO_ADD_OFFSET of processing baselines >= 04.00
NAME = re.compile(
    r"^(S2[ABCD])_MSIL(1C|2A)_(\d{8}T\d{6})_N(\d{4})_R(\d{3})_T(\w{5})_(\d{8}T\d{6})(\.SAFE)?$"
)
STAMP = "%Y%m%dT%H%M%S"


def parse_name(name: str) -> dict[str, Any] | None:
    """Fields of a Sentinel-2 product name, or None if it is not one."""
    m = NAME.match(name)
    if not m:
        return None
    sat, level, sensing, baseline, orbit, tile, processed, _ = m.groups()
    return {
        "name": name.removesuffix(".SAFE"),
        "satellite": sat,
        "level": level,
        "sensing": datetime.strptime(sensing, STAMP).replace(tzinfo=timezone.utc),
        "baseline": f"{baseline[:2]}.{baseline[2:]}",
        "orbit": int(orbit),
        "tile": tile,
        "processed": processed,
    }


def choose_scenes(products: list[dict[str, Any]], baseline_min: str) -> list[dict[str, Any]]:
    """
    One product per tile and acquisition date: baselines >= baseline_min ("NN.NN"), then the
    highest baseline and the latest processing. Sorted by sensing time and tile.
    """
    best: dict[tuple[str, str], dict[str, Any]] = {}
    for p in products:
        if p["baseline"] < baseline_min:
            continue
        key = (p["tile"], p["sensing"].strftime("%Y%m%d"))
        if key not in best or (p["baseline"], p["processed"]) > (
            best[key]["baseline"],
            best[key]["processed"],
        ):
            best[key] = p
    return sorted(best.values(), key=lambda p: (p["sensing"], p["tile"]))


def harmonize(dn: Any) -> Any:
    """Harmonized values (the scale of COPERNICUS/S2_HARMONIZED): DN - 1000, at least 1 where DN > 0."""
    import numpy as np

    out = np.clip(dn.astype(np.int32) - HARMONIZATION_OFFSET, 1, 65535)
    return np.where(dn > 0, out, 0).astype(np.uint16)


def snow(b3: Any, b8: Any, b11: Any, rule: dict[str, Any]) -> Any:
    """The recipe's snow rule on harmonized values (NDSI of B3 and B11, B8 and B3 thresholds)."""
    import numpy as np

    b3f, b11f = b3.astype(np.float32), b11.astype(np.float32)
    with np.errstate(divide="ignore", invalid="ignore"):
        ndsi = (b3f - b11f) / (b3f + b11f)
    return (ndsi >= rule["ndsi_min"]) & (b8 > rule["b8_min"]) & (b3 > rule["b3_min"])


def cloud_mask(
    scl: Any,
    classes: tuple[int, ...] | list[int],
    dilation_px: int,
    probability: Any = None,
    probability_min: int | None = None,
) -> Any:
    """
    Pixels of the masked scene classes, and with probability_min those whose cloud probability (%)
    is at least that, dilated by dilation_px (on the SCL's own 20 m grid).
    """
    import numpy as np
    from scipy import ndimage

    masked = np.isin(scl, classes)
    if probability_min is not None:
        masked |= probability >= probability_min
    if dilation_px:
        masked = ndimage.binary_dilation(masked, structure=np.ones((3, 3)), iterations=dilation_px)
    return masked


def nan_median(stack: Any) -> tuple[Any, Any]:
    """Per-pixel median of the non-zero values of a (dates, rows, cols) uint16 stack, and the count."""
    import warnings

    import numpy as np

    values = stack.astype(np.float32)
    values[stack == 0] = np.nan
    count = np.isfinite(values).sum(axis=0).astype(np.uint16)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN pixels: no clear date
        med = np.nanmedian(values, axis=0)
    out = np.where(count > 0, np.floor(med), 0).astype(np.uint16)  # truncated, as toUint16()
    return out, count


class Session:
    """Authenticated requests to the Data Space, with the token renewed before it expires."""

    def __init__(self) -> None:
        auth = None
        for machine in NETRC_MACHINES:
            auth = netrc.netrc().authenticators(machine)
            if auth:
                break
        if not auth:
            raise RuntimeError(
                f"no Data Space credentials in ~/.netrc ({', '.join(NETRC_MACHINES)})"
            )
        self._user, self._password = auth[0], auth[2]
        self._token, self._expires = "", 0.0
        self._lock = threading.Lock()

    def headers(self) -> dict[str, str]:
        with self._lock:
            if time.time() > self._expires - 60:
                r = requests.post(
                    TOKEN_URL,
                    data={
                        "grant_type": "password",
                        "client_id": "cdse-public",
                        "username": self._user,
                        "password": self._password,
                    },
                    timeout=60,
                )
                r.raise_for_status()
                token = r.json()
                self._token = token["access_token"]
                self._expires = time.time() + float(token.get("expires_in", 600))
            return {"Authorization": f"Bearer {self._token}"}

    def get(self, url: str, stream: bool = False, tries: int = 6) -> requests.Response:
        for attempt in range(tries):
            try:
                r = requests.get(url, headers=self.headers(), stream=stream, timeout=300)
            except requests.RequestException:
                if attempt == tries - 1:
                    raise
                time.sleep(10 * (attempt + 1))
                continue
            if r.status_code == 401:
                self._expires = 0.0
            if r.status_code in (401, 429) or r.status_code >= 500:
                if attempt == tries - 1:
                    r.raise_for_status()
                time.sleep(10 * (attempt + 1))
                continue
            r.raise_for_status()
            return r
        raise RuntimeError(f"no response from {url}")


def search(
    product_type: str, path: int, bounds: tuple[float, ...], start: datetime, end: datetime
) -> list[dict[str, Any]]:
    """Products of a type (S2MSI1C, S2MSI2A) of a relative orbit over the bounds in [start, end)."""
    lon0, lat0, lon1, lat1 = bounds
    poly = f"POLYGON(({lon0} {lat0},{lon1} {lat0},{lon1} {lat1},{lon0} {lat1},{lon0} {lat0}))"
    stamp = "%Y-%m-%dT%H:%M:%S.000Z"
    flt = (
        "Collection/Name eq 'SENTINEL-2' and "
        "Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' and "
        f"att/OData.CSC.StringAttribute/Value eq '{product_type}') and "
        "Attributes/OData.CSC.IntegerAttribute/any(att:att/Name eq 'relativeOrbitNumber' and "
        f"att/OData.CSC.IntegerAttribute/Value eq {path}) and "
        f"OData.CSC.Intersects(area=geography'SRID=4326;{poly}') and "
        f"ContentDate/Start ge {start.strftime(stamp)} and ContentDate/Start lt {end.strftime(stamp)}"
    )
    out, skip = [], 0
    while True:
        params = {"$filter": flt, "$top": 1000, "$skip": skip, "$select": "Id,Name,Online"}
        for attempt in range(6):
            r = requests.get(CATALOGUE, params=params, timeout=120)
            if r.status_code == 200:
                break
            time.sleep(10 * (attempt + 1))
        r.raise_for_status()
        page = r.json().get("value", [])
        for p in page:
            fields = parse_name(p["Name"])
            if fields:
                out.append({**fields, "id": p["Id"], "online": p.get("Online", True)})
        if len(page) < 1000:
            return out
        skip += 1000


def find_scenes(
    path: int, bounds: tuple[float, ...], start: datetime, end: datetime, baseline_min: str
) -> dict[str, Any]:
    """The window's chosen Level-1C scenes, each with its Level-2A match, and those left out."""
    l1c = choose_scenes(search("S2MSI1C", path, bounds, start, end), baseline_min)
    l2a = {
        (p["tile"], p["sensing"].strftime(STAMP)): p
        for p in choose_scenes(search("S2MSI2A", path, bounds, start, end), baseline_min)
    }
    scenes, left_out = [], {}
    for p in l1c:
        match = l2a.get((p["tile"], p["sensing"].strftime(STAMP)))
        if not match:
            left_out[p["name"]] = "no Level-2A product of the acquisition"
        elif not (p["online"] and match["online"]):
            left_out[p["name"]] = "offline"
        else:
            scenes.append({"l1c": p, "l2a": match})
    return {"scenes": scenes, "left_out": left_out}


def recipe_mask(recipe: str) -> dict[str, Any]:
    """The recipe's cloud mask: its cdse_mask, or the scene classes of the *cdse recipes."""
    return store.COMPOSITE_RECIPES[recipe].get("cdse_mask") or {
        "scl": store.CDSE_SCL_MASKED,
        "dilation_px": store.CDSE_SCL_DILATION_PX,
    }


def mask_signature(recipe: str) -> str:
    """Folder name of the date layers: they depend on the masks only, so recipes can share them."""
    d = store.COMPOSITE_RECIPES[recipe]
    snow_rule = d.get("snow_mask")
    mask = recipe_mask(recipe)
    tag = f"scl{''.join(map(str, mask['scl']))}d{mask['dilation_px']}"
    if mask.get("cloud_probability_min") is not None:
        tag += f"_prb{mask['cloud_probability_min']}"
    if snow_rule:
        tag += f"_snow{snow_rule['ndsi_min']}-{snow_rule['b8_min']}-{snow_rule['b3_min']}"
    return f"dates_{tag}"


def scene_dir(work: Path, scene: dict[str, Any]) -> Path:
    p = scene["l1c"]
    return work / "raw" / f"{p['tile']}_{p['sensing'].strftime(STAMP)}"


def _nodes(session: Session, product: dict[str, Any], *parts: str) -> list[str]:
    url = f"{DOWNLOAD}({product['id']})/Nodes({product['name']}.SAFE)"
    for part in parts:
        url += f"/Nodes({part})"
    return [n["Name"] for n in session.get(url + "/Nodes").json()["result"]]


# A connection broken while a file streams (the Data Space resets long downloads now and then):
# Session.get retries the request, not the body, so the whole file is fetched again.
STREAM_ERRORS = (
    requests.exceptions.ChunkedEncodingError,
    requests.exceptions.ConnectionError,
    requests.exceptions.Timeout,
)


def _fetch(
    session: Session, product: dict[str, Any], parts: list[str], dst: Path, tries: int = 6
) -> None:
    url = f"{DOWNLOAD}({product['id']})/Nodes({product['name']}.SAFE)"
    for part in parts:
        url += f"/Nodes({part})"
    part_file = dst.with_suffix(dst.suffix + ".part")
    for attempt in range(tries):
        try:
            with session.get(url + "/$value", stream=True) as r, open(part_file, "wb") as fh:
                for chunk in r.iter_content(1 << 22):
                    fh.write(chunk)
            break
        except STREAM_ERRORS:
            if attempt == tries - 1:
                raise
            time.sleep(10 * (attempt + 1))
    os.replace(part_file, dst)


def download_scene(
    session: Session, work: Path, scene: dict[str, Any], probability: bool = False
) -> str:
    """
    B03, B08, B11 of the Level-1C product and the SCL of the Level-2A one, and with probability the
    Level-2A cloud probability (MSK_CLDPRB_20m); skips files present.
    """
    folder = scene_dir(work, scene)
    done = folder / "done.json"
    status = "present"
    if not done.exists():
        _download_bands(session, folder, scene)
        status = "downloaded"
    if probability and not (folder / "CLDPRB.jp2").exists():
        l2a = scene["l2a"]
        granule = _nodes(session, l2a, "GRANULE")[0]
        files = _nodes(session, l2a, "GRANULE", granule, "QI_DATA")
        name = next(f for f in files if f.startswith("MSK_CLDPRB_20m"))
        _fetch(session, l2a, ["GRANULE", granule, "QI_DATA", name], folder / "CLDPRB.jp2")
        status = "downloaded"
    return status


def _download_bands(session: Session, folder: Path, scene: dict[str, Any]) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    done = folder / "done.json"
    l1c, l2a = scene["l1c"], scene["l2a"]
    granule = _nodes(session, l1c, "GRANULE")[0]
    files = _nodes(session, l1c, "GRANULE", granule, "IMG_DATA")
    for band in L1C_BANDS:
        name = next(f for f in files if f.endswith(f"_{band}.jp2"))
        if not (folder / f"{band}.jp2").exists():
            _fetch(session, l1c, ["GRANULE", granule, "IMG_DATA", name], folder / f"{band}.jp2")
    granule2 = _nodes(session, l2a, "GRANULE")[0]
    files2 = _nodes(session, l2a, "GRANULE", granule2, "IMG_DATA", "R20m")
    name = next(f for f in files2 if f.endswith("_SCL_20m.jp2"))
    if not (folder / "SCL.jp2").exists():
        _fetch(session, l2a, ["GRANULE", granule2, "IMG_DATA", "R20m", name], folder / "SCL.jp2")
    done.write_text(json.dumps({"l1c": l1c["name"], "l2a": l2a["name"]}) + "\n")


def _warp(src: Any, grid: dict[str, Any], bounds: tuple[float, ...]) -> Any:
    """A raster (path or dataset) warped by nearest neighbour onto the grid's window at bounds."""
    from osgeo import gdal

    ds = gdal.Warp(
        "",
        src,
        format="MEM",
        dstSRS=f"EPSG:{grid['epsg']}",
        outputBounds=bounds,
        xRes=composites.GRID_M,
        yRes=composites.GRID_M,
        resampleAlg="near",
    )
    return ds.GetRasterBand(1).ReadAsArray()


def _tile_window(path: Path, grid: dict[str, Any]) -> tuple[int, int, int, int] | None:
    """The grid's (col0, row0, cols, rows) covering a tile (in any UTM zone), or None outside."""
    import math

    from osgeo import gdal
    from pyproj import Transformer

    ds = gdal.Open(str(path))
    gt = ds.GetGeoTransform()
    src = int(ds.GetSpatialRef().GetAuthorityCode(None))
    xs = [gt[0], gt[0] + ds.RasterXSize * gt[1]]
    ys = [gt[3], gt[3] + ds.RasterYSize * gt[5]]
    if src != grid["epsg"]:
        t = Transformer.from_crs(src, grid["epsg"], always_xy=True)
        corners = [t.transform(x, y) for x in xs for y in ys]
        xs, ys = [c[0] for c in corners], [c[1] for c in corners]
    g = composites.GRID_M
    c0 = max(0, math.floor((min(xs) - grid["x0"]) / g))
    c1 = min(grid["cols"], math.ceil((max(xs) - grid["x0"]) / g))
    r0 = max(0, math.floor((grid["y0"] - max(ys)) / g))
    r1 = min(grid["rows"], math.ceil((grid["y0"] - min(ys)) / g))
    if c1 <= c0 or r1 <= r0:
        return None
    return c0, r0, c1 - c0, r1 - r0


def date_layer(task: dict[str, Any]) -> str:
    """Make one date's masked layer on the grid from its tiles (in a worker process)."""
    import numpy as np
    from osgeo import gdal

    gdal.UseExceptions()
    out, grid, recipe = Path(task["out"]), task["grid"], task["recipe"]
    if out.exists():
        return "present"
    definition = store.COMPOSITE_RECIPES[recipe]
    rule = definition.get("snow_mask")
    mask = recipe_mask(recipe)
    probability_min = mask.get("cloud_probability_min")
    layer = np.zeros((grid["rows"], grid["cols"]), dtype=np.uint16)
    g = composites.GRID_M
    for folder in map(Path, task["scene_dirs"]):
        window = _tile_window(folder / "B08.jp2", grid)
        if window is None:
            continue
        c0, r0, nc, nr = window
        bounds = (
            grid["x0"] + c0 * g,
            grid["y0"] - (r0 + nr) * g,
            grid["x0"] + (c0 + nc) * g,
            grid["y0"] - r0 * g,
        )
        b8 = harmonize(_warp(str(folder / "B08.jp2"), grid, bounds))
        keep = b8 > 0
        scl_ds = gdal.Open(str(folder / "SCL.jp2"))
        probability = None
        if probability_min is not None:  # same 20 m grid as the SCL
            probability = gdal.Open(str(folder / "CLDPRB.jp2")).ReadAsArray()
        clouds = cloud_mask(
            scl_ds.GetRasterBand(1).ReadAsArray(),
            mask["scl"],
            mask["dilation_px"],
            probability,
            probability_min,
        )
        mem = gdal.GetDriverByName("MEM").Create(
            "", scl_ds.RasterXSize, scl_ds.RasterYSize, 1, gdal.GDT_Byte
        )
        mem.SetGeoTransform(scl_ds.GetGeoTransform())
        mem.SetProjection(scl_ds.GetProjection())
        mem.GetRasterBand(1).WriteArray(clouds.astype(np.uint8) + 1)  # 1 clear, 2 masked
        mem.GetRasterBand(1).SetNoDataValue(0)
        cloud10 = _warp(mem, grid, bounds)
        keep &= cloud10 == 1
        if rule:
            b3 = harmonize(_warp(str(folder / "B03.jp2"), grid, bounds))
            b11 = harmonize(_warp(str(folder / "B11.jp2"), grid, bounds))
            keep &= ~snow(b3, b8, b11, rule)
        sub = layer[r0 : r0 + nr, c0 : c0 + nc]
        fill = keep & (sub == 0)  # tiles of one acquisition: the first clear value
        sub[fill] = b8[fill]
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".part.tif")
    ds = gdal.GetDriverByName("GTiff").Create(
        str(tmp),
        grid["cols"],
        grid["rows"],
        1,
        gdal.GDT_UInt16,
        options=["COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_SAFER"],
    )
    ds.SetGeoTransform((grid["x0"], g, 0, grid["y0"], 0, -g))
    from osgeo import osr

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(grid["epsg"])
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(0)
    band.WriteArray(layer)
    ds = None
    os.replace(tmp, out)
    return "made"


def aoi_clip(grid: dict[str, Any], bounds: tuple[float, ...], row0: int, rows: int) -> Any:
    """Pixels of grid rows [row0, row0 + rows) whose centres lie in the lon/lat AOI rectangle."""
    import numpy as np
    from pyproj import Transformer

    g = composites.GRID_M
    x = grid["x0"] + (np.arange(grid["cols"]) + 0.5) * g
    y = grid["y0"] - (np.arange(row0, row0 + rows) + 0.5) * g
    X, Y = np.meshgrid(x, y)
    lon, lat = Transformer.from_crs(grid["epsg"], 4326, always_xy=True).transform(X, Y)
    lon0, lat0, lon1, lat1 = bounds
    return (lon >= lon0) & (lon <= lon1) & (lat >= lat0) & (lat <= lat1)


def median_composite(
    layers: list[Path], grid: dict[str, Any], bounds: tuple[float, ...], out: Path, count_out: Path
) -> None:
    """The median of the date layers (uint16, nodata 0, clipped to the AOI) and the date count."""
    import numpy as np
    from osgeo import gdal, osr

    gdal.UseExceptions()
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(grid["epsg"])
    g = composites.GRID_M
    outs = {}
    for path, dtype in ((out, gdal.GDT_UInt16), (count_out, gdal.GDT_UInt16)):
        path.parent.mkdir(parents=True, exist_ok=True)
        ds = gdal.GetDriverByName("GTiff").Create(
            str(path.with_suffix(".part.tif")),
            grid["cols"],
            grid["rows"],
            1,
            dtype,
            options=["COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES", "BIGTIFF=IF_SAFER"],
        )
        ds.SetGeoTransform((grid["x0"], g, 0, grid["y0"], 0, -g))
        ds.SetProjection(srs.ExportToWkt())
        outs[path] = ds
    outs[out].GetRasterBand(1).SetNoDataValue(0)
    sources = [gdal.Open(str(p)) for p in layers]
    step = 256
    for row in range(0, grid["rows"], step):
        rows = min(step, grid["rows"] - row)
        stack = np.stack(
            [s.GetRasterBand(1).ReadAsArray(0, row, grid["cols"], rows) for s in sources]
        )
        med, count = nan_median(stack)
        inside = aoi_clip(grid, bounds, row, rows)
        med[~inside] = 0
        count[~inside] = 0
        outs[out].GetRasterBand(1).WriteArray(med, 0, row)
        outs[count_out].GetRasterBand(1).WriteArray(count, 0, row)
    sources.clear()
    for path in list(outs):
        outs[path].FlushCache()
        outs[path] = None  # closes the file
        os.replace(path.with_suffix(".part.tif"), path)


def event_setup(root: Path, event: str, reference: str) -> dict[str, Any]:
    """Event metadata, origin, AOI (as the pipeline makes it) and the reference composites' grid."""
    from aria_coseis.aoi import load_aoi_from_json, make_aoi
    from aria_coseis.optical.farfield import finite_fault_url

    meta = json.loads((root / "events.json").read_text()).get(event)
    if not meta or not meta.get("event_id"):
        raise ValueError(f"{event} has no USGS record in events.json")
    url = finite_fault_url(meta["event_id"])
    if url:
        aoi = load_aoi_from_json(url).buffer(composites.AOI_BUFFER_DEG).envelope
    else:
        aoi = make_aoi([meta["longitude"], meta["latitude"]])
    pairs, _ = runs.composite_pairs(root, event, reference)
    grids = {
        json.dumps(composites.raster_grid(p[k]), sort_keys=True)
        for p in pairs
        for k in ("ref", "sec")
    }
    if len(grids) != 1:
        raise ValueError(f"{event}: the {reference} composites are not on one grid")
    grid = json.loads(grids.pop())
    if grid["epsg"] != composites.utm_epsg(meta["longitude"], meta["latitude"]):
        raise ValueError(f"{event}: the {reference} grid is not in the epicentre's UTM zone")
    return {"meta": meta, "origin": composites.origin_time(meta), "aoi": aoi, "grid": grid}


def plan_event(
    root: Path, event: str, recipe: str, reference: str, setup: dict[str, Any] | None = None
) -> list[dict[str, Any]]:
    """The recipe's composites for an event, with the scenes of each window to export."""
    setup = setup or event_setup(root, event, reference)
    baseline_min = store.COMPOSITE_RECIPES[recipe]["processing_baseline_min"]
    items = composites.planned_composites(root, event, recipe, setup["origin"], reference)
    for item in items:
        if item["action"] == "export":
            item.update(
                find_scenes(
                    item["path"], setup["aoi"].bounds, item["start"], item["end"], baseline_min
                )
            )
    return items


def build_event(
    root: Path, event: str, recipe: str, reference: str, code: dict[str, Any], workers: int
) -> dict[str, str]:
    """Download, make the date layers and register the recipe's composites for one event."""
    setup = event_setup(root, event, reference)
    work = root / "work" / "cdse" / event
    store.write_recipes(root)
    outcome: dict[str, str] = {}
    items = plan_event(root, event, recipe, reference, setup)
    session = Session()
    for item in items:
        name = item["name"]
        if item["action"] == "registered":
            outcome[name] = "already registered"
            continue
        if item["action"] == "reuse":
            composites.register_reused(root, event, item, recipe, code)
            outcome[name] = f"registered (same file as {item['source']})"
            continue
        scenes = item["scenes"]
        print(
            f"  P{item['path']} {item['stage']}: {len(scenes)} scenes ({len(item['left_out'])} left out)"
        )
        with ThreadPoolExecutor(DOWNLOAD_THREADS) as pool:
            probability = recipe_mask(recipe).get("cloud_probability_min") is not None
            for status in pool.map(
                lambda s: download_scene(session, work, s, probability), scenes
            ):
                pass
        layer_dir = mask_signature(recipe)
        by_date: dict[str, list[str]] = defaultdict(list)
        for s in scenes:
            by_date[s["l1c"]["sensing"].strftime("%Y%m%d")].append(str(scene_dir(work, s)))
        tasks = [
            {
                "out": str(work / f"P{item['path']}" / layer_dir / f"{date}.tif"),
                "grid": setup["grid"],
                "recipe": recipe,
                "scene_dirs": dirs,
            }
            for date, dirs in sorted(by_date.items())
        ]
        with ProcessPoolExecutor(workers) as pool:
            list(pool.map(date_layer, tasks))
        staged = work / "staged" / name
        count = work / "counts" / name.replace(".tif", "_count.tif")
        median_composite(
            [Path(t["out"]) for t in tasks], setup["grid"], setup["aoi"].bounds, staged, count
        )
        record = composites._record(
            item,
            recipe,
            scenes=[s["l1c"]["name"] for s in scenes],
            cloud_mask_products=[s["l2a"]["name"] for s in scenes],
            left_out=item["left_out"],
            dates=len(tasks),
            valid_fraction=composites.valid_fraction(staged),
            registered_by=code,
        )
        store.add_composite(root, event, staged, record, source="Copernicus Data Space (OData)")
        outcome[name] = f"registered ({len(scenes)} scenes, {len(tasks)} dates)"
        print(f"  {name}: {outcome[name]}", flush=True)
    return outcome


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="aria-coseis-cdse",
        description="Make composites of a *cdse recipe from the Copernicus Data Space, on the grid "
        "of the events' med90 composites.",
    )
    parser.add_argument(
        "command", choices=["plan", "build"], help="plan only reads; build downloads"
    )
    choices = sorted(
        k
        for k, v in store.COMPOSITE_RECIPES.items()
        if v.get("builder") == "aria_coseis.optical.cdse"
    )
    parser.add_argument("--recipe", required=True, choices=choices)
    parser.add_argument(
        "--event", action="append", required=True, help="Event folder (repeatable)"
    )
    parser.add_argument("--reference", default="med90", help="Recipe giving the grid and paths")
    parser.add_argument("--store", help="Store root (default: COSEIS_OPTICAL_DIR)")
    parser.add_argument("--workers", type=int, default=4, help="Processes for the date layers")
    args = parser.parse_args(argv)
    root = store.require_store(args.store).resolve()
    print(f"Recipe {args.recipe}: {store.COMPOSITE_RECIPES[args.recipe]['description']}")
    if args.command == "plan":
        for event in args.event:
            print(f"\n{event}")
            for item in plan_event(root, event, args.recipe, args.reference):
                what = item["action"]
                if what == "export":
                    base = sorted({s["l1c"]["baseline"] for s in item["scenes"]})
                    what = (
                        f"export: {len(item['scenes'])} scenes on "
                        f"{len({s['l1c']['sensing'].date() for s in item['scenes']})} dates, baselines "
                        f"{base}, {len(item['left_out'])} left out"
                    )
                print(f"  P{item['path']:<4} {item['stage']:<4} {item['name']}: {what}")
        return
    # Composites are permanent and not versioned by code, so they are made from committed code
    code = runs.code_signature(store.code_provenance())["aria_coseis"]
    if code.get("dirty") or not code.get("commit"):
        parser.error(f"build needs committed aria_coseis code (got {code})")
    failed = []
    for event in args.event:
        print(f"\n{event}")
        try:
            for name, what in build_event(
                root, event, args.recipe, args.reference, code, args.workers
            ).items():
                print(f"  {name}: {what}")
        except Exception as e:  # noqa: BLE001 - report the event and go on with the others
            print(f"  FAILED: {type(e).__name__}: {e}")
            failed.append(event)
    if failed:
        sys.exit(f"{len(failed)} of {len(args.event)} events failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
