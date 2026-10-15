"""
Versioned store for optical composites and autoRIFT products: `python -m aria_coseis.optical.store`.

Layout under config.OPTICAL_DIR:

    README.md
    events.json                          USGS metadata for every event folder name
    composites/recipes.json              how each composite recipe is made
    composites/<event>/composites.json   one record per composite (recipe, window, checksum, source)
    composites/<event>/<name>.tif        GEE exports; read-only and never modified
    products/<version>/VERSION.json      the code, parameters and recipes that define the version
    products/<version>/catalog.csv       one row per event: status, paths used and skipped (and why)
    products/<version>/<event>/          the products of that version
    work/<version>/<event>/              intermediates; safe to delete

Rules: a composite is written once, and a different recipe or time window is a different file. A
product version is made by one code state with one set of parameters; changing either makes a new
version. Versions never share folders.

No GDAL here, so the store can be managed from the SAR environment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import sys
from collections.abc import Callable
from datetime import datetime, timedelta, timezone
from importlib.util import find_spec
from pathlib import Path
from typing import Any

import requests

from aria_coseis import config

USGS_EVENT_QUERY = "https://earthquake.usgs.gov/fdsnws/event/1/query"

# How each composite recipe is made. Ids are permanent: once composites exist for a recipe its
# definition must not change (add a new recipe instead).
COMPOSITE_RECIPES: dict[str, dict[str, Any]] = {
    "med90": {
        "description": (
            "Per-pixel median of the cloud-free scenes of one relative orbit, 90 days before and "
            "90 days after the event (the composites of the first optical products)"
        ),
        "collection": "COPERNICUS/S2_HARMONIZED",
        "level": "TOA",
        "band": "B8",
        "cloud_mask": "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED cs_cdf >= 0.65",
        "snow_mask": None,
        "reducer": "median",
        "grouping": "SENSING_ORBIT_NUMBER",
        "pre_window": "[event - 90 d, event)",
        "post_window": "[event, event + 90 d)",
        "crs": "UTM zone of the epicenter",
        "scale_m": 10,
        "dtype": "uint16",
    },
}

# Recipes made by aria_coseis.optical.composites, on the grid of the event's med90 composites and
# for the same paths. "window_days" gives each exported stage's [start, end) from the origin time
# t: ["time", d] is t + d days, ["date", d] is 00:00 UTC of the day t + d days falls on (so med90 is
# PRE [["date", -90], ["time", 0]], POST [["time", 0], ["date", 90]]). "stage_from_recipe" names
# the recipe whose composite of that stage is used unchanged (registered again, hard link), or
# [recipe, stage] to use another stage of it. A "control" recipe has no earthquake between its
# PRE and POST: its offsets measure the error of that kind of pair.
_STORE_RECIPE = {
    "collection": "COPERNICUS/S2_HARMONIZED",
    "level": "TOA",
    "band": "B8",
    "cloud_mask": "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED cs_cdf >= 0.65",
    "cloud_score_min": 0.65,
    "reducer": "median",
    "grouping": "SENSING_ORBIT_NUMBER",
    "crs": "UTM zone of the epicenter",
    "grid": (
        "the grid of the event's med90 composites: the AOI rectangle (USGS finite-fault model "
        "buffered by 0.15 degrees, or a 1 x 1 degree box) in UTM, snapped outward to 10 m"
    ),
    "clip": "to the AOI rectangle (longitude/latitude), as med90",
    "scale_m": 10,
    "dtype": "uint16",
}
_SNOW_MASK = {
    "rule": (
        "per scene, before the median: snow where NDSI = (B3 - B11) / (B3 + B11) >= 0.4, B8 > 1100 "
        "and B3 > 1000 (TOA reflectance x 10^4; the thresholds of the MODIS snow algorithm, "
        "SNOMAP); snow pixels are left out"
    ),
    "ndsi_min": 0.4,
    "b8_min": 1100,
    "b3_min": 1000,
}
COMPOSITE_RECIPES.update(
    {
        "seas90": {
            **_STORE_RECIPE,
            "description": (
                "Season-matched PRE: the med90 POST window one year (365 d) earlier, so that PRE "
                "and POST share the season and sun elevation. POST is the med90 POST composite"
            ),
            "snow_mask": None,
            "pre_window": "[event - 365 d, 00:00 UTC of event + 90 d - 365 d)",
            "post_window": "the med90 POST composite",
            "window_days": {"PRE": [["time", -365], ["date", -275]]},
            "stage_from_recipe": {"POST": "med90"},
        },
        "med90snow": {
            **_STORE_RECIPE,
            "description": "med90 windows, with snow masked in every scene before the median",
            "snow_mask": _SNOW_MASK,
            "pre_window": "[00:00 UTC of event - 90 d, event)",
            "post_window": "[event, 00:00 UTC of event + 90 d)",
            "window_days": {
                "PRE": [["date", -90], ["time", 0]],
                "POST": [["time", 0], ["date", 90]],
            },
        },
        "seas90snow": {
            **_STORE_RECIPE,
            "description": (
                "seas90 with snow masked in every scene before the median. POST is the med90snow "
                "POST composite"
            ),
            "snow_mask": _SNOW_MASK,
            "pre_window": "[event - 365 d, 00:00 UTC of event + 90 d - 365 d)",
            "post_window": "the med90snow POST composite",
            "window_days": {"PRE": [["time", -365], ["date", -275]]},
            "stage_from_recipe": {"POST": "med90snow"},
        },
        "med30": {
            **_STORE_RECIPE,
            "description": (
                "Same-year short windows: the cloud-free scenes in the 30 days before and the 30 "
                "days after the event. Less change of season and sun than med90, more cloud gaps"
            ),
            "snow_mask": None,
            "pre_window": "[00:00 UTC of event - 30 d, event)",
            "post_window": "[event, 00:00 UTC of event + 30 d)",
            "window_days": {
                "PRE": [["date", -30], ["time", 0]],
                "POST": [["time", 0], ["date", 30]],
            },
        },
        "seas90ctrl": {
            **_STORE_RECIPE,
            "description": (
                "No-earthquake control for seas90: PRE is the seas90 PRE window one more year "
                "earlier (two years before the event), POST is the seas90 PRE composite (one year "
                "before). Its offsets measure the error of a one-year, same-season pair"
            ),
            "snow_mask": None,
            "pre_window": "[event - 730 d, 00:00 UTC of event + 90 d - 730 d)",
            "post_window": "the seas90 PRE composite",
            "window_days": {"PRE": [["time", -730], ["date", -640]]},
            "stage_from_recipe": {"POST": ["seas90", "PRE"]},
            "control": True,
        },
    }
)

# The first products ("v0"), made 2026-10-01 to 2026-10-13 in data/GEE_Optical_Downloads.
V0_OLD_SCRIPT_EVENTS = (
    "m_75_elbistan_earthquake_kahramanmaras_earthquake_sequence",
    "m_78_pazarcik_earthquake_kahramanmaras_earthquake_sequence",
)
V0_UTM_FIX_EVENTS = (
    "m_60_73_km_ssw_of_kaktovik_alaska",
    "m_62_47_km_nw_of_mosquito_lake_alaska",
    "m_63_48_km_nw_of_mosquito_lake_alaska",
    "m_64_89_km_sw_of_kaktovik_alaska",
    "m_64_chukotskiy_avtonomnyy_okrug_russia",
    "m_70_2025_hubbard_glacier_earthquake",
    "m_77_2025_mandalay_burma_myanmar_earthquake",
)
V0_ARIA_CODE = {
    "develop-8b7b7f7": (
        "scripts/batch_autorift.py at develop 8b7b7f7, run 2026-10-01 (processing identical to "
        "aria_coseis 405ae83)"
    ),
    "aria-405ae83": "aria_coseis.optical.autorift at 405ae83 (main), run 2026-10-07 to 2026-10-13",
    "aria-971b216": (
        "aria_coseis.optical.autorift with the projection fix (971b216; issue #41, PR #42), run 2026-10-13 "
        "from a frozen copy of the identical commit bb57e9a"
    ),
}
V0_VERSION = {
    "version": "v0",
    "description": (
        "First 90 m optical products, registered from data/GEE_Optical_Downloads (frozen). Made "
        "before the store existed, by three code states (catalog column aria_code)."
    ),
    "frozen_source": "data/GEE_Optical_Downloads (production clone)",
    "composite_recipe": "med90",
    "code": {
        "aria_coseis": V0_ARIA_CODE,
        "hyp3_autorift": (
            "cmspeed/hyp3-autorift tag coseis-v0 (8034761): develop ab4d9b5 plus the working-tree "
            "edits of 2026-07-22 that every v0 product used"
        ),
        "autoRIFT": "conda-forge autorift 2.0.0 (unmodified)",
    },
    "parameters": {
        "grid_m": 90,
        "chip_sizes_m": [360, 720],
        "search_limit_px": 4,
        "effective_search_px": 6,
        "oversample_ratio": 128,
        "prefilter": "Wallis fill (aria) then autoRIFT 5x5 high-pass",
        "geogrid_inputs": "ITS_LIVE autorift_solidearth_0120m parameter tiles",
        "reference_velocity": "ITS_LIVE vx0/vy0 used (no NullReferenceVelocity)",
    },
    "products": {
        "<event>_final_mosaic_sentinel-2_toa_offset_m.tif": (
            "bands dx_m, dy_m, interp_mask, chip_size_x_px; multi-path mosaics: mean of the paths, "
            "nodata 0; single path: a copy, NaN without a nodata tag"
        ),
        "<event>_Path<N>_sentinel-2_toa_autorift_offset_m.tif": "per path, NaN without a nodata tag",
    },
    "known_issues": [
        "glacier reference velocities set the search centres (mis-signed north-south)",
        "season and sun differences between the composites bias offsets by metres",
        "no geolocation-bias correction; results not bit-reproducible (Wallis fill noise)",
    ],
}

V0_COMPOSITE = re.compile(
    r"^(?P<event>.+)_(?P<sensor>S2)_Path(?P<path>\d+)_(?P<level>TOA|SR)_(?P<band>B\w+?)_"
    r"(?P<stage>PRE|POST)_(?P<start>\d{4}-\d{2}-\d{2}_\d{6})UTC_to_(?P<end>\d{4}-\d{2}-\d{2}_\d{6})UTC"
    r"\.tif$"
)
V0_MOSAIC_SUFFIX = "_final_mosaic_sentinel-2_toa_offset_m.tif"
V0_PATH_PRODUCT = re.compile(r"^.+_Path(?P<path>\d+)_sentinel-2_toa_autorift_offset_m\.tif$")

README = """\
# COSEIS optical data store

Composites and autoRIFT products, kept apart by recipe and by product version. Managed with
`python -m aria_coseis.optical.store` (aria_coseis); `info` lists what is here.

- `events.json`: USGS metadata for every event folder name.
- `composites/recipes.json`: how each composite recipe is made.
- `composites/<event>/`: the composites (read-only) and `composites.json`, one record per file
  (recipe, window, checksum, where it came from).
- `products/<version>/`: one folder per product version. `VERSION.json` holds the code, parameters
  and recipes that define it; `catalog.csv` has one row per event (status, paths used and skipped,
  and why).
- `work/<version>/`: intermediates. Safe to delete.

Rules:
- Never edit or overwrite a composite. A new recipe or time window is a new file.
- A product version is made by one code state with one set of parameters; changing either makes a
  new version. Experiments use versions named `exp-*`.
- v0 is the first set of products, registered from `data/GEE_Optical_Downloads` in the production
  clone, which is frozen. The v0 files here are hard links to it.
"""


def store_root(root: str | os.PathLike[str] | None = None) -> Path:
    """The store's root directory (config.OPTICAL_DIR unless given)."""
    return Path(root if root is not None else config.OPTICAL_DIR)


def composite_name(
    sensor: str,
    path: str | int,
    level: str,
    band: str,
    stage: str,
    recipe: str,
    start: datetime,
    end: datetime,
) -> str:
    """File name of a composite, e.g. S2_P121_TOA_B8_PRE_med90_20221108T000000Z_20230206T011734Z.tif."""
    if recipe not in COMPOSITE_RECIPES:
        raise ValueError(f"unknown composite recipe {recipe!r}")
    stamp = "%Y%m%dT%H%M%SZ"
    return (
        f"{sensor}_P{int(path)}_{level}_{band}_{stage}_{recipe}_"
        f"{start.strftime(stamp)}_{end.strftime(stamp)}.tif"
    )


def sha256(path: str | os.PathLike[str], chunk: int = 1 << 24) -> str:
    """SHA-256 hex digest of a file."""
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        while block := f.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path, default: Any) -> Any:
    return json.loads(path.read_text()) if path.exists() else default


def _write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n")
    tmp.replace(path)


def _link_or_copy(src: Path, dst: Path, allow_copy: bool) -> str:
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(src, dst)
        return "hardlink"
    except OSError:
        if not allow_copy:
            raise
        shutil.copy2(src, dst)
        return "copy"


def write_recipes(root: Path) -> None:
    """Write composites/recipes.json; refuse to change the definition of an existing recipe."""
    path = root / "composites" / "recipes.json"
    existing = _read_json(path, {})
    for rid, definition in COMPOSITE_RECIPES.items():
        if rid in existing and existing[rid] != definition:
            raise ValueError(f"recipe {rid!r} already exists with a different definition")
    _write_json(path, {**existing, **COMPOSITE_RECIPES})


def add_composite(
    root: Path,
    event: str,
    src: Path,
    record: dict[str, Any],
    allow_copy: bool = False,
    source: str | None = None,
) -> Path:
    """
    Register a composite file under composites/<event>/ (hard link, or copy if allowed).

    record must hold sensor, path, level, band, stage, recipe, start and end (datetimes). Adding
    the same file again is a no-op; a different file under an existing name raises. source is
    recorded as where the file came from (default: src).
    """
    name = composite_name(
        record["sensor"],
        record["path"],
        record["level"],
        record["band"],
        record["stage"],
        record["recipe"],
        record["start"],
        record["end"],
    )
    folder = root / "composites" / event
    dst = folder / name
    index_path = folder / "composites.json"
    index = _read_json(index_path, {})
    if name in index and dst.exists() and os.path.samefile(src, dst):
        return dst  # already registered (hard link to the same file)
    digest = sha256(src)
    if dst.exists():
        known = index.get(name, {}).get("sha256") or sha256(dst)
        if known != digest:
            raise FileExistsError(f"{dst} exists with different content")
        if name in index:
            return dst
        how = "existing"
    else:
        how = _link_or_copy(src, dst, allow_copy)
    os.chmod(dst, stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH)
    entry = {
        k: (v.strftime("%Y-%m-%dT%H:%M:%SZ") if isinstance(v, datetime) else v)
        for k, v in record.items()
    }
    entry.update(sha256=digest, bytes=dst.stat().st_size, stored_as=how, source=source or str(src))
    entry["registered"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    index[name] = entry
    _write_json(index_path, index)
    return dst


def git_state(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Commit, describe string and dirty flag of the git checkout containing path."""

    def git(*args: str) -> str:
        return subprocess.run(
            ["git", "-C", str(path), *args], capture_output=True, text=True, check=True
        ).stdout.strip()

    try:
        state: dict[str, Any] = {
            "repo": git("rev-parse", "--show-toplevel"),
            "commit": git("rev-parse", "HEAD"),
            "describe": git("describe", "--tags", "--always", "--dirty"),
            "dirty": bool(git("status", "--porcelain", "--untracked-files=no")),
        }
    except (subprocess.CalledProcessError, FileNotFoundError, NotADirectoryError):
        return {"repo": None}
    if state["dirty"]:
        # Identifies the uncommitted changes, so two different dirty states never look the same
        diff = git("diff", "HEAD", "--no-color")
        state["diff_sha256"] = hashlib.sha256(diff.encode()).hexdigest()
    return state


def module_provenance(module: str) -> dict[str, Any]:
    """Where a top-level module is installed and its git state, without importing it."""
    spec = find_spec(module)
    if spec is None or spec.origin is None:
        return {"installed": False}
    location = Path(spec.origin).parent
    out: dict[str, Any] = {"installed": True, "location": str(location)}
    out.update(git_state(location))
    return out


def code_provenance() -> dict[str, Any]:
    """The code a run uses: aria_coseis, hyp3_autorift and autoRIFT, and the Python environment."""
    from importlib.metadata import PackageNotFoundError, version

    out = {name: module_provenance(name) for name in ("aria_coseis", "hyp3_autorift", "autoRIFT")}
    # autoRIFT comes from conda-forge (distribution geo_autoRIFT), not from a git checkout
    for module, dist in (("hyp3_autorift", "hyp3_autorift"), ("autoRIFT", "geo_autoRIFT")):
        try:
            out[module]["distribution_version"] = version(dist)
        except PackageNotFoundError:
            pass
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        packages = sorted((Path(prefix) / "conda-meta").glob("autorift-*.json"))
        if packages:
            out["autoRIFT"]["conda_package"] = packages[0].stem
    out["python"] = sys.version.split()[0]
    out["environment"] = os.environ.get("CONDA_PREFIX") or sys.prefix
    return out


def write_environment_listing(dst: Path) -> None:
    """Save `conda list --explicit` (or `pip freeze`) of the running environment."""
    prefix = os.environ.get("CONDA_PREFIX")
    cmd = ["conda", "list", "--explicit", "--prefix", prefix] if prefix else None
    try:
        text = (
            subprocess.run(cmd, capture_output=True, text=True, check=True).stdout if cmd else ""
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        text = ""
    if not text:
        text = subprocess.run(
            [sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True
        ).stdout
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(text)


def require_store(root: str | os.PathLike[str] | None = None) -> Path:
    """The store root; raises if it has not been set up (avoids creating a second store)."""
    store = store_root(root)
    if not (store / "README.md").exists():
        raise FileNotFoundError(
            f"no optical data store at {store.resolve()}; set COSEIS_OPTICAL_DIR (production: "
            "/u/trappist-r0/colespeed/work/coseis_optical) or pass --store"
        )
    return store


def start_version(
    root: Path, version: str, description: str, definition: dict[str, Any]
) -> dict[str, Any]:
    """
    Create products/<version>/VERSION.json, or check that an existing version has the same
    definition (code, parameters, recipe). Returns the version record.

    Versions not named exp-* must be made from committed code: a dirty checkout in the
    definition's "code" (or, for a remosaicked version, "mosaic_code") raises.
    """
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", version) or version == "v0":
        raise ValueError(f"invalid version name {version!r}")
    if not version.startswith("exp-"):
        dirty = [
            k
            for key in ("code", "mosaic_code")
            for k, v in definition.get(key, {}).items()
            if isinstance(v, dict) and v.get("dirty")
        ]
        if dirty:
            raise ValueError(
                f"version {version} must be made from committed code; uncommitted changes in "
                f"{', '.join(dirty)} (use an exp-* version for experiments)"
            )
    path = root / "products" / version / "VERSION.json"
    definition = json.loads(json.dumps(definition, sort_keys=True))
    existing = _read_json(path, None)
    if existing is not None:
        if existing.get("definition") != definition:
            changed = sorted(
                k
                for k in set(existing.get("definition", {})) | set(definition)
                if existing.get("definition", {}).get(k) != definition.get(k)
            )
            raise ValueError(
                f"version {version} was made with a different {', '.join(changed)}; "
                "use a new version name"
            )
        return existing
    record = {
        "version": version,
        "description": description,
        "created": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "definition": definition,
    }
    _write_json(path, record)
    return record


def write_version(root: Path, definition: dict[str, Any]) -> Path:
    """Create products/<version>/VERSION.json; an existing version must have the same definition."""
    path = root / "products" / definition["version"] / "VERSION.json"
    existing = _read_json(path, None)
    if existing is not None:
        if existing != json.loads(json.dumps(definition, sort_keys=True)):
            raise FileExistsError(
                f"version {definition['version']} exists with another definition"
            )
        return path
    _write_json(path, definition)
    return path


def query_usgs_event(
    origin: datetime, magnitude: float, window_s: int = 10
) -> dict[str, Any] | None:
    """The USGS event closest to an origin time (within window_s) with the given magnitude."""
    params = {
        "format": "geojson",
        "starttime": (origin - timedelta(seconds=window_s)).strftime("%Y-%m-%dT%H:%M:%S"),
        "endtime": (origin + timedelta(seconds=window_s)).strftime("%Y-%m-%dT%H:%M:%S"),
        "minmagnitude": magnitude - 0.15,
        "maxmagnitude": magnitude + 0.15,
    }
    response = requests.get(USGS_EVENT_QUERY, params=params, timeout=60)
    response.raise_for_status()
    features = response.json().get("features", [])
    if not features:
        return None
    t0 = origin.timestamp() * 1000

    def distance(f: dict[str, Any]) -> tuple[float, float]:
        props = f["properties"]
        return (abs(props["time"] - t0), abs((props.get("mag") or 0) - magnitude))

    return _event_record(min(features, key=distance))


def query_usgs_event_by_id(event_id: str) -> dict[str, Any] | None:
    """USGS metadata for an event id (the preferred id of a merged event is returned)."""
    response = requests.get(
        USGS_EVENT_QUERY, params={"format": "geojson", "eventid": event_id}, timeout=60
    )
    if response.status_code == 404:
        return None
    response.raise_for_status()
    return _event_record(response.json())


def _event_record(feature: dict[str, Any]) -> dict[str, Any]:
    props, coords = feature["properties"], feature["geometry"]["coordinates"]
    return {
        "event_id": feature["id"],
        "title": props.get("title"),
        "time_utc": datetime.fromtimestamp(props["time"] / 1000, timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%S.%fZ"
        )[:-4]
        + "Z",
        "magnitude": props.get("mag"),
        "magnitude_type": props.get("magType"),
        "longitude": coords[0],
        "latitude": coords[1],
        "depth_km": coords[2],
        "url": props.get("url"),
    }


def resolve_events(
    root: Path,
    origins: dict[str, tuple[datetime, float]],
    query: Callable[[datetime, float], dict[str, Any] | None] = query_usgs_event,
    known_ids: dict[str, str] | None = None,
    query_by_id: Callable[[str], dict[str, Any] | None] = query_usgs_event_by_id,
) -> dict[str, Any]:
    """
    USGS metadata per event folder, cached in events.json.

    Looked up by origin time and magnitude ({event: (origin, magnitude)}); events without an
    origin are looked up by the event id their manifest recorded (known_ids), if any.
    """
    path = root / "events.json"
    events = _read_json(path, {})
    for name in sorted(set(origins) | set(known_ids or {})):
        if events.get(name, {}).get("event_id"):
            continue
        if name in origins:
            origin, magnitude = origins[name]
            found = query(origin, magnitude)
            events[name] = found or {"event_id": None, "note": "no USGS match for origin time"}
            events[name]["origin_from_composites"] = origin.strftime("%Y-%m-%dT%H:%M:%SZ")
        else:
            found = query_by_id(known_ids[name])
            events[name] = found or {"event_id": None, "note": "manifest event id not found"}
            events[name]["event_id_from_manifest"] = known_ids[name]
    _write_json(path, events)
    return events


def _v0_aria_code(event: str) -> str:
    if event in V0_OLD_SCRIPT_EVENTS:
        return "develop-8b7b7f7"
    if event in V0_UTM_FIX_EVENTS:
        return "aria-971b216"
    return "aria-405ae83"


def _utc(stamp: str) -> datetime:
    return datetime.strptime(stamp, "%Y-%m-%d_%H%M%S").replace(tzinfo=timezone.utc)


def register_event_composites(
    root: Path, event_dir: Path, recipe: str, allow_copy: bool = False
) -> tuple[datetime, float] | None:
    """
    Register the composites of one GEE download folder (<downloads>/<event>/, files named by
    aria_coseis.optical.gee) under a recipe, with the acquisition times its manifest lists.
    Returns the event's (origin time, magnitude), read from the composite names, if any.
    """
    event = event_dir.name
    magnitude = int(event.split("_")[1]) / 10
    origin = None
    for tif in sorted(event_dir.glob("*.tif")):
        m = V0_COMPOSITE.match(tif.name)
        if not m or m["event"] != event:
            continue
        start, end = _utc(m["start"]), _utc(m["end"])
        record = {
            "sensor": m["sensor"],
            "path": int(m["path"]),
            "level": m["level"],
            "band": m["band"],
            "stage": m["stage"],
            "recipe": recipe,
            "start": start,
            "end": end,
        }
        add_composite(root, event, tif, record, allow_copy)
        origin = (end if m["stage"] == "PRE" else start, magnitude)
    manifests = sorted(event_dir.glob("*_autorift_manifest.json"))
    if manifests and origin is not None:
        manifest = json.loads(manifests[0].read_text())
        scene_times = {
            k: manifest[k] for k in ("pre_dates_used", "post_dates_used") if k in manifest
        }
        if scene_times:
            index_path = root / "composites" / event / "composites.json"
            index = _read_json(index_path, {})
            index["_scene_times_all_paths"] = scene_times
            _write_json(index_path, index)
    return origin


def register_downloads(
    downloads: str | os.PathLike[str],
    root: str | os.PathLike[str] | None = None,
    recipe: str = "med90",
    events: list[str] | None = None,
    allow_copy: bool = False,
    resolve: bool = True,
) -> list[str]:
    """
    Register new GEE downloads (<downloads>/<event>/ folders) under a recipe and look the events
    up on USGS. Returns the events registered. For the frozen v0 folder use import_v0.
    """
    store = require_store(root)
    if recipe not in COMPOSITE_RECIPES:
        raise ValueError(f"unknown composite recipe {recipe!r}")
    if "window_days" in COMPOSITE_RECIPES[recipe]:
        raise ValueError(f"{recipe} composites are made by aria_coseis.optical.composites")
    write_recipes(store)
    origins = {}
    for event_dir in sorted(p for p in Path(downloads).iterdir() if p.is_dir()):
        if events and event_dir.name not in events:
            continue
        origin = register_event_composites(store, event_dir, recipe, allow_copy)
        if origin is not None:
            origins[event_dir.name] = origin
    if resolve:
        resolve_events(store, origins)
    return sorted(origins)


def catalogue_origins(eq_info: str | os.PathLike[str]) -> dict[str, tuple[datetime, float]]:
    """{event folder name: (origin time, magnitude)} from an earthquake_info_*.json list."""
    from aria_coseis.utils import to_snake_case

    out = {}
    for e in json.loads(Path(eq_info).read_text()):
        if "time" in e and "magnitude" in e.get("epicenter", {}):
            origin = datetime.strptime(e["time"][:19], "%Y-%m-%d %H:%M:%S")
            out[to_snake_case(e["title"])] = (
                origin.replace(tzinfo=timezone.utc),
                float(e["epicenter"]["magnitude"]),
            )
    return out


def import_v0(
    downloads: str | os.PathLike[str],
    root: str | os.PathLike[str] | None = None,
    allow_copy: bool = False,
    resolve: bool = True,
    query: Callable[[datetime, float], dict[str, Any] | None] = query_usgs_event,
    eq_info: str | os.PathLike[str] | None = None,
    query_by_id: Callable[[str], dict[str, Any] | None] = query_usgs_event_by_id,
) -> list[dict[str, Any]]:
    """
    Register the first products (v0) and their composites from a GEE_Optical_Downloads folder.

    Composites become recipe med90 under composites/<event>/; mosaics and per-path offsets are
    linked into products/v0/<event>/ under their original names. Nothing in downloads is changed
    except that registered composites become read-only. Safe to run again.

    Events are matched to USGS by the origin time in the composite names (PRE ends and POST starts
    at the event). Events without composites use the event id in their manifest, or else an exact
    title match in eq_info, if given.
    """
    downloads = Path(downloads)
    store = store_root(root)
    store.mkdir(parents=True, exist_ok=True)
    if not (store / "README.md").exists():
        (store / "README.md").write_text(README)
    write_recipes(store)
    write_version(store, V0_VERSION)
    rows: list[dict[str, Any]] = []
    origins: dict[str, tuple[datetime, float]] = {}
    manifest_ids: dict[str, str] = {}
    for manifest_path in sorted(downloads.glob("*/*_autorift_manifest.json")):
        event_dir = manifest_path.parent
        event = event_dir.name
        manifest = json.loads(manifest_path.read_text())
        magnitude = int(event.split("_")[1]) / 10
        if manifest.get("event_id"):
            manifest_ids[event] = manifest["event_id"]
        origin = register_event_composites(store, event_dir, "med90", allow_copy)
        if origin is not None:
            origins[event] = origin

        product_dir = store / "products" / "v0" / event
        mosaic = event_dir / f"{event}{V0_MOSAIC_SUFFIX}"
        linked = []
        for src in [mosaic, *sorted(event_dir.glob("*_sentinel-2_toa_autorift_offset_m.tif"))]:
            if src.exists():
                dst = product_dir / src.name
                if not dst.exists():
                    _link_or_copy(src, dst, allow_copy)
                linked.append(src.name)
        skipped = manifest.get("skipped_track_pairs", {})
        rows.append(
            {
                "event": event,
                "magnitude": magnitude,
                "status": manifest.get("status"),
                "paths_used": " ".join(sorted(manifest.get("track_pairs", {}), key=int)),
                "paths_skipped": "; ".join(
                    f"{p}: {v.get('reason', '') if isinstance(v, dict) else v}"
                    for p, v in sorted(skipped.items())
                ),
                "processed_utc": datetime.fromtimestamp(
                    mosaic.stat().st_mtime, timezone.utc
                ).strftime("%Y-%m-%dT%H:%M:%SZ")
                if mosaic.exists()
                else "",
                "aria_code": _v0_aria_code(event) if mosaic.exists() else "",
                "mosaic": f"{event}/{mosaic.name}" if mosaic.exists() else "",
            }
        )
    if eq_info is not None:
        known = catalogue_origins(eq_info)
        for row in rows:
            name = row["event"]
            if name not in origins and name not in manifest_ids and name in known:
                origins[name] = known[name]
    by_id = {k: v for k, v in manifest_ids.items() if k not in origins}
    events = resolve_events(store, origins, query, by_id, query_by_id) if resolve else {}
    for row in rows:
        row["event_id"] = events.get(row["event"], {}).get("event_id") or ""
    write_catalog(store / "products" / "v0" / "catalog.csv", rows)
    return rows


CATALOG_COLUMNS = (
    "event",
    "event_id",
    "magnitude",
    "status",
    "paths_used",
    "paths_skipped",
    "processed_utc",
    "aria_code",
    "mosaic",
)


def update_catalog(path: Path, row: dict[str, Any]) -> None:
    """Add or replace one event's row in a version's catalog.csv."""
    rows = []
    if path.exists():
        with open(path, newline="") as f:
            rows = [r for r in csv.DictReader(f) if r["event"] != row["event"]]
    write_catalog(path, sorted([*rows, row], key=lambda r: r["event"]))


def write_catalog(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CATALOG_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def info(root: str | os.PathLike[str] | None = None) -> str:
    """Summary of the store: recipes, composites per recipe, and product versions."""
    store = store_root(root)
    lines = [f"Store: {store}"]
    recipes = _read_json(store / "composites" / "recipes.json", {})
    counts: dict[str, int] = {}
    for index_path in (store / "composites").glob("*/composites.json"):
        for name, entry in _read_json(index_path, {}).items():
            if not name.startswith("_"):
                counts[entry["recipe"]] = counts.get(entry["recipe"], 0) + 1
    for rid in sorted(recipes):
        lines.append(f"  recipe {rid}: {counts.get(rid, 0)} composites")
    for version_path in sorted((store / "products").glob("*/VERSION.json")):
        version = json.loads(version_path.read_text())
        catalog = version_path.parent / "catalog.csv"
        n = sum(1 for _ in open(catalog)) - 1 if catalog.exists() else 0
        lines.append(f"  version {version['version']}: {n} events - {version['description'][:70]}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="python -m aria_coseis.optical.store", description=__doc__.split("\n")[1]
    )
    parser.add_argument("--store", default=None, help="Store root (default: COSEIS_OPTICAL_DIR).")
    sub = parser.add_subparsers(dest="command", required=True)
    imp = sub.add_parser("import-v0", help="Register the v0 products and composites.")
    imp.add_argument("downloads", type=Path, help="The GEE_Optical_Downloads folder.")
    imp.add_argument("--copy", action="store_true", help="Copy if hard links are not possible.")
    imp.add_argument("--no-resolve", action="store_true", help="Skip the USGS event lookup.")
    imp.add_argument(
        "--eq_info", type=Path, help="earthquake_info_*.json for events without composites."
    )
    reg = sub.add_parser("register", help="Register new GEE downloads under a recipe.")
    reg.add_argument("downloads", type=Path, help="Folder of <event>/ download folders.")
    reg.add_argument("--recipe", default="med90", help="Composite recipe (default med90).")
    reg.add_argument("--event", action="append", help="Only this event folder (repeatable).")
    reg.add_argument("--copy", action="store_true", help="Copy if hard links are not possible.")
    sub.add_parser("info", help="List recipes, composites and product versions.")
    args = parser.parse_args(argv)
    if args.command == "import-v0":
        rows = import_v0(
            args.downloads,
            args.store,
            allow_copy=args.copy,
            resolve=not args.no_resolve,
            eq_info=args.eq_info,
        )
        missing = [r["event"] for r in rows if not r["event_id"]]
        print(f"Registered {len(rows)} events; without a USGS event id: {len(missing)}")
        for name in missing:
            print(f"  {name}")
    elif args.command == "register":
        if (args.downloads / "README_FROZEN_v0.md").exists():
            raise SystemExit(f"{args.downloads} is the frozen v0 folder; use import-v0")
        done = register_downloads(
            args.downloads, args.store, args.recipe, args.event, allow_copy=args.copy
        )
        print(f"Registered composites of {len(done)} event(s) under recipe {args.recipe}")
    print(info(args.store))


if __name__ == "__main__":
    main()
