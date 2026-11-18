"""
COSEIS_SAR products in the archive: names, archive paths, orbit types, and placement.

Forward mode makes each product twice, from the same job (the same SLCs):
- urgent response: right after the post-event SLCs arrive, with restituted orbits, named
  S1-COSEIS_SAR_UR-... and placed in URGENT_RESPONSE/;
- production: once precise orbits exist (about 3 weeks later), named S1-COSEIS_SAR-... as
  DockerizedTopsApp makes it, and placed in HISTORIC_EVENTS/.
The same SLCs give the same dates, corner and scene hash, so the two names differ only in the
product type. Archive layout: <ARCHIVE_DIR>/<set>/<year>/<event folder>/<product>/<product>.{nc,png}
"""

from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from aria_coseis import config

URGENT_SET = "URGENT_RESPONSE"
PRODUCTION_SET = "HISTORIC_EVENTS"
PRODUCTION_PREFIX = "S1-COSEIS_SAR-"
URGENT_PREFIX = "S1-COSEIS_SAR_UR-"

# Run as a script in the topsApp environment (it needs h5py; see its docstring)
REPACK_SCRIPT = Path(__file__).with_name("repack.py")

# ISCE2 working folders deleted once a product exists (raw SLCs are deleted too)
INTERMEDIATE_DIRS = (
    "geom_reference",
    "ion",
    "fine_interferogram",
    "fine_offsets",
    "fine_coreg",
    "mask",
    "reference",
    "secondary",
    "PICKLE",
    "aux_cal",
    "orbits",
)


def urgent_name(production_name: str) -> str:
    """S1-COSEIS_SAR-... -> S1-COSEIS_SAR_UR-..."""
    if not production_name.startswith(PRODUCTION_PREFIX):
        raise ValueError(f"Not a COSEIS_SAR product name: {production_name}")
    return URGENT_PREFIX + production_name[len(PRODUCTION_PREFIX) :]


def production_name(urgent: str) -> str:
    """S1-COSEIS_SAR_UR-... -> S1-COSEIS_SAR-..."""
    if not urgent.startswith(URGENT_PREFIX):
        raise ValueError(f"Not an urgent-response product name: {urgent}")
    return PRODUCTION_PREFIX + urgent[len(URGENT_PREFIX) :]


def event_folder(title: str, event_id: str) -> str:
    """The event's folder name: its snake-case title and USGS event ID."""
    return f"{title}_{event_id}"


def event_year(event_time: int | float) -> str:
    """Year of the event's origin time (USGS milliseconds since the epoch, UTC)."""
    return str(datetime.fromtimestamp(event_time / 1000, tz=timezone.utc).year)


def archive_rel_dir(product_set: str, year: str, folder: str, name: str) -> str:
    """A product's folder relative to the archive root (queue entries store these)."""
    return os.path.join(product_set, year, folder, name)


def find_product(processing_dir: str) -> str | None:
    """The COSEIS_SAR netCDF that topsApp wrote in processing_dir (in <name>/<name>.nc)."""
    found = sorted(glob.glob(os.path.join(processing_dir, f"{PRODUCTION_PREFIX}*", "*.nc")))
    found = [p for p in found if Path(p).stem == Path(p).parent.name]
    if len(found) > 1:
        raise RuntimeError(f"Several products in {processing_dir}: {found}")
    return found[0] if found else None


def orbits_used(processing_dir: str) -> list[str]:
    """Orbit files topsApp used (from its loc_data.json, which outlives the orbits/ folder)."""
    path = os.path.join(processing_dir, "loc_data.json")
    if not os.path.exists(path):
        return []
    with open(path) as f:
        loc = json.load(f)
    files = loc.get("reference_orbits", []) + loc.get("secondary_orbits", [])
    return sorted({os.path.basename(p) for p in files})


def all_precise(orbit_files: list[str]) -> bool:
    """True if there are orbit files and all are precise (AUX_POEORB)."""
    return bool(orbit_files) and all("_AUX_POEORB_" in f for f in orbit_files)


def place_file(src: str, dst: str) -> None:
    """
    Put src at dst in the archive: a netCDF is repacked with its extra layers compressed, other
    files are copied; either is checked against src before it is renamed into place (repack.py).
    """
    cmd = ["conda", "run", "-n", config.TOPSAPP_ENV, "python", "-I", str(REPACK_SCRIPT), src, dst]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        detail = (result.stdout.strip().splitlines() or [result.stderr.strip()])[-1]
        raise RuntimeError(f"Placing {os.path.basename(src)} failed: {detail}")


def place_product(src_nc: str, archive_dir: str, name: str) -> str:
    """
    Place a product (its .nc and, if present, .png) in archive_dir/<name>/ under the given name.
    Returns the product folder. Refuses to overwrite a product that is already there.
    """
    if os.path.exists(archive_dir):
        raise FileExistsError(f"{archive_dir} already exists")
    for ext in (".nc", ".png"):
        src = os.path.splitext(src_nc)[0] + ext
        if ext == ".png" and not os.path.exists(src):
            continue
        place_file(src, os.path.join(archive_dir, name + ext))
    return archive_dir


def clean_processing_dir(processing_dir: str) -> None:
    """Delete raw SLCs and ISCE2 intermediate folders once the product exists."""
    for slc_zip in glob.glob(os.path.join(processing_dir, "S1[A-D]*.zip")):
        os.remove(slc_zip)
    for slc_safe in glob.glob(os.path.join(processing_dir, "S1[A-D]*.SAFE")):
        shutil.rmtree(slc_safe, ignore_errors=True)
    for idir in INTERMEDIATE_DIRS:
        dir_path = os.path.join(processing_dir, idir)
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path, ignore_errors=True)
