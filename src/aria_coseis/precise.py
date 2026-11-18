"""
Production products from precise orbits, and the processing machine's email outbox.

Forward processing runs as soon as the post-event SLCs exist, when only restituted orbits are
available, so its products are urgent-response products (URGENT_RESPONSE/, S1-COSEIS_SAR_UR-...).
Each one is queued here (one JSON per product in config.precise_queue_dir()). The precise runner
(`python -m aria_coseis --precise`, ops/run_coseis_precise.sh) works through the queue: once
precise orbits (AUX_POEORB, about 3 weeks after acquisition) exist for every SLC of the job and the
event meets the historic significance criteria, it reruns the same job (topsApp takes precise
orbits when they exist) and places the result in HISTORIC_EVENTS/ under the production name.

Queue statuses: pending -> done | not_eligible | failed. A failed entry is left for a human: fix
the cause and set its status back to "pending" to retry.

The processing machine cannot send email: the precise runner writes notifications to
config.OUTBOX_STAGING_DIR, ops/run_coseis_precise.sh moves them to config.OUTBOX_DIR (in git) and
pushes, and the GitHub Action emails them (send_outbox) and deletes them.
"""

from __future__ import annotations

import glob
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

from aria_coseis import config
from aria_coseis.notify import send_email
from aria_coseis.sar.products import (
    PRODUCTION_SET,
    URGENT_SET,
    all_precise,
    archive_rel_dir,
    clean_processing_dir,
    event_folder,
    event_year,
    find_product,
    orbits_used,
    place_product,
    urgent_name,
)
from aria_coseis.sar.topsapp import run_dockerized_topsApp
from aria_coseis.significance import check_significance
from aria_coseis.usgs import parse_geojson
from aria_coseis.utils import convert_time

# ASF's orbit service: redirects to the orbit file it would serve for a scene (precise if it
# exists, else restituted), the same lookup DockerizedTopsApp makes (s1_orbits.fetch_for_scene)
ORBIT_API = "https://s1-orbits.asf.alaska.edu/scene"


def now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_json(path: str, data: Any) -> None:
    """Write JSON atomically (a hidden temporary file, then rename)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = os.path.join(os.path.dirname(path), f".{os.path.basename(path)}.tmp")
    with open(tmp, "w") as f:
        json.dump(data, f, indent=4)
    os.replace(tmp, path)


def orbit_type(scene: str) -> str | None:
    """'POEORB' (precise) or 'RESORB' (restituted) for the orbit ASF would serve; None if none."""
    res = requests.get(f"{ORBIT_API}/{scene}", allow_redirects=False, timeout=60)
    if res.status_code == 404:
        return None
    if res.is_redirect:
        target = res.headers.get("Location", "")
    else:
        res.raise_for_status()
        target = res.url
    for kind in ("POEORB", "RESORB"):
        if f"AUX_{kind}" in target:
            return kind
    return None


def meets_criteria(event_id: str) -> bool | None:
    """
    Whether the event meets the historic significance criteria, with its current USGS parameters.
    None if that cannot be decided now (USGS or the coastline service unreachable).
    """
    try:
        res = requests.get(
            config.USGS_api_alltime,
            params={"eventid": event_id, "format": "geojson"},
            timeout=60,
        )
        if res.status_code in (404, 409):  # deleted or merged into another event
            print(f"    USGS has no event {event_id} (HTTP {res.status_code}).")
            return False
        res.raise_for_status()
        earthquakes = parse_geojson({"features": [res.json()]})
        date = convert_time(earthquakes[0]["time"]).strftime("%Y-%m-%d")
        return bool(check_significance(earthquakes, date, mode="historic"))
    except (requests.RequestException, AttributeError, KeyError, ValueError) as e:
        print(f"    Could not check the significance of {event_id}: {e}")
        return None


def queue_urgent_product(
    nc_file: str,
    processing_dir: str,
    event_id: str,
    event_data: dict[str, Any],
    track_key: str,
    job: dict[str, Any],
    dates: str,
) -> str:
    """
    Place a forward-mode product in URGENT_RESPONSE/ under its urgent name and queue it for a
    precise-orbit rerun. Returns the urgent product's folder.
    """
    name = Path(nc_file).stem
    year = event_year(event_data["time"])
    folder = event_folder(event_data["title"], event_id)
    urgent_rel = archive_rel_dir(URGENT_SET, year, folder, urgent_name(name))
    production_rel = archive_rel_dir(PRODUCTION_SET, year, folder, name)

    urgent_dir = place_product(
        nc_file, os.path.join(config.ARCHIVE_DIR, urgent_rel), urgent_name(name)
    )
    entry = {
        "status": "pending",
        "event_id": event_id,
        "title": event_data["title"],
        "time": event_data["time"],
        "track": track_key,
        "dates": dates,
        "urgent_product": urgent_rel,
        "production_product": production_rel,
        "urgent_orbits": orbits_used(processing_dir),
        "job": job,
        "log": [f"{now_utc()} queued"],
    }
    write_json(os.path.join(config.precise_queue_dir(), f"{urgent_name(name)}.json"), entry)
    return urgent_dir


def load_queue() -> list[tuple[str, dict[str, Any]]]:
    """(path, entry) for every queued product, oldest name first."""
    queue = []
    for path in sorted(glob.glob(os.path.join(config.precise_queue_dir(), "*.json"))):
        with open(path) as f:
            queue.append((path, json.load(f)))
    return queue


def notify(kind: str, entry: dict[str, Any], location: str, message: str) -> None:
    """Leave an email notification for the GitHub Action (see the module docstring)."""
    note = {
        "kind": kind,
        "title": entry["title"],
        "track": entry["track"],
        "dates": entry["dates"],
        "location": location,
        "message": message,
        "time": now_utc(),
    }
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    name = f"{stamp}_{os.path.basename(entry['urgent_product'])}.json"
    write_json(os.path.join(config.OUTBOX_STAGING_DIR, name), note)


def set_status(entry: dict[str, Any], status: str, message: str) -> None:
    entry["status"] = status
    entry["log"].append(f"{now_utc()} {status}: {message}")
    print(f"    {status}: {message}")


def make_production(entry: dict[str, Any]) -> None:
    """Rerun the entry's job with precise orbits and place the product in HISTORIC_EVENTS/."""
    job = entry["job"]
    params = job["job_parameters"]
    production_dir = os.path.join(config.ARCHIVE_DIR, entry["production_product"])
    expected = os.path.basename(entry["production_product"])

    # The urgent run may already have had precise orbits (processed late): same product
    if all_precise(entry["urgent_orbits"]):
        urgent_dir = os.path.join(config.ARCHIVE_DIR, entry["urgent_product"])
        urgent_nc = os.path.join(urgent_dir, os.path.basename(urgent_dir) + ".nc")
        place_product(urgent_nc, production_dir, expected)
        set_status(entry, "done", "urgent product already had precise orbits; copied")
        return

    # Scratch folder named like forward mode's: <title>/<track>/precise/<track>_<older>_<newer>
    flight_dir, path_num = entry["track"].split("_")
    track = f"{flight_dir}{int(path_num):03d}"
    newer, older = expected.split("-")[6].split("_")
    work_dir = os.path.join(
        config.root_dir, entry["title"], track, "precise", f"{track}_{older}_{newer}"
    )
    print(f"    Rerunning with precise orbits in {work_dir}")
    run_dockerized_topsApp(job, work_dir, cores=config.PRECISE_CORES)

    nc_file = find_product(work_dir)
    if nc_file is None:
        raise RuntimeError(f"topsApp finished but wrote no product in {work_dir}")
    if Path(nc_file).stem != expected:
        raise RuntimeError(f"product {Path(nc_file).stem} is not the expected {expected}")
    orbits = orbits_used(work_dir)
    if not all_precise(orbits):
        raise RuntimeError(f"topsApp did not use precise orbits for every SLC: {orbits}")
    clean_processing_dir(work_dir)
    place_product(nc_file, production_dir, expected)
    set_status(
        entry,
        "done",
        f"placed {expected} ({len(params['granules'])} + "
        f"{len(params['secondary_granules'])} SLCs, orbits {', '.join(orbits)})",
    )


def process_entry(entry: dict[str, Any]) -> bool:
    """Advance one pending entry. Returns True if its status changed."""
    production_dir = os.path.join(config.ARCHIVE_DIR, entry["production_product"])
    if os.path.exists(production_dir):
        set_status(entry, "done", "production product already in the archive")
        return True

    params = entry["job"]["job_parameters"]
    scenes = params["granules"] + params["secondary_granules"]
    if not all_precise(entry["urgent_orbits"]):
        try:
            types = {scene: orbit_type(scene) for scene in scenes}
        except requests.RequestException as e:
            print(f"    Could not check the orbits (trying again next run): {e}")
            return False
        missing = [scene for scene, kind in types.items() if kind != "POEORB"]
        if missing:
            queued = datetime.strptime(entry["log"][0].split()[0], "%Y-%m-%dT%H:%M:%SZ")
            waited = datetime.now(timezone.utc) - queued.replace(tzinfo=timezone.utc)
            if waited.days >= config.PRECISE_GIVE_UP_DAYS:
                message = (
                    f"no precise orbits after {waited.days} days for {len(missing)} of "
                    f"{len(scenes)} SLCs"
                )
                set_status(entry, "failed", message)
                notify("failed", entry, production_dir, message)
                return True
            print(f"    Waiting for precise orbits ({len(missing)} of {len(scenes)} SLCs).")
            return False

    eligible = meets_criteria(entry["event_id"])
    if eligible is None:
        return False
    if not eligible:
        set_status(entry, "not_eligible", "event does not meet the significance criteria")
        return True

    try:
        make_production(entry)
    except Exception as e:
        set_status(entry, "failed", str(e))
        notify("failed", entry, production_dir, str(e))
        return True
    notify("production", entry, production_dir, entry["log"][-1].split(" ", 1)[1])
    return True


def check_precise_queue() -> None:
    """Work through the queue of urgent-response products (one precise-runner run)."""
    queue = [(path, entry) for path, entry in load_queue() if entry["status"] == "pending"]
    print(f"{len(queue)} urgent-response products await a precise-orbit rerun.")
    for path, entry in queue:
        print(f"  {os.path.basename(entry['urgent_product'])}")
        try:
            changed = process_entry(entry)
        except Exception as e:  # one bad entry must not stop the others
            print(f"    Unexpected error (trying again next run): {type(e).__name__}: {e}")
            continue
        if changed:
            write_json(path, entry)


def outbox_ready() -> bool:
    """True if outbox notifications can be emailed (recipients and credentials configured)."""
    return bool(
        config.SECONDARY_RECIPIENTS
        and os.environ.get("GMAIL_USER")
        and os.environ.get("GMAIL_APP_PSWD")
    )


def send_outbox() -> None:
    """GitHub Action: email the processing machine's notifications in one message, then delete them."""
    files = sorted(glob.glob(os.path.join(config.OUTBOX_DIR, "*.json")))
    if not files:
        return
    if not outbox_ready():
        print(f"WARNING: {len(files)} outbox notifications kept: email is not configured.")
        return
    notes = []
    for path in files:
        with open(path) as f:
            notes.append(json.load(f))

    placed = sum(n["kind"] == "production" for n in notes)
    failed = len(notes) - placed
    subject = f"PRODUCTION PRODUCTS (PRECISE ORBITS): {placed} placed"
    if failed:
        subject += f", {failed} need attention"

    body = """
    <html>
        <body style="font-family: Arial, sans-serif; color: #333; line-height: 1.6; margin: 0; padding: 20px;">
        <div style="max-width: 600px; margin: 0 auto; border: 1px solid #e0e0e0; padding: 20px; border-radius: 8px;">
            <h2 style="color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; margin-top: 0;">Production Products Update</h2>
            <p>Urgent-response products rerun with precise orbits for the production archive:</p>
    """
    for n in notes:
        ok = n["kind"] == "production"
        color = "#27ae60" if ok else "#e74c3c"
        status = "Placed" if ok else "Needs attention"
        body += f"""
            <div style="background-color: #f8f9fa; padding: 15px; margin-bottom: 15px; border-left: 5px solid {color}; border-radius: 4px;">
            <p style="margin: 0 0 5px;"><strong>Event:</strong> {n["title"]}</p>
            <p style="margin: 0 0 5px;"><strong>Track:</strong> {n["track"]}</p>
            <p style="margin: 0 0 5px;"><strong>Dates:</strong> {n["dates"]}</p>
            <p style="margin: 0 0 5px;"><strong>Status:</strong> <span style="color: {color}; font-weight: bold;">{status}</span>: {n["message"]}</p>
            <p style="margin: 0;"><strong>Location:</strong> <br>
                <code style="background: #e9ecef; padding: 4px; display: block; margin-top: 5px; word-wrap: break-word; font-size: 12px; color: #c0392b;">{n["location"]}</code>
            </p>
            </div>
        """
    body += """
            <p style="font-size: 12px; color: #95a5a6; border-top: 1px solid #e0e0e0; padding-top: 15px; margin-top: 20px; text-align: center;">
            This is an automated message from the COSEIS pipeline.
            </p>
        </div>
        </body>
    </html>
    """
    print("Sending production-products email to secondary recipients...")
    send_email(subject, body, recipients=config.SECONDARY_RECIPIENTS)
    for path in files:
        os.remove(path)
