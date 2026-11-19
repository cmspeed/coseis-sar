"""
Production products from precise orbits, and the processing machine's email outbox.

Forward processing runs as soon as the post-event SLCs exist, when only restituted orbits are
available, so its products are urgent-response products (URGENT_RESPONSE/, S1-COSEIS_SAR_UR-...).
Each one is queued here (one JSON per urgent product in config.precise_queue_dir()). The precise
runner (`python -m aria_coseis --precise`, ops/run_coseis_precise.sh) turns the queue into production
products for HISTORIC_EVENTS/, made the way historic mode makes them:
- the owner is the largest event (ties: the later one) among the events that used the urgent product
  and meet the historic significance criteria, with their current USGS parameters;
- the pair is the one historic mode selects for the owner on that track (its AOI, coverage rule and
  SLC search), so it can differ from the urgent product's: forward mode picks the pre-event date
  when the event is found and takes the first post-event pass;
- it runs once precise orbits (AUX_POEORB, about 3 weeks after acquisition) exist for every SLC,
  and a production job that another queue entry already produced is not run again.

Pairs shared by several events (foreshocks, aftershocks, sequences) follow the archive's rules: one
urgent copy per product name, in the folder of the largest event that used it; one production
product per pair, owned by the largest eligible event.

Queue statuses: pending -> done | not_eligible | failed (not_eligible goes back to pending if another
event joins). A failed entry is left for a human: fix the cause and set its status back to
"pending" to retry. Entries with any other status (e.g. "backfill", for urgent products placed by
hand whose production products were made by a historic run) are left alone.

The processing machine cannot send email: the precise runner writes notifications to
config.OUTBOX_STAGING_DIR, ops/run_coseis_precise.sh moves them to config.OUTBOX_DIR (in git) and
pushes, and the GitHub Action emails them (send_outbox) and deletes them.
"""

from __future__ import annotations

import contextlib
import fcntl
import glob
import json
import os
import shutil
from collections.abc import Callable, Iterator
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
    production_name,
    urgent_name,
)
from aria_coseis.sar.topsapp import run_dockerized_topsApp
from aria_coseis.significance import check_significance
from aria_coseis.usgs import parse_geojson
from aria_coseis.utils import to_snake_case

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


@contextlib.contextmanager
def queue_lock() -> Iterator[None]:
    """Both runners read and rewrite queue entries: one at a time (a lock file on local disk)."""
    with open(config.QUEUE_LOCK_FILE, "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def update_entry(path: str, change: Callable[[dict[str, Any]], None]) -> dict[str, Any]:
    """Apply change to the entry as it is on disk now, under the queue lock, and save it."""
    with queue_lock():
        with open(path) as f:
            entry = json.load(f)
        change(entry)
        write_json(path, entry)
    return entry


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


def event_info(event_id: str) -> dict[str, Any] | None:
    """
    The event as USGS describes it now (parse_geojson's format), plus "historic": whether it meets
    the historic significance criteria. {} if USGS has no such event (deleted or merged); None if
    that cannot be decided now (USGS or the coastline service unreachable).
    """
    try:
        res = requests.get(
            config.USGS_api_alltime,
            params={"eventid": event_id, "format": "geojson"},
            timeout=60,
        )
        if res.status_code in (404, 409):
            print(f"    USGS has no event {event_id} (HTTP {res.status_code}).")
            return {}
        res.raise_for_status()
        (eq,) = parse_geojson({"features": [res.json()]})
        date = datetime.fromtimestamp(eq["time"] / 1000, tz=timezone.utc).strftime("%Y-%m-%d")
        eq["historic"] = bool(check_significance([dict(eq)], date, mode="historic"))
        return eq
    except (requests.RequestException, AttributeError, KeyError, ValueError) as e:
        print(f"    Could not check event {event_id}: {e}")
        return None


def outranks(a: dict[str, Any], b: dict[str, Any]) -> bool:
    """True if event a owns over b: a larger magnitude, or the same magnitude and later."""
    return (a.get("mag") or 0, a["time"]) > (b.get("mag") or 0, b["time"])


def largest(events: list[dict[str, Any]]) -> dict[str, Any]:
    owner = events[0]
    for e in events[1:]:
        if outranks(e, owner):
            owner = e
    return owner


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
    production product. If another event already queued the same product (the same SLCs), it is not
    placed again: this event is added to the entry, and the copy moves to this event's folder if this
    event outranks the others. Returns the urgent product's folder.
    """
    urgent = urgent_name(Path(nc_file).stem)
    queue_file = os.path.join(config.precise_queue_dir(), f"{urgent}.json")
    info = event_info(event_id) or {}
    event = {
        "event_id": event_id,
        "title": event_data["title"],
        "time": event_data["time"],
        "mag": info.get("mag"),
    }

    def folder_of(e: dict[str, Any]) -> str:
        return archive_rel_dir(
            URGENT_SET, event_year(e["time"]), event_folder(e["title"], e["event_id"]), urgent
        )

    with queue_lock():
        if not os.path.exists(queue_file):
            rel = folder_of(event)
            place_product(nc_file, os.path.join(config.ARCHIVE_DIR, rel), urgent)
            entry = {
                "status": "pending",
                "urgent_product": rel,
                "events": [event],
                "track": track_key,
                "dates": dates,
                "urgent_orbits": orbits_used(processing_dir),
                "job": job,
                "production": {},
                "log": [f"{now_utc()} queued for {event_id}"],
            }
            write_json(queue_file, entry)
            return os.path.join(config.ARCHIVE_DIR, rel)

        with open(queue_file) as f:
            entry = json.load(f)
        if event_id not in {e["event_id"] for e in entry["events"]}:
            entry["events"].append(event)
            entry["log"].append(f"{now_utc()} also used by {event_id}")
            if entry["status"] == "not_eligible":
                entry["status"] = "pending"  # the new event may be significant
        rel = folder_of(largest(entry["events"]))
        if rel != entry["urgent_product"]:
            # The event that outranks the others keeps the copy (the archive's ownership rule)
            src = os.path.join(config.ARCHIVE_DIR, entry["urgent_product"])
            dst = os.path.join(config.ARCHIVE_DIR, rel)
            if os.path.exists(dst):
                raise FileExistsError(f"{dst} already exists")
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            os.rename(src, dst)
            with contextlib.suppress(OSError):
                os.rmdir(os.path.dirname(src))  # the other event's folder, if now empty
            entry["log"].append(f"{now_utc()} moved to {os.path.dirname(rel)}")
            entry["urgent_product"] = rel
        write_json(queue_file, entry)
        return os.path.join(config.ARCHIVE_DIR, rel)


def load_queue() -> list[tuple[str, dict[str, Any]]]:
    """(path, entry) for every queued product, oldest name first."""
    queue = []
    for path in sorted(glob.glob(os.path.join(config.precise_queue_dir(), "*.json"))):
        with open(path) as f:
            queue.append((path, json.load(f)))
    return queue


def notify(kind: str, entry: dict[str, Any], location: str, message: str) -> None:
    """Leave an email notification for the GitHub Action (see the module docstring)."""
    production = entry.get("production", {})
    note = {
        "kind": kind,
        "title": production.get("owner_title") or entry["events"][0]["title"],
        "track": entry["track"],
        "dates": production.get("dates") or entry["dates"],
        "location": location,
        "message": message,
        "time": now_utc(),
    }
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
    name = f"{stamp}_{os.path.basename(entry['urgent_product'])}.json"
    write_json(os.path.join(config.OUTBOX_STAGING_DIR, name), note)


def log(entry: dict[str, Any], message: str) -> None:
    entry["log"].append(f"{now_utc()} {message}")
    print(f"    {message}")


def set_status(entry: dict[str, Any], status: str, message: str) -> None:
    entry["status"] = status
    log(entry, f"{status}: {message}")


def fail(entry: dict[str, Any], location: str, message: str) -> None:
    set_status(entry, "failed", message)
    notify("failed", entry, location, message)


def scene_key(job: dict[str, Any]) -> list[list[str]]:
    p = job["job_parameters"]
    return [sorted(p["granules"]), sorted(p["secondary_granules"])]


def pair_dates(job: dict[str, Any]) -> tuple[str, str]:
    p = job["job_parameters"]
    return min(s[17:25] for s in p["secondary_granules"]), min(s[17:25] for s in p["granules"])


def contains(job: dict[str, Any], event_time: int | float) -> bool:
    """True if the event's origin falls between the pair's pre- and post-event acquisitions."""
    p = job["job_parameters"]
    origin = datetime.fromtimestamp(event_time / 1000, tz=timezone.utc).strftime("%Y%m%dT%H%M%S")
    return (
        max(s[17:32] for s in p["secondary_granules"])
        < origin
        < min(s[17:32] for s in p["granules"])
    )


def historic_job(owner: dict[str, Any], track_key: str) -> dict[str, Any] | None:
    """
    The coseismic job historic mode makes for the owner event on this track (30 m, as the archive),
    or None if it makes none (the track misses the owner's AOI, or there is no pair yet).
    """
    from aria_coseis.pipeline import event_aoi
    from aria_coseis.sar.pairing import find_reference_and_secondary_pairs
    from aria_coseis.sar.search import get_path_and_frame_numbers, get_SLCs

    flight_direction, path_number = track_key.split("_")
    aoi = event_aoi(owner)
    tracks, _ = get_path_and_frame_numbers(aoi, owner["time"])
    if not any(fd == flight_direction and int(pn) == int(path_number) for fd, pn in tracks):
        return None
    slcs = get_SLCs(
        flight_direction, path_number, aoi.wkt, owner["time"], processing_mode="historic"
    )
    jobs = find_reference_and_secondary_pairs(
        slcs,
        owner["time"],
        flight_direction,
        path_number,
        to_snake_case(owner["title"]),
        aoi,
        owner["id"],
        "coseismic",
        job_list=True,
        resolution=30,
    )
    return jobs[0] if jobs else None


def wait_or_give_up(entry: dict[str, Any], what: str) -> None:
    """Keep waiting for what, or flag the entry once it has waited too long."""
    queued = datetime.strptime(entry["log"][0].split()[0], "%Y-%m-%dT%H:%M:%SZ")
    days = (datetime.now(timezone.utc) - queued.replace(tzinfo=timezone.utc)).days
    if days >= config.PRECISE_GIVE_UP_DAYS:
        fail(entry, entry["urgent_product"], f"{what} after {days} days")
    else:
        print(f"    Waiting: {what}.")


def make_production(entry: dict[str, Any], job: dict[str, Any], production_dir: str) -> str:
    """
    Make the production product for job and place it in production_dir (relative to the archive).
    The urgent product is copied instead if it is the same pair and already had precise orbits.
    Returns the product's folder relative to the archive.
    """
    if scene_key(job) == scene_key(entry["job"]) and all_precise(entry["urgent_orbits"]):
        urgent_dir = os.path.join(config.ARCHIVE_DIR, entry["urgent_product"])
        name = production_name(os.path.basename(urgent_dir))
        rel = os.path.join(production_dir, name)
        place_product(
            os.path.join(urgent_dir, os.path.basename(urgent_dir) + ".nc"),
            os.path.join(config.ARCHIVE_DIR, rel),
            name,
        )
        log(entry, f"copied the urgent product (same pair, precise orbits) as {name}")
        return rel

    # Scratch folder named like historic mode's: <title>/<track>/precise/<track>_<pre>_<post>
    folder, track = job["name"].rsplit("-", 1)
    pre, post = pair_dates(job)
    work_dir = os.path.join(config.root_dir, folder, track, "precise", f"{track}_{pre}_{post}")
    if find_product(work_dir) is None and os.path.exists(work_dir):
        shutil.rmtree(work_dir)  # an earlier attempt that left no product
    print(f"    Running with precise orbits in {work_dir}")
    run_dockerized_topsApp(job, work_dir, cores=config.PRECISE_CORES)

    nc_file = find_product(work_dir)
    if nc_file is None:
        raise RuntimeError(f"topsApp finished but wrote no product in {work_dir}")
    orbits = orbits_used(work_dir)
    if not all_precise(orbits):
        raise RuntimeError(f"topsApp did not use precise orbits for every SLC: {orbits}")
    clean_processing_dir(work_dir)
    name = Path(nc_file).stem
    rel = os.path.join(production_dir, name)
    place_product(nc_file, os.path.join(config.ARCHIVE_DIR, rel), name)
    log(entry, f"placed {name} (orbits {', '.join(orbits)})")
    return rel


def process_entry(entry: dict[str, Any], others: list[dict[str, Any]]) -> None:
    """Advance one pending entry (others: the rest of the queue, as on disk)."""
    production = entry.setdefault("production", {})

    # Owner: the largest event the pair contains that meets the historic criteria. Contained: the
    # events that used this urgent product, and events queued on the same track (with other
    # frames) whose origin falls between the pair's acquisitions
    events = {e["event_id"] for e in entry["events"]}
    for other in others:
        if other["track"] == entry["track"]:
            events |= {e["event_id"] for e in other["events"] if contains(entry["job"], e["time"])}
    infos = []
    for e in sorted(events):
        info = event_info(e)
        if info is None:
            return  # try again next run
        if info and info["historic"]:
            infos.append(info)
    if not infos:
        set_status(entry, "not_eligible", "no event meets the significance criteria")
        return
    owner = largest(infos)
    if production.get("owner") not in (None, owner["id"]):
        production.pop("job", None)  # the owner changed: select the pair again
    production.update(owner=owner["id"], owner_title=to_snake_case(owner["title"]))

    # The pair historic mode selects for the owner on this track
    if "job" not in production:
        try:
            job = historic_job(owner, entry["track"])
        except requests.RequestException as e:
            print(f"    Could not search ASF (trying again next run): {e}")
            return
        if job is None:
            wait_or_give_up(entry, "no historic-mode pair for the owner on this track")
            return
        production["job"] = job
        log(
            entry,
            f"historic-mode pair for {owner['id']}: {job['name']} {'_'.join(pair_dates(job))}",
        )
    job = production["job"]
    pre, post = pair_dates(job)
    production["dates"] = f"{pre} (Pre) - {post} (Post)"

    # Made already through another urgent product (the same pair of the same owner)
    for other in others:
        o = other.get("production", {})
        if o.get("product") and o.get("job") and scene_key(o["job"]) == scene_key(job):
            production["product"] = o["product"]
            set_status(entry, "done", f"same production product as {other['urgent_product']}")
            return

    # Precise orbits for every SLC of the production job
    if not (scene_key(job) == scene_key(entry["job"]) and all_precise(entry["urgent_orbits"])):
        p = job["job_parameters"]
        scenes = p["granules"] + p["secondary_granules"]
        try:
            missing = [s for s in scenes if orbit_type(s) != "POEORB"]
        except requests.RequestException as e:
            print(f"    Could not check the orbits (trying again next run): {e}")
            return
        if missing:
            wait_or_give_up(entry, f"no precise orbits for {len(missing)} of {len(scenes)} SLCs")
            return

    production_dir = os.path.join(
        PRODUCTION_SET,
        event_year(owner["time"]),
        event_folder(production["owner_title"], owner["id"]),
    )
    try:
        production["product"] = make_production(entry, job, production_dir)
    except Exception as e:
        fail(entry, os.path.join(config.ARCHIVE_DIR, production_dir), str(e))
        return
    message = entry["log"][-1].split(" ", 1)[1]
    set_status(entry, "done", f"production product {production['product']}")
    notify("production", entry, os.path.join(config.ARCHIVE_DIR, production["product"]), message)


def check_precise_queue() -> None:
    """Work through the queue of urgent-response products (one precise-runner run)."""
    pending = [path for path, entry in load_queue() if entry["status"] == "pending"]
    print(f"{len(pending)} urgent-response products await a production product.")
    for path in pending:
        with open(path) as f:
            entry = json.load(f)
        print(f"  {os.path.basename(entry['urgent_product'])}")
        others = [e for p, e in load_queue() if p != path]
        before = json.dumps(entry, sort_keys=True)
        try:
            process_entry(entry, others)
        except Exception as e:  # one bad entry must not stop the others
            print(f"    Unexpected error (trying again next run): {type(e).__name__}: {e}")
            continue
        if json.dumps(entry, sort_keys=True) == before:
            continue

        def merge(current: dict[str, Any], mine: dict[str, Any] = entry) -> None:
            """Save this run's changes, keeping events the urgent runner added meanwhile."""
            known = {e["event_id"] for e in mine["events"]}
            added = [e for e in current["events"] if e["event_id"] not in known]
            urgent_product = current["urgent_product"]  # the urgent runner may have moved it
            current.update(mine)
            current["urgent_product"] = urgent_product
            current["events"] = mine["events"] + added
            if added and current["status"] == "not_eligible":
                current["status"] = "pending"  # an event joined: it may be significant

        update_entry(path, merge)


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
            <p>Production products made with precise orbits after urgent-response products:</p>
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
