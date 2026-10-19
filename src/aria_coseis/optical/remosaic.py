"""
A product version made from other versions' per-path products with another mosaic:
`python -m aria_coseis.optical.remosaic`.

    python -m aria_coseis.optical.remosaic --source exp-g90 --version exp-g90-aligned [--event E]
    python -m aria_coseis.optical.remosaic --source exp-g90-med30snow --source exp-g90-seas90snow \
        --version exp-g90-ens --far_field_km 25

autoRIFT is not rerun and the sources are not changed. With several sources (an ensemble), every
(source, path) is a member of the aligned_weighted mosaic: the same path from two composite
recipes is aligned and weighted like two overlapping paths (runs.ensemble_definition). With
--far_field_km, the median beyond that distance from the ruptures is subtracted afterwards
(aria_coseis.optical.farfield) and the ruptures are saved next to the deliverable.

For each event every source made (status COMPLETE), products/<version>/<event>/ gets:

    paths/                hard links to the source's per-path files (their names keep the
                          source version)
    E_<version>_offsets.tif  the deliverable, mosaicked with the new method
    run_info.json         the source run's inputs and per-path records, the per-path files'
                          checksums, the source run_info's checksum and the alignment

VERSION.json is the source definition with the new deliverable, "derived_from" and
"mosaic_code" (runs.derived_definition), or for an ensemble the members' definitions
("ensemble_of"). Needs GDAL and scipy, not hyp3_autorift; the far-field step needs the network
(USGS finite-fault models).
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import time
import traceback
from pathlib import Path
from typing import Any

from aria_coseis.optical import runs, store


def link(src: Path, dst: Path) -> None:
    """Hard-link src to dst (copy if they are on different file systems)."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def remosaic_event(
    root: Path,
    version_record: dict[str, Any],
    event: str,
    provenance: dict[str, Any],
    force: bool = False,
) -> str:
    """
    Remosaic one event from the source version(s); returns its status. Errors are recorded. With
    several sources (an ensemble) every (source, path) is a member of the mosaic, labelled
    <recipe>/P<path>.
    """
    version = version_record["version"]
    definition = version_record["definition"]
    sources = runs.source_versions(definition)
    ensemble = "ensemble_of" in definition
    source_runs = {}
    for source in sources:
        source_run_path = root / "products" / source / event / "run_info.json"
        if not source_run_path.exists():
            print(f"Skipping {event}: not in {source}")
            return "NOT_IN_SOURCE"
        source_run = json.loads(source_run_path.read_text())
        if source_run.get("status") != "COMPLETE":
            print(f"Skipping {event}: {source} is {source_run.get('status')}")
            return "SOURCE_NOT_COMPLETE"
        source_runs[source] = (source_run, source_run_path)
    out_dir = root / "products" / version / event
    run_path = out_dir / "run_info.json"
    if run_path.exists() and not force:
        status = json.loads(run_path.read_text()).get("status")
        if status == "COMPLETE":
            print(f"Skipping {event}: {version} is {status} (use --force to redo)")
            return status
    if out_dir.exists():
        shutil.rmtree(out_dir)  # a rerun replaces the event's products entirely
    out_dir.mkdir(parents=True)

    events_path = root / "events.json"
    events = json.loads(events_path.read_text()) if events_path.exists() else {}
    first_run = source_runs[sources[0]][0]
    event_meta = events.get(event) or first_run.get("event_metadata", {})
    origin = {
        s: {"run_info_sha256": store.sha256(path), "finished_utc": r.get("finished_utc")}
        for s, (r, path) in source_runs.items()
    }
    run: dict[str, Any] = {
        "event": event,
        "version": version,
        "status": "PROCESSING",
        "started_utc": runs.utc_now(),
        "derived_from": {"versions": origin}
        if ensemble
        else {"version": sources[0], **origin[sources[0]]},
        "event_metadata": event_meta,
        "inputs": [],
        "skipped": {},
        "paths": [],
        "code": provenance,
        "source_code": {s: r.get("code") for s, (r, _) in source_runs.items()}
        if ensemble
        else first_run.get("code"),
    }
    try:
        for source, (source_run, _) in source_runs.items():
            recipe = (definition["ensemble_of"][source] if ensemble else definition)[
                "composite_recipe"
            ]
            tag = {"version": source, "recipe": recipe} if ensemble else {}
            run["inputs"] += [{**i, **tag} for i in source_run.get("inputs", [])]
            for p_, reason in source_run.get("skipped", {}).items():
                run["skipped"][p_] = f"{recipe}: {reason}" if ensemble else reason
            for record in source_run.get("paths", []):
                record = {**record, **tag}
                if ensemble:
                    record["member"] = f"{recipe}/P{record['path']}"
                if record.get("offset_m"):
                    src = root / "products" / source / event / record["offset_m"]
                    link(src, out_dir / record["offset_m"])
                    record["offset_m_sha256"] = store.sha256(out_dir / record["offset_m"])
                run["paths"].append(record)
        written = runs.write_event_deliverable(
            version_record, event, event_meta, run, out_dir, catalog=events
        )
        run["status"] = "COMPLETE" if written else "NO_DATA"
    except Exception as e:  # noqa: BLE001 - record the failure and go on with the next event
        run["status"] = "FAILED"
        run["error"] = f"{type(e).__name__}: {e}"
        run["traceback"] = traceback.format_exc()
        print(f"FAILED {event}: {run['error']}")
    run["finished_utc"] = runs.utc_now()
    run_path.write_text(json.dumps(run, indent=1) + "\n")
    mosaic_commit = definition["mosaic_code"]["aria_coseis"].get("commit", "")
    store.update_catalog(
        root / "products" / version / "catalog.csv",
        runs.catalog_row(event, event_meta, run, mosaic_commit or ""),
    )
    alignment = run.get("deliverable_grid", {}).get("alignment")
    shift = run.get("far_field", {}).get("shift_m")
    note = f", far-field shift removed {shift} m" if shift else ""
    if alignment:
        offsets = {c: a["offsets_removed_m"] for c, a in alignment.items()}
        print(f"{event}: {run['status']}, offsets removed (m): {offsets}{note}")
    else:
        print(f"{event}: {run['status']}{note}")
    return run["status"]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="aria-coseis-remosaic",
        description="Make a product version from another version's per-path products with "
        "another mosaic (no autoRIFT rerun; the source version is not changed).",
    )
    parser.add_argument(
        "--source",
        action="append",
        required=True,
        help="Version whose per-path files are used (repeat for an ensemble of versions)",
    )
    parser.add_argument("--version", required=True, help="New version, e.g. exp-g90-aligned")
    parser.add_argument("--description", default="", help="What the version is for (first run)")
    parser.add_argument("--store", help="Store root (default: COSEIS_OPTICAL_DIR)")
    parser.add_argument(
        "--event", action="append", help="Event folder name (repeatable); default: all events"
    )
    parser.add_argument(
        "--mosaic",
        choices=sorted(runs.DELIVERABLE_BANDS),
        default=runs.DEFAULT_MOSAIC,
        help=f"How paths are combined (default {runs.DEFAULT_MOSAIC})",
    )
    parser.add_argument(
        "--far_field_km",
        type=float,
        default=0.0,
        help="Subtract the median beyond this distance from the ruptures (default 0: none)",
    )
    parser.add_argument("--force", action="store_true", help="Redo events already remosaicked")
    args = parser.parse_args(argv)

    root = store.require_store(args.store).resolve()
    if args.version in args.source:
        parser.error("--version must differ from --source")
    records = []
    for source in args.source:
        source_path = root / "products" / source / "VERSION.json"
        if not source_path.exists():
            parser.error(f"no version {source} in {root}")
        records.append(json.loads(source_path.read_text()))
    provenance = store.code_provenance()
    if len(records) > 1:
        definition = runs.ensemble_definition(records, args.mosaic, provenance, args.far_field_km)
    else:
        definition = runs.derived_definition(
            records[0], args.mosaic, provenance, args.far_field_km
        )
    description = args.description or (
        f"{' + '.join(args.source)} per-path products, mosaicked with {args.mosaic}"
        + (
            f", far-field median beyond {args.far_field_km:g} km removed"
            if args.far_field_km
            else ""
        )
    )
    record = store.start_version(root, args.version, description, definition)
    env_listing = root / "products" / args.version / "environment.txt"
    if not env_listing.exists():
        store.write_environment_listing(env_listing)

    events = args.event or sorted(
        p.parent.name for p in (root / "products" / args.source[0]).glob("*/run_info.json")
    )
    started = time.time()
    statuses = {e: remosaic_event(root, record, e, provenance, args.force) for e in events}
    summary = {s: list(statuses.values()).count(s) for s in sorted(set(statuses.values()))}
    print(f"\nDone in {time.time() - started:.0f} s: {summary}")


if __name__ == "__main__":
    main()
