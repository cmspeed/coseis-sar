"""
A product version made from another version's per-path products with another mosaic:
`python -m aria_coseis.optical.remosaic`.

    python -m aria_coseis.optical.remosaic --source exp-g90 --version exp-g90-aligned [--event E]

autoRIFT is not rerun and the source version is not changed. For each event the source made
(status COMPLETE), products/<version>/<event>/ gets:

    paths/                hard links to the source's per-path files (their names keep the
                          source version)
    E_<version>_offsets.tif  the deliverable, mosaicked with the new method
    run_info.json         the source run's inputs and per-path records, the per-path files'
                          checksums, the source run_info's checksum and the alignment

VERSION.json is the source definition with the new deliverable, "derived_from" and
"mosaic_code" (runs.derived_definition). Needs GDAL and scipy, not hyp3_autorift.
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
    """Remosaic one event of the source version; returns its status. Errors are recorded."""
    version = version_record["version"]
    source = version_record["definition"]["derived_from"]["version"]
    source_dir = root / "products" / source / event
    source_run_path = source_dir / "run_info.json"
    if not source_run_path.exists():
        print(f"Skipping {event}: not in {source}")
        return "NOT_IN_SOURCE"
    source_run = json.loads(source_run_path.read_text())
    if source_run.get("status") != "COMPLETE":
        print(f"Skipping {event}: {source} is {source_run.get('status')}")
        return "SOURCE_NOT_COMPLETE"
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
    event_meta = events.get(event) or source_run.get("event_metadata", {})
    run: dict[str, Any] = {
        "event": event,
        "version": version,
        "status": "PROCESSING",
        "started_utc": runs.utc_now(),
        "derived_from": {
            "version": source,
            "run_info_sha256": store.sha256(source_run_path),
            "finished_utc": source_run.get("finished_utc"),
        },
        "event_metadata": event_meta,
        "inputs": source_run.get("inputs", []),
        "skipped": source_run.get("skipped", {}),
        "paths": [],
        "code": provenance,
        "source_code": source_run.get("code"),
    }
    try:
        for record in source_run.get("paths", []):
            record = dict(record)
            if record.get("offset_m"):
                link(source_dir / record["offset_m"], out_dir / record["offset_m"])
                record["offset_m_sha256"] = store.sha256(out_dir / record["offset_m"])
            run["paths"].append(record)
        written = runs.write_event_deliverable(version_record, event, event_meta, run, out_dir)
        run["status"] = "COMPLETE" if written else "NO_DATA"
    except Exception as e:  # noqa: BLE001 - record the failure and go on with the next event
        run["status"] = "FAILED"
        run["error"] = f"{type(e).__name__}: {e}"
        run["traceback"] = traceback.format_exc()
        print(f"FAILED {event}: {run['error']}")
    run["finished_utc"] = runs.utc_now()
    run_path.write_text(json.dumps(run, indent=1) + "\n")
    mosaic_commit = version_record["definition"]["mosaic_code"]["aria_coseis"].get("commit", "")
    store.update_catalog(
        root / "products" / version / "catalog.csv",
        runs.catalog_row(event, event_meta, run, mosaic_commit or ""),
    )
    alignment = run.get("deliverable_grid", {}).get("alignment")
    if alignment:
        offsets = {c: a["offsets_removed_m"] for c, a in alignment.items()}
        print(f"{event}: {run['status']}, path offsets removed (m): {offsets}")
    else:
        print(f"{event}: {run['status']}")
    return run["status"]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="aria-coseis-remosaic",
        description="Make a product version from another version's per-path products with "
        "another mosaic (no autoRIFT rerun; the source version is not changed).",
    )
    parser.add_argument("--source", required=True, help="Version whose per-path files are used")
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
    parser.add_argument("--force", action="store_true", help="Redo events already remosaicked")
    args = parser.parse_args(argv)

    root = store.require_store(args.store).resolve()
    if args.version == args.source:
        parser.error("--version must differ from --source")
    source_path = root / "products" / args.source / "VERSION.json"
    if not source_path.exists():
        parser.error(f"no version {args.source} in {root}")
    provenance = store.code_provenance()
    definition = runs.derived_definition(
        json.loads(source_path.read_text()), args.mosaic, provenance
    )
    description = args.description or (
        f"{args.source} per-path products, mosaicked with {args.mosaic}"
    )
    record = store.start_version(root, args.version, description, definition)
    env_listing = root / "products" / args.version / "environment.txt"
    if not env_listing.exists():
        store.write_environment_listing(env_listing)

    events = args.event or sorted(
        p.parent.name for p in (root / "products" / args.source).glob("*/run_info.json")
    )
    started = time.time()
    statuses = {e: remosaic_event(root, record, e, provenance, args.force) for e in events}
    summary = {s: list(statuses.values()).count(s) for s in sorted(set(statuses.values()))}
    print(f"\nDone in {time.time() - started:.0f} s: {summary}")


if __name__ == "__main__":
    main()
