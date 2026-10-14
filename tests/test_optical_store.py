"""The optical data store: composites by recipe, product versions, and the v0 import."""

import csv
import json
import os
import stat
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import pytest

import coseis_api as api

store = api.optical_store

EVENT = "m_78_pazarcik_earthquake_kahramanmaras_earthquake_sequence"
FAILED = "m_61_96_km_ssw_of_altata_mexico"
ORIGIN = datetime(2023, 2, 6, 1, 17, 34, tzinfo=timezone.utc)


def write(path: Path, content: bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


@pytest.fixture
def downloads(tmp_path: Path) -> Path:
    """A miniature data/GEE_Optical_Downloads: one finished event (two paths, one skipped), one failed."""
    root = tmp_path / "GEE_Optical_Downloads"
    ev = root / EVENT
    for path in ("121", "78"):
        write(
            ev
            / f"{EVENT}_S2_Path{path}_TOA_B8_PRE_2022-11-08_000000UTC_to_2023-02-06_011734UTC.tif",
            f"pre{path}".encode(),
        )
        write(
            ev
            / f"{EVENT}_S2_Path{path}_TOA_B8_POST_2023-02-06_011734UTC_to_2023-05-07_000000UTC.tif",
            f"post{path}".encode(),
        )
        write(ev / f"{EVENT}_Path{path}_sentinel-2_toa_autorift_offset_m.tif", b"offset")
        write(ev / f"{EVENT}_Path{path}_sentinel-2_toa_autorift_velocity.tif", b"velocity")
    write(ev / f"{EVENT}_final_mosaic_sentinel-2_toa_offset_m.tif", b"mosaic")
    manifest = {
        "event_title": EVENT,
        "status": "PROCESSED_AUTORIFT_COMPLETE",
        "track_pairs": {"121": {}, "78": {}},
        "skipped_track_pairs": {"21": {"reason": "empty composite"}},
        "pre_dates_used": ["2022-11-09T08:39:26"],
        "post_dates_used": ["2023-02-06T08:19:26"],
    }
    write(ev / f"{EVENT}_sentinel-2_toa_autorift_manifest.json", json.dumps(manifest).encode())
    failed = {
        "event_title": FAILED,
        "status": "FAILED_NO_DATA",
        "track_pairs": {},
        "event_id": "usb000sjt0",
    }
    write(
        root / FAILED / f"{FAILED}_sentinel-2_toa_autorift_manifest.json",
        json.dumps(failed).encode(),
    )
    return root


def fake_query(calls: List[Any]) -> Any:
    def query(origin: datetime, magnitude: float) -> Optional[Dict[str, Any]]:
        calls.append((origin, magnitude))
        return {
            "event_id": "us6000jllz",
            "title": "M 7.8 - Pazarcik",
            "time_utc": "2023-02-06T01:17:34.86Z",
        }

    return query


def test_composite_name() -> None:
    name = store.composite_name(
        "S2", 121, "TOA", "B8", "PRE", "med90", datetime(2022, 11, 8), ORIGIN
    )
    assert name == "S2_P121_TOA_B8_PRE_med90_20221108T000000Z_20230206T011734Z.tif"
    with pytest.raises(ValueError):
        store.composite_name("S2", 121, "TOA", "B8", "PRE", "nope", ORIGIN, ORIGIN)


def test_import_v0_registers_composites_products_and_catalog(
    downloads: Path, tmp_path: Path
) -> None:
    root = tmp_path / "store"
    calls: List[Any] = []
    by_id: List[str] = []
    rows = store.import_v0(
        downloads,
        root,
        query=fake_query(calls),
        query_by_id=lambda i: by_id.append(i) or {"event_id": i},
    )

    # Composites: hard links under recipe names, read-only, with a record each
    comp = root / "composites" / EVENT
    pre = comp / "S2_P121_TOA_B8_PRE_med90_20221108T000000Z_20230206T011734Z.tif"
    source = next(downloads.glob(f"{EVENT}/*Path121*_PRE_*.tif"))
    assert os.path.samefile(pre, source)
    assert not os.stat(pre).st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
    index = json.loads((comp / "composites.json").read_text())
    assert len([k for k in index if not k.startswith("_")]) == 4
    assert index[pre.name]["recipe"] == "med90"
    assert index[pre.name]["end"] == "2023-02-06T01:17:34Z"
    assert index[pre.name]["sha256"] == store.sha256(source)
    assert index["_scene_times_all_paths"]["pre_dates_used"] == ["2022-11-09T08:39:26"]
    assert json.loads((root / "composites" / "recipes.json").read_text())["med90"]["band"] == "B8"

    # Products: mosaic and per-path offsets linked under their original names; velocity left out
    v0 = root / "products" / "v0"
    assert sorted(p.name for p in (v0 / EVENT).iterdir()) == sorted(
        [f"{EVENT}_final_mosaic_sentinel-2_toa_offset_m.tif"]
        + [f"{EVENT}_Path{p}_sentinel-2_toa_autorift_offset_m.tif" for p in ("121", "78")]
    )
    assert json.loads((v0 / "VERSION.json").read_text())["composite_recipe"] == "med90"

    # Catalog: one row per event, with the code that made it and the USGS id
    with open(v0 / "catalog.csv") as f:
        catalog = {r["event"]: r for r in csv.DictReader(f)}
    assert catalog[EVENT]["paths_used"] == "78 121"
    assert catalog[EVENT]["paths_skipped"] == "21: empty composite"
    assert catalog[EVENT]["aria_code"] == "develop-8b7b7f7"
    assert catalog[EVENT]["event_id"] == "us6000jllz"
    assert catalog[FAILED]["status"] == "FAILED_NO_DATA" and catalog[FAILED]["mosaic"] == ""
    assert calls == [(ORIGIN, 7.8)]  # origin from the composite names
    assert by_id == ["usb000sjt0"] and catalog[FAILED]["event_id"] == "usb000sjt0"  # manifest id
    assert len(rows) == 2

    # Again: nothing changes and USGS is not asked twice
    store.import_v0(downloads, root, query=fake_query(calls), query_by_id=by_id.append)
    assert len(calls) == 1 and len(by_id) == 1
    assert json.loads((comp / "composites.json").read_text())[pre.name] == index[pre.name]


def test_a_different_file_under_an_existing_name_is_refused(
    downloads: Path, tmp_path: Path
) -> None:
    root = tmp_path / "store"
    store.import_v0(downloads, root, resolve=False)
    other = write(tmp_path / "other.tif", b"something else")
    record = {
        "sensor": "S2",
        "path": 121,
        "level": "TOA",
        "band": "B8",
        "stage": "PRE",
        "recipe": "med90",
        "start": datetime(2022, 11, 8),
        "end": ORIGIN,
    }
    with pytest.raises(FileExistsError):
        store.add_composite(root, EVENT, other, record)


def test_versions_and_recipes_cannot_be_redefined(tmp_path: Path) -> None:
    store.write_version(tmp_path, {"version": "exp-a", "description": "first"})
    store.write_version(tmp_path, {"version": "exp-a", "description": "first"})  # same: fine
    with pytest.raises(FileExistsError):
        store.write_version(tmp_path, {"version": "exp-a", "description": "changed"})
    (tmp_path / "composites").mkdir()
    (tmp_path / "composites" / "recipes.json").write_text(json.dumps({"med90": {"band": "B4"}}))
    with pytest.raises(ValueError):
        store.write_recipes(tmp_path)


def test_code_provenance_records_the_aria_coseis_checkout() -> None:
    prov = store.code_provenance()
    aria = prov["aria_coseis"]
    assert aria["installed"]
    if aria.get("repo"):  # a git checkout (always, in development and CI)
        assert len(aria["commit"]) == 40 and isinstance(aria["dirty"], bool)
