"""Store-based optical runs: input pairs, version definitions, catalog rows and the deliverable."""

import json
from pathlib import Path
from typing import Any, Dict

import pytest

import coseis_api as api

runs = api.optical_runs
store = api.optical_store

EVENT = "m_71_ridgecrest_earthquake_sequence"
CODE = {
    "aria_coseis": {
        "installed": True,
        "location": "/a",
        "repo": "/a",
        "commit": "a" * 40,
        "dirty": False,
    },
    "hyp3_autorift": {
        "installed": True,
        "location": "/h",
        "repo": "/h",
        "commit": "b" * 40,
        "dirty": False,
    },
    "autoRIFT": {
        "installed": True,
        "location": "/c",
        "repo": None,
        "distribution_version": "2.0.0",
    },
}


def composite_index(root: Path, entries: Dict[str, Dict[str, Any]]) -> None:
    folder = root / "composites" / EVENT
    folder.mkdir(parents=True)
    (folder / "composites.json").write_text(json.dumps(entries))


def entry(path: int, stage: str, recipe: str = "med90") -> Dict[str, Any]:
    window = ["2019-04-07T00:00:00Z", "2019-07-06T03:19:53Z"]
    if stage == "POST":
        window = ["2019-07-06T03:19:53Z", "2019-10-04T00:00:00Z"]
    return {
        "path": path,
        "stage": stage,
        "recipe": recipe,
        "start": window[0],
        "end": window[1],
        "sha256": f"{stage}{path}",
    }


def test_composite_pairs(tmp_path: Path) -> None:
    composite_index(
        tmp_path,
        {
            "S2_P27_PRE": entry(27, "PRE"),
            "S2_P27_POST": entry(27, "POST"),
            "S2_P70_PRE": entry(70, "PRE"),
            "S2_P70_PRE_again": entry(70, "PRE"),
            "S2_P70_POST": entry(70, "POST"),
            "S2_P27_PRE_anniv": entry(27, "PRE", recipe="anniv90"),
            "_scene_times_all_paths": {"pre_dates_used": []},
        },
    )
    pairs, skipped = runs.composite_pairs(tmp_path, EVENT, "med90")
    assert [p["path"] for p in pairs] == [27]
    assert pairs[0]["ref_name"] == "S2_P27_PRE" and pairs[0]["sec_sha256"] == "POST27"
    assert pairs[0]["ref"] == str(tmp_path / "composites" / EVENT / "S2_P27_PRE")
    assert skipped == {"70": "2 PRE and 1 POST composites"}
    assert runs.composite_pairs(tmp_path, "no_such_event", "med90") == ([], {})


def definition(profile_name: str = "g90", code: Dict[str, Any] = CODE) -> Dict[str, Any]:
    profile = api.AUTORIFT_PROFILES[profile_name]
    geom = api.pixel_geometry(profile, 10.0)
    return runs.version_definition(profile, "med90", geom, api.driver_options(profile, geom), code)


def test_version_definition_identifies_code_not_its_location() -> None:
    d = definition()
    assert d["code"]["aria_coseis"] == {"commit": "a" * 40, "dirty": False}
    assert d["code"]["autoRIFT"] == {
        "installed": True,
        "distribution_version": "2.0.0",
        "conda_package": None,
    }
    assert d["driver_options"]["NullReferenceVelocity"] is True
    moved = {**CODE, "aria_coseis": {**CODE["aria_coseis"], "location": "/elsewhere"}}
    assert definition(code=moved) == d


def test_a_version_keeps_its_definition(tmp_path: Path) -> None:
    record = store.start_version(tmp_path, "v1-90m", "first", definition())
    assert store.start_version(tmp_path, "v1-90m", "again", definition()) == record
    with pytest.raises(ValueError, match="different driver_options, profile"):
        store.start_version(tmp_path, "v1-90m", "", definition("legacy90"))
    other_code = {**CODE, "hyp3_autorift": {**CODE["hyp3_autorift"], "commit": "c" * 40}}
    with pytest.raises(ValueError, match="different code"):
        store.start_version(tmp_path, "v1-90m", "", definition(code=other_code))


def test_only_experiments_may_use_uncommitted_code(tmp_path: Path) -> None:
    dirty = {**CODE, "aria_coseis": {**CODE["aria_coseis"], "dirty": True, "diff_sha256": "d"}}
    with pytest.raises(ValueError, match="committed code"):
        store.start_version(tmp_path, "v1-90m", "", definition(code=dirty))
    assert store.start_version(tmp_path, "exp-g90", "", definition(code=dirty))
    for bad in ("v0", "../x", ""):
        with pytest.raises(ValueError, match="invalid version"):
            store.start_version(tmp_path, bad, "", definition())


def test_catalog_row() -> None:
    run = {
        "status": "COMPLETE",
        "finished_utc": "2026-10-15T00:00:00Z",
        "skipped": {"70": "2 PRE and 1 POST composites"},
        "paths": [
            {"path": 27, "offset_m": "paths/x.tif"},
            {"path": 33, "skipped": "empty composite"},
        ],
        "deliverable": f"{EVENT}_exp-g90_offsets.tif",
    }
    row = runs.catalog_row(EVENT, {"event_id": "ci38457511", "magnitude": 7.1}, run, "a" * 40)
    assert row["paths_used"] == "27"
    assert row["paths_skipped"] == "33: empty composite; 70: 2 PRE and 1 POST composites"
    assert row["aria_code"] == "aaaaaaa"
    assert row["mosaic"] == f"{EVENT}/{EVENT}_exp-g90_offsets.tif"


def write_path_product(path: Path, x0: float, dx: float, dy: float, nan_col: int) -> str:
    np = pytest.importorskip("numpy")
    gdal = pytest.importorskip("osgeo.gdal")
    osr = pytest.importorskip("osgeo.osr")
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32611)
    ds = gdal.GetDriverByName("GTiff").Create(str(path), 4, 3, 4, gdal.GDT_Float32)
    ds.SetGeoTransform((x0, 90.0, 0.0, 3_960_000.0, 0.0, -90.0))
    ds.SetProjection(srs.ExportToWkt())
    for band, value in ((1, dx), (2, dy), (3, 0.0), (4, 36.0)):
        arr = np.full((3, 4), value, np.float32)
        arr[:, nan_col] = np.nan
        ds.GetRasterBand(band).WriteArray(arr)
    ds = None
    return str(path)


def test_deliverable_is_the_mean_of_the_paths(tmp_path: Path) -> None:
    np = pytest.importorskip("numpy")
    gdal = pytest.importorskip("osgeo.gdal")
    a = write_path_product(tmp_path / "a.tif", 400_000.0, 1.0, 2.0, nan_col=0)
    b = write_path_product(tmp_path / "b.tif", 400_180.0, 3.0, -2.0, nan_col=3)  # 2 nodes east
    out = tmp_path / "out.tif"
    stats = runs.write_deliverable([a, b], str(out), {"COSEIS_VERSION": "exp-test"})
    assert stats == {"nodes": [3, 6], "valid_nodes": 12, "paths": 2}
    ds = gdal.Open(str(out))
    assert ds.GetGeoTransform()[0] == 400_000.0 and ds.RasterCount == 2
    assert [ds.GetRasterBand(i).GetDescription() for i in (1, 2)] == ["dx_m", "dy_m"]
    assert np.isnan(ds.GetRasterBand(1).GetNoDataValue())
    assert ds.GetMetadataItem("COSEIS_VERSION") == "exp-test"
    dx = ds.GetRasterBand(1).ReadAsArray()[0]
    # columns: a only (a's col 0 is NaN), a only, both, both, b only, b's NaN col
    np.testing.assert_allclose(dx, [np.nan, 1.0, 2.0, 2.0, 3.0, np.nan])
    off = write_path_product(tmp_path / "c.tif", 400_045.0, 0.0, 0.0, nan_col=0)
    with pytest.raises(ValueError, match="offset from the grid"):
        runs.write_deliverable([a, off], str(tmp_path / "bad.tif"), {})
