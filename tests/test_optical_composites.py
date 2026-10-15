"""Composites of new recipes for store events: windows, grid, plan, reuse and cropping."""

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict

import pytest

import coseis_api as api

composites = api.optical_composites
store = api.optical_store

EVENT = "m_71_ridgecrest_earthquake_sequence"
ORIGIN = datetime(2019, 7, 6, 3, 19, 53, tzinfo=timezone.utc)
MED90 = {
    ("PRE", 27): "S2_P27_TOA_B8_PRE_med90_20190407T000000Z_20190706T031953Z.tif",
    ("POST", 27): "S2_P27_TOA_B8_POST_med90_20190706T031953Z_20191004T000000Z.tif",
}


def fake_store(root: Path) -> Path:
    """A store with Ridgecrest's med90 P27 pair (small dummy files)."""
    folder = root / "composites" / EVENT
    folder.mkdir(parents=True)
    index: Dict[str, Any] = {}
    for (stage, path), name in MED90.items():
        (folder / name).write_bytes(stage.encode())
        start, end = name[:-4].split("_")[-2:]
        index[name] = {
            "path": path,
            "stage": stage,
            "recipe": "med90",
            "start": datetime.strptime(start, "%Y%m%dT%H%M%SZ").strftime("%Y-%m-%dT%H:%M:%SZ"),
            "end": datetime.strptime(end, "%Y%m%dT%H%M%SZ").strftime("%Y-%m-%dT%H:%M:%SZ"),
            "sha256": store.sha256(folder / name),
        }
    (folder / "composites.json").write_text(json.dumps(index))
    return folder


def test_windows_reproduce_med90_and_shift_a_year() -> None:
    # med90snow has med90's windows, so its names differ from med90's only in the recipe
    for stage in ("PRE", "POST"):
        start, end = composites.stage_window("med90snow", stage, ORIGIN)
        name = store.composite_name("S2", 27, "TOA", "B8", stage, "med90snow", start, end)
        assert name == MED90[(stage, 27)].replace("med90", "med90snow")
    # seas90 PRE is the med90 POST window one year (365 d) earlier
    start, end = composites.stage_window("seas90", "PRE", ORIGIN)
    assert (start, end) == (
        datetime(2018, 7, 6, 3, 19, 53, tzinfo=timezone.utc),
        datetime(2018, 10, 4, tzinfo=timezone.utc),
    )
    with pytest.raises(ValueError, match="unknown window bound"):
        composites.window_bound(ORIGIN, ["week", 1])


def test_snapped_grid_matches_the_med90_composites() -> None:
    # Earth Engine's UTM bounds of the AOIs and the grids of the med90 composites (2026-10-15)
    ridgecrest = composites.snapped_grid(375549.4, 3897152.7, 497384.4, 4038269.5)
    assert ridgecrest == {"x0": 375540, "y0": 4038270, "cols": 12185, "rows": 14112}
    pazarcik = composites.snapped_grid(220404.8, 3971069.0, 487989.5, 4251921.4)
    assert pazarcik == {"x0": 220400, "y0": 4251930, "cols": 26759, "rows": 28087}
    assert composites.utm_epsg(-117.6, 35.77) == 32611
    assert composites.utm_epsg(119.85, -0.26) == 32750


def test_new_recipes_are_complete_and_med90_is_unchanged() -> None:
    recipes = store.COMPOSITE_RECIPES
    assert recipes["med90"]["pre_window"] == "[event - 90 d, event)"
    assert recipes["med90"]["snow_mask"] is None and "window_days" not in recipes["med90"]
    made_here = [r for r, d in recipes.items() if "window_days" in d]
    assert sorted(made_here) == ["med90snow", "seas90", "seas90snow"]
    for recipe in made_here:
        d = recipes[recipe]
        reused = d.get("stage_from_recipe", {})
        assert set(d["window_days"]) | set(reused) == {"PRE", "POST"}
        assert all(source in recipes for source in reused.values())
        assert d["cloud_score_min"] == 0.65 and d["scale_m"] == 10
    assert recipes["seas90snow"]["snow_mask"] == recipes["med90snow"]["snow_mask"]
    assert recipes["seas90"]["window_days"] == recipes["seas90snow"]["window_days"]


def test_plan_exports_new_windows_and_reuses_stages(tmp_path: Path) -> None:
    fake_store(tmp_path)
    plan = composites.planned_composites(tmp_path, EVENT, "seas90", ORIGIN)
    assert [(i["stage"], i["action"]) for i in plan] == [("PRE", "export"), ("POST", "reuse")]
    assert plan[0]["name"] == "S2_P27_TOA_B8_PRE_seas90_20180706T031953Z_20181004T000000Z.tif"
    assert plan[1]["source"] == MED90[("POST", 27)]
    assert plan[1]["name"] == MED90[("POST", 27)].replace("med90", "seas90")
    with pytest.raises(ValueError, match="make med90snow first"):
        composites.planned_composites(tmp_path, EVENT, "seas90snow", ORIGIN)
    later = composites.planned_composites(tmp_path, EVENT, "seas90snow", ORIGIN, strict=False)
    assert later[1]["action"] == "reuse later"
    assert later[1]["source"] == MED90[("POST", 27)].replace("med90", "med90snow")
    with pytest.raises(ValueError, match="not made by this module"):
        composites.planned_composites(tmp_path, EVENT, "med90", ORIGIN)


def test_a_reused_stage_is_the_same_file(tmp_path: Path) -> None:
    folder = fake_store(tmp_path)
    plan = composites.planned_composites(tmp_path, EVENT, "seas90", ORIGIN)
    code = {"commit": "a" * 40, "dirty": False}
    dst = composites.register_reused(tmp_path, EVENT, plan[1], "seas90", code)
    assert dst.stat().st_ino == (folder / MED90[("POST", 27)]).stat().st_ino
    index = json.loads((folder / "composites.json").read_text())
    assert index[dst.name]["recipe"] == "seas90"
    assert index[dst.name]["same_file_as"] == MED90[("POST", 27)]
    assert index[dst.name]["registered_by"] == code
    again = composites.planned_composites(tmp_path, EVENT, "seas90", ORIGIN)
    assert again[1]["action"] == "registered"


def test_pipeline_downloads_cannot_be_registered_as_store_recipes(tmp_path: Path) -> None:
    (tmp_path / "README.md").write_text("store\n")
    with pytest.raises(ValueError, match="made by aria_coseis.optical.composites"):
        store.register_downloads(tmp_path, root=tmp_path, recipe="seas90", resolve=False)


def test_crop_to_grid(tmp_path: Path) -> None:
    np = pytest.importorskip("numpy")
    gdal = pytest.importorskip("osgeo.gdal")
    osr = pytest.importorskip("osgeo.osr")
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32611)
    path = tmp_path / "c.tif"
    ds = gdal.GetDriverByName("GTiff").Create(str(path), 12, 10, 1, gdal.GDT_UInt16)
    ds.SetGeoTransform((375530.0, 10.0, 0.0, 4038280.0, 0.0, -10.0))
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(np.arange(120, dtype=np.uint16).reshape(10, 12))
    ds = None
    got = composites.raster_grid(path)
    grid = {"x0": 375540.0, "y0": 4038270.0, "cols": 10, "rows": 8, "epsg": 32611}
    off_lattice = {**grid, "x0": 375545.0}
    assert not composites.crop_to_grid(path, got, off_lattice)
    assert composites.crop_to_grid(path, got, grid)
    assert composites.raster_grid(path) == grid
    cropped = gdal.Open(str(path))  # keep the dataset while reading its band
    assert cropped.GetRasterBand(1).ReadAsArray()[0, 0] == 13  # row 1, column 1 of the original
    cropped = None
    assert os.listdir(tmp_path) == ["c.tif"]
