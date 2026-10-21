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


def test_mean_mosaic_spec_is_unchanged() -> None:
    # exp-g90 and exp-legacy90-repro were defined with this deliverable
    spec = runs.deliverable_spec("mean")
    assert spec["bands"] == ["dx_m", "dy_m"]
    assert spec["mosaic"] == "mean of the paths with a value at each node"
    assert runs.mosaic_method({"deliverable": spec}) == "mean"
    aligned = definition()["deliverable"]
    assert aligned["bands"] == ["dx_m", "dy_m", "sigma_h_m"]
    assert runs.mosaic_method({"deliverable": aligned}) == "aligned_weighted"


def test_derived_definition(tmp_path: Path) -> None:
    source = store.start_version(tmp_path, "exp-g90", "", definition())
    mosaic_code = {**CODE, "aria_coseis": {**CODE["aria_coseis"], "commit": "e" * 40}}
    derived = runs.derived_definition(source, "aligned_weighted", mosaic_code)
    assert derived["code"] == source["definition"]["code"]  # the code that made the paths
    assert derived["profile"] == source["definition"]["profile"]
    assert derived["derived_from"]["version"] == "exp-g90"
    assert derived["mosaic_code"] == {"aria_coseis": {"commit": "e" * 40, "dirty": False}}
    record = store.start_version(tmp_path, "v1-test", "", derived)
    with pytest.raises(ValueError, match="remosaicked itself; remosaic exp-g90"):
        runs.derived_definition(record, "aligned_weighted", CODE)
    with pytest.raises(ValueError, match="not made by a store run"):
        runs.derived_definition({"version": "v0", "code": {}}, "aligned_weighted", CODE)
    dirty = {**CODE, "aria_coseis": {**CODE["aria_coseis"], "dirty": True, "diff_sha256": "d"}}
    with pytest.raises(ValueError, match="committed code"):
        store.start_version(
            tmp_path, "v1-other", "", runs.derived_definition(source, "aligned_weighted", dirty)
        )


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
    stats = runs.write_deliverable({27: a, 70: b}, str(out), {"COSEIS_VERSION": "exp-test"})
    assert stats == {"nodes": [3, 6], "valid_nodes": 12, "paths": 2}
    ds = gdal.Open(str(out))
    assert ds.GetGeoTransform()[0] == 400_000.0 and ds.RasterCount == 2
    assert [ds.GetRasterBand(i).GetDescription() for i in (1, 2)] == ["dx_m", "dy_m"]
    assert np.isnan(ds.GetRasterBand(1).GetNoDataValue())
    assert ds.GetMetadataItem("COSEIS_VERSION") == "exp-test"
    assert ds.GetMetadataItem("MOSAIC_METHOD") == "mean"
    dx = ds.GetRasterBand(1).ReadAsArray()[0]
    # columns: a only (a's col 0 is NaN), a only, both, both, b only, b's NaN col
    np.testing.assert_allclose(dx, [np.nan, 1.0, 2.0, 2.0, 3.0, np.nan])
    off = write_path_product(tmp_path / "c.tif", 400_045.0, 0.0, 0.0, nan_col=0)
    with pytest.raises(ValueError, match="offset from the grid"):
        runs.write_deliverable({27: a, 70: off}, str(tmp_path / "bad.tif"), {})


def write_noisy_path(path: Path, x0: float, shift: float, seed: int) -> str:
    """A 50 x 40-node path product: dx = shift + noise, dy = -shift + noise (0.1 m)."""
    np = pytest.importorskip("numpy")
    gdal = pytest.importorskip("osgeo.gdal")
    osr = pytest.importorskip("osgeo.osr")
    rng = np.random.default_rng(seed)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(32611)
    ds = gdal.GetDriverByName("GTiff").Create(str(path), 40, 50, 4, gdal.GDT_Float32)
    ds.SetGeoTransform((x0, 90.0, 0.0, 3_960_000.0, 0.0, -90.0))
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray((shift + rng.normal(0, 0.1, (50, 40))).astype(np.float32))
    ds.GetRasterBand(2).WriteArray((-shift + rng.normal(0, 0.1, (50, 40))).astype(np.float32))
    ds.GetRasterBand(3).WriteArray(np.zeros((50, 40), np.float32))
    ds.GetRasterBand(4).WriteArray(np.full((50, 40), 4.0, np.float32))
    ds = None
    return str(path)


def test_aligned_deliverable_has_uncertainty_and_path_offsets(tmp_path: Path) -> None:
    np = pytest.importorskip("numpy")
    gdal = pytest.importorskip("osgeo.gdal")
    pytest.importorskip("scipy")
    a = write_noisy_path(tmp_path / "a.tif", 400_000.0, 0.0, seed=1)
    b = write_noisy_path(tmp_path / "b.tif", 400_000.0 + 10 * 90, 0.5, seed=2)  # 1500 shared
    out = tmp_path / "out.tif"
    stats = runs.write_deliverable({27: a, 70: b}, str(out), {}, "aligned_weighted")
    assert stats["valid_nodes"] == 50 * 50
    assert stats["alignment"]["dx"]["reference_path"] == 27
    assert stats["alignment"]["dx"]["overlaps"][0]["paths"] == [27, 70]
    assert stats["alignment"]["dx"]["offsets_removed_m"]["70"] == pytest.approx(0.5, abs=0.02)
    assert stats["alignment"]["dy"]["offsets_removed_m"]["70"] == pytest.approx(-0.5, abs=0.02)
    ds = gdal.Open(str(out))
    names = [ds.GetRasterBand(i).GetDescription() for i in (1, 2, 3)]
    assert names == ["dx_m", "dy_m", "sigma_h_m"] and ds.RasterCount == 3
    assert ds.GetRasterBand(3).GetUnitType() == "m"
    assert ds.GetMetadataItem("MOSAIC_METHOD") == "aligned_weighted"
    assert ds.GetMetadataItem("REFERENCE_PATH") == "27"
    removed = json.loads(ds.GetMetadataItem("PATH_OFFSETS_REMOVED_M"))
    assert set(removed) == {"dx", "dy"} and removed["dx"]["27"] == 0.0
    dx = ds.GetRasterBand(1).ReadAsArray()
    assert abs(np.median(dx[:, 40:]) - np.median(dx[:, :10])) < 0.05  # no step between paths
    assert np.isfinite(ds.GetRasterBand(3).ReadAsArray()).all()


def make_source_version(root: Path) -> Path:
    """A store with version exp-src (mean mosaic) made for EVENT from two path products."""
    (root / "README.md").write_text("store\n")
    (root / "events.json").write_text(json.dumps({EVENT: {"event_id": "ci38457511"}}))
    record = store.start_version(root, "exp-src", "", definition())
    record["definition"]["deliverable"] = runs.deliverable_spec("mean")
    (root / "products" / "exp-src" / "VERSION.json").write_text(json.dumps(record))
    out_dir = root / "products" / "exp-src" / EVENT
    (out_dir / "paths").mkdir(parents=True)
    paths = []
    for n, (path, x0, shift) in enumerate(((27, 400_000.0, 0.0), (70, 400_900.0, 0.5))):
        name = runs.path_product_name(EVENT, "exp-src", path)
        write_noisy_path(out_dir / "paths" / name, x0, shift, seed=n)
        paths.append({"path": path, "offset_m": f"paths/{name}"})
    run = {"event": EVENT, "status": "COMPLETE", "inputs": [], "skipped": {}, "paths": paths}
    assert runs.write_event_deliverable(record, EVENT, {}, run, out_dir)
    (out_dir / "run_info.json").write_text(json.dumps(run))
    return out_dir


def test_remosaic_makes_a_new_version_and_leaves_the_source(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    gdal = pytest.importorskip("osgeo.gdal")
    pytest.importorskip("scipy")
    source_dir = make_source_version(tmp_path)
    before = {p: (p.stat().st_mtime_ns, p.stat().st_size) for p in source_dir.rglob("*")}
    monkeypatch.setattr(store, "code_provenance", lambda: CODE)
    monkeypatch.setattr(store, "write_environment_listing", lambda dst: dst.write_text("env\n"))
    remosaic = api.optical_remosaic
    args = ["--store", str(tmp_path), "--source", "exp-src", "--version", "exp-src-aligned"]
    remosaic.main(args)

    assert {p: (p.stat().st_mtime_ns, p.stat().st_size) for p in source_dir.rglob("*")} == before
    out_dir = tmp_path / "products" / "exp-src-aligned" / EVENT
    run = json.loads((out_dir / "run_info.json").read_text())
    assert run["status"] == "COMPLETE" and run["derived_from"]["version"] == "exp-src"
    for record in run["paths"]:
        linked, original = out_dir / record["offset_m"], source_dir / record["offset_m"]
        assert linked.stat().st_ino == original.stat().st_ino
        assert record["offset_m_sha256"] == store.sha256(original)
    ds = gdal.Open(str(out_dir / runs.deliverable_name(EVENT, "exp-src-aligned")))
    assert ds.RasterCount == 3 and ds.GetMetadataItem("DERIVED_FROM_VERSION") == "exp-src"
    assert ds.GetMetadataItem("USGS_EVENT_ID") == "ci38457511"
    version = json.loads((tmp_path / "products" / "exp-src-aligned" / "VERSION.json").read_text())
    assert version["definition"]["derived_from"]["version"] == "exp-src"
    catalog = (tmp_path / "products" / "exp-src-aligned" / "catalog.csv").read_text()
    assert f"{EVENT},ci38457511" in catalog and "COMPLETE" in catalog

    remosaic.main(args)  # done events are skipped
    assert json.loads((out_dir / "run_info.json").read_text()) == run
    with pytest.raises(SystemExit):
        remosaic.main(["--store", str(tmp_path), "--source", "v0", "--version", "exp-x"])


def make_recipe_version(root: Path, version: str, recipe: str, shift: float, seed: int) -> None:
    """A version of EVENT made from one recipe: paths 27 and 70 (P70 is +0.5 m from P27)."""
    record = store.start_version(root, version, "", {**definition(), "composite_recipe": recipe})
    out_dir = root / "products" / version / EVENT
    (out_dir / "paths").mkdir(parents=True)
    paths, inputs = [], []
    for n, (path, x0, offset) in enumerate(((27, 400_000.0, 0.0), (70, 400_900.0, 0.5))):
        name = runs.path_product_name(EVENT, version, path)
        write_noisy_path(out_dir / "paths" / name, x0, offset + shift, seed=seed + n)
        paths.append({"path": path, "offset_m": f"paths/{name}"})
        inputs.append(
            {
                "path": path,
                "ref_name": f"PRE_{recipe}_{path}",
                "sec_name": f"POST_{recipe}_{path}",
                "ref_sha256": "a",
                "sec_sha256": "b",
                "pre_window": ["2019-04-07T00:00:00Z", "2019-07-06T03:19:53Z"],
                "post_window": ["2019-07-06T03:19:53Z", "2019-10-04T00:00:00Z"],
            }
        )
    run = {"event": EVENT, "status": "COMPLETE", "inputs": inputs, "skipped": {}, "paths": paths}
    (out_dir / "run_info.json").write_text(json.dumps(run))
    assert record["version"] == version


def test_ensemble_definition_rules(tmp_path: Path) -> None:
    a = store.start_version(
        tmp_path, "exp-a", "", {**definition(), "composite_recipe": "med30snow"}
    )
    b = store.start_version(
        tmp_path, "exp-b", "", {**definition(), "composite_recipe": "seas90snow"}
    )
    ens = runs.ensemble_definition([a, b], "aligned_weighted", CODE, 25)
    assert set(ens["ensemble_of"]) == {"exp-a", "exp-b"}
    assert ens["deliverable"]["far_field_correction"]["distance_km"] == 25
    assert runs.source_versions(ens) == ["exp-a", "exp-b"]
    assert "far_field_correction" not in runs.deliverable_spec("aligned_weighted")
    other = store.start_version(tmp_path, "exp-c", "", definition("legacy90"))
    with pytest.raises(ValueError, match="share one profile"):
        runs.ensemble_definition([a, other], "aligned_weighted", CODE, 25)
    with pytest.raises(ValueError, match="an aligned mosaic"):
        runs.ensemble_definition([a, b], "mean", CODE, 0)
    plain = runs.ensemble_definition([a, b], "aligned_mean", CODE, 25)
    assert "equal weights" in plain["deliverable"]["mosaic"]["description"]
    total = runs.ensemble_definition([a, b], "aligned_total", CODE, 25)
    spec = total["deliverable"]
    assert spec["mosaic"]["structured_error"]["block_m"] == 10_000
    assert "structured" in spec["band_descriptions"]["sigma_h_m"]
    with pytest.raises(ValueError, match="far-field distance"):
        runs.ensemble_definition([a, b], "aligned_total", CODE, 0)
    dirty = {**CODE, "aria_coseis": {**CODE["aria_coseis"], "dirty": True, "diff_sha256": "d"}}
    a_dirty = {**a, "definition": {**a["definition"], "code": runs.code_signature(dirty)}}
    with pytest.raises(ValueError, match="committed code"):
        store.start_version(
            tmp_path,
            "v1-ens",
            "",
            runs.ensemble_definition([a_dirty, b], "aligned_weighted", CODE, 25),
        )


def test_ensemble_of_two_recipes_with_the_far_field_removed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    np = pytest.importorskip("numpy")
    gdal = pytest.importorskip("osgeo.gdal")
    pytest.importorskip("scipy")
    pytest.importorskip("shapely")
    pyproj = pytest.importorskip("pyproj")
    from shapely.geometry import box

    lon, lat = pyproj.Transformer.from_crs(32611, 4326, always_xy=True).transform(
        402_250, 3_957_750
    )
    (tmp_path / "README.md").write_text("store\n")
    catalog = {
        EVENT: {
            "event_id": "ci38457511",
            "magnitude": 7.1,
            "time_utc": "2019-07-06T03:19:53Z",
            "longitude": lon,
            "latitude": lat,
        }
    }
    (tmp_path / "events.json").write_text(json.dumps(catalog))
    make_recipe_version(tmp_path, "exp-a", "med30snow", 0.2, seed=1)
    make_recipe_version(tmp_path, "exp-b", "seas90snow", 0.5, seed=5)  # 0.3 m from exp-a
    monkeypatch.setattr(store, "code_provenance", lambda: CODE)
    monkeypatch.setattr(store, "write_environment_listing", lambda dst: dst.write_text("env\n"))
    rupture = box(400_000, 3_955_500, 400_300, 3_960_000)  # along the west edge
    monkeypatch.setattr(
        api.optical_farfield,
        "rupture_geometry",
        lambda quakes, epsg: (rupture, [{"folder": q["folder"], "circle_km": 0} for q in quakes]),
    )
    api.optical_remosaic.main(
        [
            "--store",
            str(tmp_path),
            "--source",
            "exp-a",
            "--source",
            "exp-b",
            "--version",
            "exp-ens",
            "--far_field_km",
            "0.5",
        ]
    )

    out_dir = tmp_path / "products" / "exp-ens" / EVENT
    run = json.loads((out_dir / "run_info.json").read_text())
    assert run["status"] == "COMPLETE", run.get("error")
    members = sorted(p["member"] for p in run["paths"])
    assert members == ["med30snow/P27", "med30snow/P70", "seas90snow/P27", "seas90snow/P70"]
    expected = {
        "med30snow/P27": 0.2,
        "med30snow/P70": 0.7,
        "seas90snow/P27": 0.5,
        "seas90snow/P70": 1.0,
    }
    alignment = run["deliverable_grid"]["alignment"]["dx"]
    ref = expected[alignment["reference_path"]]
    for member, removed in alignment["offsets_removed_m"].items():
        assert removed == pytest.approx(expected[member] - ref, abs=0.03)
    assert run["far_field"]["shift_m"] == pytest.approx([ref, -ref], abs=0.05)
    assert [r["folder"] for r in run["far_field"]["ruptures"]] == [EVENT]
    assert (out_dir / f"{EVENT}_exp-ens_ruptures.geojson").exists()
    ds = gdal.Open(str(out_dir / f"{EVENT}_exp-ens_offsets.tif"))
    assert abs(np.nanmedian(ds.GetRasterBand(1).ReadAsArray())) < 0.05
    assert ds.GetMetadataItem("COMPOSITE_RECIPE") == "med30snow + seas90snow"
    assert ds.GetMetadataItem("DERIVED_FROM_VERSION") == "exp-a exp-b"
    assert json.loads(ds.GetMetadataItem("FAR_FIELD_SHIFT_REMOVED_M")) == pytest.approx(
        [ref, -ref], abs=0.05
    )
    version = json.loads((tmp_path / "products" / "exp-ens" / "VERSION.json").read_text())
    assert set(version["definition"]["ensemble_of"]) == {"exp-a", "exp-b"}


def test_ensemble_with_total_error_weights(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    np = pytest.importorskip("numpy")
    gdal = pytest.importorskip("osgeo.gdal")
    pytest.importorskip("scipy")
    pytest.importorskip("shapely")
    pyproj = pytest.importorskip("pyproj")
    from shapely.geometry import box

    lon, lat = pyproj.Transformer.from_crs(32611, 4326, always_xy=True).transform(
        402_250, 3_957_750
    )
    (tmp_path / "README.md").write_text("store\n")
    catalog = {
        EVENT: {
            "event_id": "ci38457511",
            "magnitude": 7.1,
            "time_utc": "2019-07-06T03:19:53Z",
            "longitude": lon,
            "latitude": lat,
        }
    }
    (tmp_path / "events.json").write_text(json.dumps(catalog))
    make_recipe_version(tmp_path, "exp-a", "med30snow", 0.2, seed=1)
    make_recipe_version(tmp_path, "exp-b", "seas90snow", 0.5, seed=5)
    monkeypatch.setattr(store, "code_provenance", lambda: CODE)
    monkeypatch.setattr(store, "write_environment_listing", lambda dst: dst.write_text("env\n"))
    rupture = box(400_000, 3_955_500, 400_300, 3_960_000)
    monkeypatch.setattr(
        api.optical_farfield,
        "rupture_geometry",
        lambda quakes, epsg: (rupture, [{"folder": q["folder"], "circle_km": 0} for q in quakes]),
    )
    api.optical_remosaic.main(
        [
            "--store",
            str(tmp_path),
            "--source",
            "exp-a",
            "--source",
            "exp-b",
            "--version",
            "exp-tot",
            "--mosaic",
            "aligned_total",
            "--far_field_km",
            "0.5",
        ]
    )
    out_dir = tmp_path / "products" / "exp-tot" / EVENT
    run = json.loads((out_dir / "run_info.json").read_text())
    assert run["status"] == "COMPLETE", run.get("error")
    structured = run["deliverable_grid"]["alignment"]["dx"]["structured_error"]
    assert set(structured) == {"med30snow", "seas90snow"}
    assert all(r["s_m"] is None for r in structured.values())  # too small a scene: local weights
    ds = gdal.Open(str(out_dir / f"{EVENT}_exp-tot_offsets.tif"))
    assert ds.GetMetadataItem("MOSAIC_METHOD") == "aligned_total"
    assert set(json.loads(ds.GetMetadataItem("STRUCTURED_ERROR_M"))["dy"]) == {
        "med30snow",
        "seas90snow",
    }
    assert np.isfinite(ds.GetRasterBand(3).ReadAsArray()).any()
