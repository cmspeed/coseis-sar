"""Copernicus Data Space composites: product names, scene choice, masks and the median."""

import json
from datetime import datetime, timezone

import pytest

import coseis_api as api

np = pytest.importorskip("numpy")
cdse = api.optical_cdse
store = api.optical_store


def product(name: str) -> dict:
    fields = cdse.parse_name(name)
    assert fields is not None
    return fields


def test_parse_name_reads_the_fields():
    p = product("S2B_MSIL1C_20190711T183929_N0500_R070_T11SLA_20230521T191842.SAFE")
    assert p["name"] == "S2B_MSIL1C_20190711T183929_N0500_R070_T11SLA_20230521T191842"
    assert p["level"] == "1C"
    assert p["baseline"] == "05.00"
    assert p["orbit"] == 70
    assert p["tile"] == "11SLA"
    assert p["sensing"] == datetime(2019, 7, 11, 18, 39, 29, tzinfo=timezone.utc)
    assert cdse.parse_name("not a product") is None


def test_choose_scenes_keeps_one_reprocessed_product_per_tile_and_date():
    names = [
        "S2A_MSIL1C_20171120T074209_N0206_R092_T38SNC_20171120T093513",  # original: left out
        "S2B_MSIL1C_20171120T074209_N0500_R092_T38SNC_20231015T220322",
        "S2B_MSIL1C_20171120T074209_N0500_R092_T38SNC_20231101T000000",  # later processing wins
        "S2B_MSIL1C_20171120T074209_N0500_R092_T38SND_20231015T220322",
        "S2A_MSIL1C_20171125T074209_N0500_R092_T38SNC_20231015T220322",
    ]
    chosen = cdse.choose_scenes([product(n) for n in names], "05.00")
    assert [(p["tile"], p["sensing"].day, p["processed"]) for p in chosen] == [
        ("38SNC", 20, "20231101T000000"),
        ("38SND", 20, "20231015T220322"),
        ("38SNC", 25, "20231015T220322"),
    ]


def test_harmonize_removes_the_offset_and_keeps_nodata():
    dn = np.array([0, 500, 1000, 1001, 3000], dtype=np.uint16)
    assert cdse.harmonize(dn).tolist() == [0, 1, 1, 1, 2000]


def test_snow_rule_is_the_earth_engine_rule():
    rule = store.COMPOSITE_RECIPES["med90snowcdse"]["snow_mask"]
    b3 = np.array([5000, 5000, 900, 5000], dtype=np.uint16)
    b11 = np.array([1000, 4000, 100, 1000], dtype=np.uint16)
    b8 = np.array([4000, 4000, 4000, 1000], dtype=np.uint16)
    # snow; NDSI too low; B3 too dark; B8 too dark
    assert cdse.snow(b3, b8, b11, rule).tolist() == [True, False, False, False]


def test_cloud_mask_dilates_the_masked_classes():
    scl = np.full((7, 7), 4, dtype=np.uint8)
    scl[3, 3] = 9
    mask = cdse.cloud_mask(scl, store.CDSE_SCL_MASKED, 2)
    assert mask.sum() == 25 and mask[1:6, 1:6].all()
    assert cdse.cloud_mask(scl, store.CDSE_SCL_MASKED, 0).sum() == 1


def test_cloud_mask_adds_the_cloud_probability():
    scl = np.full((7, 7), 3, dtype=np.uint8)  # cloud shadow everywhere: not a masked class here
    probability = np.full((7, 7), 10, dtype=np.uint8)
    probability[3, 3] = 15
    assert cdse.cloud_mask(scl, [0, 1], 0, probability, 15).sum() == 1
    assert cdse.cloud_mask(scl, [0, 1], 1, probability, 15).sum() == 9
    assert cdse.cloud_mask(scl, [0, 1], 0, probability, 16).sum() == 0


def test_nan_median_ignores_zeros_and_counts_dates():
    stack = np.array([[[0, 10]], [[0, 20]], [[0, 25]], [[0, 0]]], dtype=np.uint16)
    med, count = cdse.nan_median(stack)
    assert med.tolist() == [[0, 20]]
    assert count.tolist() == [[0, 3]]
    med2, _ = cdse.nan_median(np.array([[[11]], [[12]]], dtype=np.uint16))
    assert med2.tolist() == [[11]]  # 11.5 truncated, as Earth Engine's toUint16


def test_cdse_recipes_mirror_the_earth_engine_ones():
    for rid in ("med30snow", "med90snow", "seas90snow"):
        ee, ds = store.COMPOSITE_RECIPES[rid], store.COMPOSITE_RECIPES[f"{rid}cdse"]
        assert ds["window_days"] == ee["window_days"]
        assert ds["snow_mask"] == ee["snow_mask"] and ds["band"] == ee["band"] == "B8"
        assert ds["processing_baseline_min"] == "05.00" and "cloud_score_min" not in ds
        assert ds["builder"] == "aria_coseis.optical.cdse" and "builder" not in ee
    assert store.COMPOSITE_RECIPES["seas90snowcdse"]["stage_from_recipe"] == {
        "POST": "med90snowcdse"
    }


def test_cdseprob_recipes_change_only_the_cloud_mask():
    for rid in ("med30snow", "med90snow", "seas90snow"):
        scl, prob = (
            store.COMPOSITE_RECIPES[f"{rid}cdse"],
            store.COMPOSITE_RECIPES[f"{rid}cdseprob"],
        )
        same = {k for k in scl if k not in ("cloud_mask", "description", "stage_from_recipe")}
        assert all(prob[k] == scl[k] for k in same)
        assert cdse.recipe_mask(f"{rid}cdseprob") == store.CDSE_PROBABILITY_MASK
        assert 3 not in cdse.recipe_mask(f"{rid}cdseprob")["scl"]
    assert store.COMPOSITE_RECIPES["seas90snowcdseprob"]["stage_from_recipe"] == {
        "POST": "med90snowcdseprob"
    }
    # the *cdse date layers keep their folder; the new mask gets its own
    assert cdse.mask_signature("med90snowcdse") == "dates_scl0138910d2_snow0.4-1100-1000"
    assert cdse.mask_signature("med90snowcdseprob") == "dates_scl01d5_prb15_snow0.4-1100-1000"


def test_choose_scenes_keeps_only_refined_products_when_asked():
    names = [
        "S2A_MSIL1C_20230314T082721_N0510_R021_T36SYG_20240820T155416",  # REFINED, earlier
        "S2A_MSIL1C_20230314T082721_N0510_R021_T36SYG_20240823T001122",  # NOT_REFINED, later
        "S2A_MSIL1C_20230324T082721_N0510_R021_T36SYG_20240901T000000",  # only NOT_REFINED
    ]
    products = [cdse.parse_name(n) for n in names]
    refined = dict(zip(names, ["REFINED", "NOT_REFINED", "NOT_REFINED"]))
    assert [p["name"] for p in cdse.choose_scenes(products, "05.00")] == names[1:]
    assert [p["name"] for p in cdse.choose_scenes(products, "05.00", refined)] == names[:1]


def test_refining_flag_is_read_from_the_start_of_the_datastrip_metadata():
    xml = b"<a>" + b"x" * 200_000 + b'<Image_Refining flag="NOT_REFINED"/>' + b"y" * 10
    chunks = [xml[i : i + 65536] for i in range(0, len(xml), 65536)]
    assert cdse.refining_flag(iter(chunks)) == "NOT_REFINED"
    assert cdse.refining_flag(iter([b"<a/>"])) == "UNKNOWN"


def test_cdseref_recipes_require_refined_products_and_keep_their_own_date_layers():
    for rid in ("med30snow", "med90snow", "seas90snow"):
        ref, base = store.COMPOSITE_RECIPES[f"{rid}cdseref"], store.COMPOSITE_RECIPES[f"{rid}cdse"]
        assert ref["require_refined"] and "require_refined" not in base
        assert (
            ref["window_days"] == base["window_days"] and ref["cloud_mask"] == base["cloud_mask"]
        )
    assert store.COMPOSITE_RECIPES["seas90snowcdseref"]["stage_from_recipe"] == {
        "POST": "med90snowcdseref"
    }
    assert cdse.mask_signature("med90snowcdseref") == (
        "dates_scl0138910d2_snow0.4-1100-1000_refined"
    )


def test_scene_dir_keeps_another_product_of_the_same_tile_and_date_apart(tmp_path):
    first = "S2A_MSIL1C_20230314T082721_N0510_R021_T36SYG_20240823T001122"
    second = "S2A_MSIL1C_20230314T082721_N0510_R021_T36SYG_20240820T155416"
    folder = cdse.scene_dir(tmp_path, {"l1c": cdse.parse_name(first)})
    assert folder.name == "36SYG_20230314T082721"
    folder.mkdir(parents=True)
    (folder / "done.json").write_text(json.dumps({"l1c": first, "l2a": "x"}))
    assert cdse.scene_dir(tmp_path, {"l1c": cdse.parse_name(first)}) == folder
    other = cdse.scene_dir(tmp_path, {"l1c": cdse.parse_name(second)})
    assert other.name == "36SYG_20230314T082721_20240820T155416"


class BrokenOnce:
    """A response whose body breaks after the first chunk, the first time it is fetched."""

    calls = 0

    def __init__(self) -> None:
        BrokenOnce.calls += 1
        self.broken = BrokenOnce.calls == 1

    def __enter__(self) -> "BrokenOnce":
        return self

    def __exit__(self, *exc: object) -> None:
        return None

    def iter_content(self, size: int):  # noqa: ANN201
        yield b"part"
        if self.broken:
            raise cdse.requests.exceptions.ChunkedEncodingError("Connection reset by peer")
        yield b"-rest"


class FakeSession:
    def get(self, url: str, stream: bool = False) -> BrokenOnce:
        return BrokenOnce()


def test_a_download_broken_mid_stream_is_fetched_again(tmp_path, monkeypatch):
    monkeypatch.setattr(cdse.time, "sleep", lambda s: None)
    BrokenOnce.calls = 0
    product = {"id": "abc", "name": "S2A_MSIL1C_X"}
    cdse._fetch(FakeSession(), product, ["GRANULE"], tmp_path / "B08.jp2")
    assert (tmp_path / "B08.jp2").read_bytes() == b"part-rest"
    assert BrokenOnce.calls == 2
    assert not (tmp_path / "B08.jp2.part").exists()


def test_build_exits_with_an_error_when_an_event_fails(tmp_path, monkeypatch):
    def fail(*args: object) -> dict:
        raise ConnectionResetError("Connection reset by peer")

    monkeypatch.setattr(cdse, "build_event", fail)
    monkeypatch.setattr(cdse.store, "require_store", lambda root: tmp_path)
    monkeypatch.setattr(cdse.store, "code_provenance", lambda: {})
    monkeypatch.setattr(
        cdse.runs, "code_signature", lambda p: {"aria_coseis": {"commit": "x", "dirty": False}}
    )
    with pytest.raises(SystemExit) as raised:
        cdse.main(["build", "--recipe", "med90snowcdse", "--event", "e1", "--event", "e2"])
    assert "2 of 2 events failed" in str(raised.value)
