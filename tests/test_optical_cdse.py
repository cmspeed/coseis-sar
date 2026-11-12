"""Copernicus Data Space composites: product names, scene choice, masks and the median."""

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
