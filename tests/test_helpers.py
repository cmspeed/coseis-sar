"""Pure helpers: no network, no files outside the working directory."""

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import coseis_api as api


def test_to_snake_case_strips_accents_and_punctuation() -> None:
    assert api.to_snake_case("M 7.8 - Pazarcık, Türkiye") == "m_78_pazarck_turkiye"
    assert (
        api.to_snake_case("M 5.6 - 94 km SW of Tamarindo, Costa Rica")
        == "m_56_94_km_sw_of_tamarindo_costa_rica"
    )
    assert api.to_snake_case("  M 6.0 -  Ende  ") == "m_60_ende"


def test_convert_time_drops_milliseconds_and_is_utc() -> None:
    assert api.convert_time(1755208701564) == datetime(
        2025, 8, 14, 21, 58, 21, tzinfo=timezone.utc
    )


def test_make_aoi_is_one_degree_box_centered_on_epicenter() -> None:
    aoi = api.make_aoi([121.5, -8.25, 10.0])
    assert aoi.bounds == (121.0, -8.75, 122.0, -7.75)
    assert aoi.area == 1.0


def test_generate_pairs_modes() -> None:
    scenes = ["a", "b", "c"]
    assert api.generate_pairs(scenes, "sequential") == [("a", "b"), ("b", "c")]
    assert api.generate_pairs(scenes, "all") == [("a", "b"), ("a", "c"), ("b", "c")]
    assert api.generate_pairs(scenes, "coseismic") == []


def test_make_job_json_shape(golden: Callable[[str, Any], None]) -> None:
    job = api.make_job_json(
        "m_71_test_event",
        "us0000test",
        "D",
        121,
        ["S1A_REF_1", "S1A_REF_2"],
        ["S1A_SEC_1"],
        30,
    )
    golden("helpers/make_job_json", job)


def test_make_optical_job_json_shape(golden: Callable[[str, Any], None]) -> None:
    job = api.make_optical_job_json(
        "m_78_test_event",
        "us0000test",
        "021",
        "2023-01-28",
        "2023-02-12",
        ["S2A_REF"],
        ["S2B_SEC"],
        pre_cc=1.234,
        post_cc=5.678,
        zone_suffix="Z37",
    )
    golden("helpers/make_optical_job_json", job)


def test_ascii_table_to_html() -> None:
    table = "+---+---+\n| a | b |\n+---+---+\n| 1 | 2 |\n+---+---+"
    html = api.ascii_table_to_html(table)
    assert html.count("<tr>") == 2
    assert ">a</th>" in html and ">2</td>" in html
    assert api.ascii_table_to_html("No passes.").startswith("<pre")


def test_parse_custom_eq_list(tmp_path: Path, golden: Callable[[str, Any], None]) -> None:
    entries = [
        {
            "title": "M 7.1 - 2025 Southern Tibetan Plateau Earthquake",
            "epicenter": {
                "latitude": 28.639,
                "longitude": 87.361,
                "depth_km": 10.0,
                "usgs_event_url": "https://earthquake.usgs.gov/earthquakes/eventpage/us6000pi9w",
            },
            "time": "2025-01-07 01:05:16 UTC",
        },
        {"title": "missing coordinates", "epicenter": {}, "time": "2025-01-07 01:05:16 UTC"},
        {"title": "bad time", "epicenter": {"latitude": 1, "longitude": 2}, "time": "yesterday"},
    ]
    path = tmp_path / "eq_list.json"
    path.write_text(json.dumps(entries))
    golden("helpers/parse_custom_eq_list", api.parse_custom_eq_list(str(path)))
