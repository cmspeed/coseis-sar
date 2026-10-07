"""check_significance: magnitude/depth/land filters by mode, and the optical-only rake filter."""

from typing import Any, Dict, List

import pytest

import coseis_api as api

# Synthetic events. "land" points are on or near land; "ocean" is mid-Atlantic, far from any coast.
LAND = (-85.9, 10.1)  # Nicoya Peninsula, Costa Rica
OCEAN = (-30.0, 0.0)


def make_event(
    event_id: str, mag: float, depth: float, where=LAND, alert: str = None
) -> Dict[str, Any]:
    return {
        "id": event_id,
        "title": f"M {mag} - synthetic {event_id}",
        "mag": mag,
        "alert": alert,
        "time": 1735689600000,
        "coordinates": [where[0], where[1], depth],
        "place": f"synthetic {event_id}",
        "url": f"https://earthquake.usgs.gov/earthquakes/eventpage/{event_id}",
    }


SYNTHETIC = [
    make_event("m65_d10", 6.5, 10.0),
    make_event("m65_d45", 6.5, 45.0),
    make_event("m60_d40", 6.0, 40.0),
    make_event("m59_d10", 5.9, 10.0),
    make_event("m55_d15", 5.5, 15.0),
    make_event("m56_d16", 5.6, 16.0),
    make_event("m65_d10_ocean", 6.5, 10.0, where=OCEAN),
    make_event("m65_d10_alert_green", 6.5, 10.0, alert="green"),
    make_event("no_mag", None, 10.0),
]


def accepted_ids(result: List[Dict[str, Any]]) -> List[str]:
    return sorted(eq["id"] for eq in (result or []))


def assert_only_coastline_requests(urls: List[str]) -> None:
    """SAR significance must not query USGS event details (rake)."""
    assert urls, "expected the coastline request"
    assert all("coastline.geojson" in url for url in urls), urls


@pytest.mark.vcr
def test_historic_sar(workdir, http_log) -> None:
    result = api.check_significance(
        [dict(e) for e in SYNTHETIC], "2025-01-01", sensor="sar", mode="historic"
    )
    # M>=6.0, depth <=40 km, near land; USGS alert level is not required
    assert accepted_ids(result) == ["m60_d40", "m65_d10", "m65_d10_alert_green"]
    assert_only_coastline_requests(http_log)


@pytest.mark.vcr
def test_forward_sar(workdir, http_log) -> None:
    result = api.check_significance(
        [dict(e) for e in SYNTHETIC], "2025-01-01", sensor="sar", mode="forward"
    )
    # (M>=5.5 and <=15 km) or (M>=6.0 and <=40 km), near land
    assert accepted_ids(result) == [
        "m55_d15",
        "m59_d10",
        "m60_d40",
        "m65_d10",
        "m65_d10_alert_green",
    ]
    assert_only_coastline_requests(http_log)


@pytest.mark.vcr
def test_historic_sar_writes_csv_and_geojson(workdir) -> None:
    api.check_significance(
        [make_event("m65_d10", 6.5, 10.0)], "2025-01-01", sensor="sar", mode="historic"
    )
    assert (workdir / "outputs" / "significant_earthquakes_2025-01-01.csv").exists()
    assert (workdir / "outputs" / "significant_earthquakes_2025-01-01.geojson").exists()


@pytest.mark.vcr
def test_optical_requires_strike_slip_rake(workdir) -> None:
    events = [
        # 2023 Kahramanmaras M7.8, strike-slip: kept
        make_event("us6000jllz", 7.8, 10.0, where=(37.014, 37.226)),
        # 2025 Tibet M7.1, normal faulting: rejected by rake
        make_event("us6000pi9w", 7.1, 10.0, where=(87.361, 28.639)),
        # Unknown event id: no rake available, rejected
        make_event("us0000none", 7.0, 10.0),
    ]
    for sensor in ("sentinel-2", "landsat"):
        result = api.check_significance(
            [dict(e) for e in events], "2025-01-01", sensor=sensor, mode="historic"
        )
        assert accepted_ids(result) == ["us6000jllz"]
        assert result[0]["rakes"]
