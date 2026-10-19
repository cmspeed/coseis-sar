"""Far-field reference: events in a product, rupture fallback, and the far-field shift."""

from typing import Any, Dict

import pytest

import coseis_api as api

np = pytest.importorskip("numpy")
farfield = api.optical_farfield

CATALOG: Dict[str, Dict[str, Any]] = {
    "m_71_ridgecrest": {
        "event_id": "ci38457511",
        "magnitude": 7.1,
        "time_utc": "2019-07-06T03:19:53Z",
        "longitude": -117.6,
        "latitude": 35.77,
    },
    "m_64_ridgecrest": {
        "event_id": "ci38443183",
        "magnitude": 6.4,
        "time_utc": "2019-07-04T17:33:49Z",
        "longitude": -117.5,
        "latitude": 35.70,
    },
    "m_60_elsewhere": {
        "event_id": "x",
        "magnitude": 6.0,
        "time_utc": "2019-07-10T00:00:00Z",
        "longitude": 10.0,
        "latitude": 45.0,
    },
    "m_65_later": {
        "event_id": "y",
        "magnitude": 6.5,
        "time_utc": "2021-01-01T00:00:00Z",
        "longitude": -117.6,
        "latitude": 35.8,
    },
}


def test_events_in_the_product_and_the_windows() -> None:
    found = farfield.events_in(
        CATALOG, (-118.4, 35.2, -117.0, 36.5), "2019-04-07T00:00:00Z", "2019-10-04T00:00:00Z"
    )
    assert [q["folder"] for q in found] == ["m_64_ridgecrest", "m_71_ridgecrest"]


def test_events_without_a_finite_fault_model_get_a_circle(monkeypatch: pytest.MonkeyPatch) -> None:
    pytest.importorskip("shapely")
    pytest.importorskip("pyproj")
    import aria_coseis.usgs as usgs

    monkeypatch.setattr(usgs, "get_ffm_geojson_url", lambda event_id: None)
    assert farfield.circle_radius_km(7.0) == pytest.approx(0.5 * 10 ** (-3.22 + 0.69 * 7.0))
    quake = {"folder": "m_64", **CATALOG["m_64_ridgecrest"]}
    geom, records = farfield.rupture_geometry([quake], 32611)
    radius_m = farfield.circle_radius_km(6.4) * 1000
    assert geom.area == pytest.approx(np.pi * radius_m**2, rel=0.01)
    assert records[0]["circle_km"] == pytest.approx(radius_m / 1000, abs=0.01)


def test_far_field_shift_ignores_the_near_field() -> None:
    rng = np.random.default_rng(0)
    rows, cols = 100, 120
    dist = np.tile(np.arange(cols) * 500.0, (rows, 1))  # distance grows eastward, 500 m per node
    dx = 0.3 + rng.normal(0, 0.1, (rows, cols))
    dy = -0.2 + rng.normal(0, 0.1, (rows, cols))
    dx[:, :20] += 3.0  # near-fault motion within 10 km
    dy[:, :20] -= 2.0
    out = farfield.far_field_shift(dx, dy, dist, 10_000)
    assert out["shift_m"] == pytest.approx([0.3, -0.2], abs=0.02)
    assert out["nodes"] == rows * (cols - 21)
    assert farfield.far_field_shift(dx, dy, dist, 1e9)["shift_m"] is None


def test_distance_to_a_geometry() -> None:
    shapely = pytest.importorskip("shapely")
    x = np.arange(10) * 100.0 + 50
    y = 1000.0 - np.arange(8) * 100.0 - 50
    d = farfield.distance_to(shapely.Point(0, 1000), x, y, step=1)
    assert d.shape == (8, 10)
    assert d[0, 0] == pytest.approx(np.hypot(50, 50))
