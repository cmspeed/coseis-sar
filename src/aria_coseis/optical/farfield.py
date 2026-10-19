"""
Far-field reference for optical offsets: the ruptures in a product, and the constant offset of the
product far from them.

Every PRE/POST composite pair is misregistered as a whole by a random amount (about 0.1-0.5 m for
90-day composites, about 1 m for 30-day ones), so every product needs an absolute reference. Without
GNSS, the reference is the far field: the median dx and dy of the nodes more than a set distance from
every rupture, subtracted from the whole product.

The ruptures are those of the catalogued events (the store's events.json) with the epicentre in the
product and the origin within its composite windows: for each, the surface projection of its USGS
finite-fault sub-faults with at least 10 % of the peak slip, or, without a finite-fault model, a
circle of half the surface-rupture length of Wells & Coppersmith (1994, all slip types) around the
epicentre.

The median also absorbs real motion beyond the distance. For a strike-slip rupture covered on both
sides that motion is antisymmetric and largely cancels; for other geometries it biases the reference.
"""

from __future__ import annotations

from typing import Any

import numpy as np

SLIP_FRACTION = 0.1
WELLS_COPPERSMITH_SRL = (-3.22, 0.69)  # log10(surface rupture length, km) = a + b M


def events_in(
    catalog: dict[str, dict[str, Any]],
    bounds_lon_lat: tuple[float, float, float, float],
    t0: str,
    t1: str,
) -> list[dict[str, Any]]:
    """Catalogued events with the epicentre inside (lon0, lat0, lon1, lat1) and origin in [t0, t1]."""
    lon0, lat0, lon1, lat1 = bounds_lon_lat
    out = []
    for folder, m in sorted(catalog.items()):
        lon, lat, t = m.get("longitude"), m.get("latitude"), m.get("time_utc") or ""
        if lon is None or lat is None or not (t0[:19] <= t[:19] <= t1[:19]):
            continue
        if lon0 <= lon <= lon1 and lat0 <= lat <= lat1:
            keys = ("event_id", "magnitude", "time_utc", "longitude", "latitude")
            out.append({"folder": folder, **{k: m.get(k) for k in keys}})
    return out


def circle_radius_km(magnitude: float) -> float:
    """Half the Wells & Coppersmith (1994) surface-rupture length for a magnitude."""
    a, b = WELLS_COPPERSMITH_SRL
    return 0.5 * 10 ** (a + b * float(magnitude))


def rupture_geometry(quakes: list[dict[str, Any]], epsg: int) -> tuple[Any, list[dict[str, Any]]]:
    """Union of the events' ruptures in the given UTM CRS, and where each came from."""
    import requests
    from pyproj import Transformer
    from shapely.geometry import Point, Polygon
    from shapely.ops import transform, unary_union

    from aria_coseis.usgs import get_ffm_geojson_url

    to_utm = Transformer.from_crs(4326, epsg, always_xy=True).transform
    parts, records = [], []
    for q in quakes:
        url = get_ffm_geojson_url(q["event_id"]) if q.get("event_id") else None
        if url:
            feats = requests.get(url, timeout=60).json()["features"]
            peak = max(f["properties"].get("slip", 0) for f in feats)
            keep = [f for f in feats if f["properties"].get("slip", 0) >= SLIP_FRACTION * peak]
            parts += [
                transform(to_utm, Polygon([c[:2] for c in f["geometry"]["coordinates"][0]]))
                for f in keep
            ]
            source = {"finite_fault": url, "subfaults": len(keep)}
        else:
            radius = circle_radius_km(q["magnitude"])
            parts.append(Point(*to_utm(q["longitude"], q["latitude"])).buffer(radius * 1000))
            source = {"circle_km": round(radius, 2)}
        records.append({**q, **source})
    return unary_union(parts), records


def distance_to(geom: Any, x: np.ndarray, y: np.ndarray, step: int = 4) -> np.ndarray:
    """Distance (m) from each node of the x, y node centres to geom, on every step-th node."""
    import shapely

    xs, ys = x[::step], y[::step]
    X, Y = np.meshgrid(xs, ys)
    d = shapely.distance(geom, shapely.points(X.ravel(), Y.ravel())).reshape(X.shape)
    return np.repeat(np.repeat(d, step, axis=0), step, axis=1)[: len(y), : len(x)]


def far_field_shift(
    dx: np.ndarray,
    dy: np.ndarray,
    distance: np.ndarray,
    min_distance_m: float,
    min_nodes: int = 1000,
) -> dict[str, Any]:
    """Median dx and dy of the nodes beyond min_distance_m (both finite); None if too few."""
    far = (distance > min_distance_m) & np.isfinite(dx) & np.isfinite(dy)
    n = int(far.sum())
    if n < min_nodes:
        return {"shift_m": None, "nodes": n}
    return {
        "shift_m": [round(float(np.median(dx[far])), 4), round(float(np.median(dy[far])), 4)],
        "nodes": n,
    }


def to_geojson(geom: Any, epsg: int, records: list[dict[str, Any]]) -> dict[str, Any]:
    """The rupture geometry in longitude/latitude, with the events it came from."""
    from pyproj import Transformer
    from shapely.geometry import mapping
    from shapely.ops import transform

    to_ll = Transformer.from_crs(epsg, 4326, always_xy=True).transform
    return {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "geometry": mapping(transform(to_ll, geom)),
                "properties": {"events": records},
            }
        ],
    }
