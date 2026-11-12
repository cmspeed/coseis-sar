"""Which acquisition dates of a Sentinel-1 track cover the AOI well enough to be paired.

A date's footprint is the area of the union of its frames inside the AOI. A date is usable when its
footprint is at least COVERAGE_FRACTION of the track's typical footprint, the median over the
track's dates in the search window. The median rather than the largest date: one unusually long
acquisition, or a date with an extra frame, would otherwise disqualify all of the track's regular
acquisitions. COVERAGE_FRACTION leaves room for the small shifts of frame edges between passes,
and still rejects a pass that is missing one of the track's usual frames.

If no date on one side of the earthquake (before or after) reaches that, the track's coverage
changed at the event (a different acquisition plan): that side's dates are measured against their
own typical footprint instead, so the track still gets a pair over the area both sides cover.

Historic pairing (sar/pairing.py) and forward tracking (tracking.py) use the same check.
"""

from __future__ import annotations

from collections.abc import Hashable, Iterable, Mapping
from datetime import datetime
from statistics import median
from typing import Any

from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union

COVERAGE_FRACTION = 0.8


def footprints(
    scenes_by_date: Mapping[Hashable, list[dict[str, Any]]], aoi: BaseGeometry
) -> dict[Hashable, float]:
    """Area of each date's frames inside the AOI (in the AOI's units), keyed like the input."""
    return {
        date: unary_union([scene["geometry"] for scene in scenes]).intersection(aoi).area
        for date, scenes in scenes_by_date.items()
    }


def typical_footprint(areas: Iterable[float]) -> float:
    """The track's typical footprint inside the AOI: the median over its dates (0 without dates)."""
    areas = list(areas)
    return median(areas) if areas else 0.0


def covers(area: float, typical: float) -> bool:
    """True if a date's footprint is at least COVERAGE_FRACTION of the track's typical footprint."""
    return typical > 0 and area >= COVERAGE_FRACTION * typical


def acquired_after(scenes: list[dict[str, Any]], origin: datetime) -> bool:
    """True if a date's first frame was acquired after the origin time (naive UTC)."""
    first = min(scene["date"] for scene in scenes)
    return datetime.strptime(first, "%Y-%m-%dT%H:%M:%SZ") > origin


def usable_dates(
    scenes_by_date: Mapping[Hashable, list[dict[str, Any]]],
    aoi: BaseGeometry,
    origin: datetime | None = None,
) -> set[Hashable]:
    """
    The dates that can be paired (see the module docstring).
    :param scenes_by_date: the track's scenes grouped by acquisition date
    :param aoi: Shapely geometry of the Area of Interest
    :param origin: earthquake origin time (naive UTC); without it, no per-side fallback
    :return: the usable keys of scenes_by_date
    """
    areas = footprints(scenes_by_date, aoi)
    typical = typical_footprint(areas.values())
    usable = {date for date, area in areas.items() if covers(area, typical)}
    if origin is None:
        return usable
    for after in (False, True):
        side = {
            date: area
            for date, area in areas.items()
            if acquired_after(scenes_by_date[date], origin) == after
        }
        if side and not usable & side.keys():
            side_typical = typical_footprint(side.values())
            usable |= {date for date, area in side.items() if covers(area, side_typical)}
    return usable
