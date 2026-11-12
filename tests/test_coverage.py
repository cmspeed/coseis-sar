"""
The coverage check that decides which acquisition dates of a track can be paired (sar/coverage.py),
on synthetic footprints, and its effect on the historic coseismic pair.

The cases are the failures of the former rule (>95% of the track's largest single-date footprint)
found in the archive audit: one date with an extra frame disqualified every regular acquisition of
the track, a small shift of frame edges pushed the nearest date just under 95%, and a track whose
coverage changed at the event got no pair at all.
"""

from datetime import datetime, timezone
from typing import Any, Dict, List

from shapely.geometry import box

import coseis_api as api

AOI = box(0, 0, 1, 1)
ORIGIN = datetime(2026, 7, 10, 12, 0, 0, tzinfo=timezone.utc)
ORIGIN_MS = int(ORIGIN.timestamp() * 1000)


def scene(date: str, width: float, frame: int = 1, x0: float = 0.0) -> Dict[str, Any]:
    """One frame acquired on `date` (YYYY-MM-DD at 10:00 UTC) covering [x0, x0 + width] of the AOI."""
    return {
        "date": f"{date}T10:00:00Z",
        "geometry": box(x0, 0, x0 + width, 1),
        "fileID": f"S1A_{date.replace('-', '')}_{frame}-SLC",
        "frameNumber": frame,
    }


def coseismic_pair(slcs: List[Dict[str, Any]]) -> List[List[str]]:
    """The historic coseismic job's [pre, post] acquisition dates (YYYYMMDD) for these scenes."""
    jobs = api.find_reference_and_secondary_pairs(
        slcs,
        ORIGIN_MS,
        "DESCENDING",
        61,
        "test_event",
        AOI,
        "test0001",
        "coseismic",
        job_list=True,
    )
    return [
        [
            sorted({g.split("_")[1] for g in job["job_parameters"]["secondary_granules"]}),
            sorted({g.split("_")[1] for g in job["job_parameters"]["granules"]}),
        ]
        for job in jobs
    ]


def test_a_date_with_an_extra_frame_does_not_disqualify_the_regular_dates() -> None:
    # Every pass covers 56% of the AOI with one frame; one pass has a second frame (100%).
    slcs = [scene(d, 0.56) for d in ("2026-06-16", "2026-06-28", "2026-07-22", "2026-08-03")]
    slcs += [scene("2026-07-04", 0.56), scene("2026-07-04", 0.44, frame=2, x0=0.56)]
    by_date = {s["date"][:10]: [x for x in slcs if x["date"] == s["date"]] for s in slcs}
    assert api.coverage.usable_dates(by_date, AOI) == set(by_date)
    assert coseismic_pair(slcs) == [[["20260704"], ["20260722"]]]


def test_a_small_shift_of_frame_edges_keeps_the_nearest_date() -> None:
    # The first pass after the origin covers 94.7% of the usual footprint (Ende M6.1 D061, Aug 2026).
    slcs = [scene(d, 1.0) for d in ("2026-06-16", "2026-06-28", "2026-07-04", "2026-07-28")]
    slcs.append(scene("2026-07-16", 0.947))
    assert coseismic_pair(slcs) == [[["20260704"], ["20260716"]]]


def test_a_partial_acquisition_is_skipped() -> None:
    # The first pass after the origin has half the usual coverage: pair with the next full pass.
    slcs = [scene(d, 1.0) for d in ("2026-06-16", "2026-06-28", "2026-07-04", "2026-07-28")]
    slcs.append(scene("2026-07-16", 0.5))
    assert coseismic_pair(slcs) == [[["20260704"], ["20260728"]]]


def test_coverage_that_changes_at_the_event_still_gives_a_pair() -> None:
    # From the event on, every pass covers half of what the earlier passes did (a different
    # acquisition plan, as on us6000k20g D011): pair across the change, over the area both cover.
    slcs = [scene(d, 1.0) for d in ("2026-06-04", "2026-06-16", "2026-06-28", "2026-07-04")]
    slcs += [scene(d, 0.52) for d in ("2026-07-16", "2026-07-28", "2026-08-09")]
    assert coseismic_pair(slcs) == [[["20260704"], ["20260716"]]]


def test_no_dates_or_no_footprint_means_nothing_is_usable() -> None:
    assert api.coverage.typical_footprint([]) == 0.0
    assert not api.coverage.covers(0.0, 0.0)
    assert api.coverage.usable_dates({}, AOI) == set()
    outside = {"2026-07-04": [scene("2026-07-04", 1.0, x0=5.0)]}
    assert api.coverage.usable_dates(outside, AOI) == set()


def test_threshold_is_a_fraction_of_the_typical_footprint() -> None:
    typical = api.coverage.typical_footprint([1.0, 1.0, 1.0, 0.2, 2.0])
    assert typical == 1.0
    assert api.coverage.covers(api.coverage.COVERAGE_FRACTION * typical, typical)
    assert not api.coverage.covers(0.99 * api.coverage.COVERAGE_FRACTION * typical, typical)
