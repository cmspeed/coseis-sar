"""Historic optical mode (USGS + scene catalogs replayed): rake filtering, FFM AOI buffering, pairing."""

import glob
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict

import pytest
import time_machine
from shapely.geometry import shape

import coseis_api as api
from conftest import read_json

NOW = datetime(2026, 10, 5, 12, 0, 0, tzinfo=timezone.utc)
KAHRAMANMARAS_DATE = (
    "2023-02-06"  # M7.8 and M7.5 strike-slip events (us6000jllz, us6000jlqa) plus aftershocks
)


def collect_outputs(workdir: Path) -> Dict[str, Any]:
    def all_json(pattern: str) -> Dict[str, Any]:
        return {
            Path(p).name.split("_20")[0] if "jobs_list" in p else Path(p).name: read_json(Path(p))
            for p in sorted(glob.glob(str(workdir / "outputs" / pattern)))
        }

    return {
        "jobs": all_json("jobs_list_*.json"),
        "aois": all_json("*_AOI.geojson"),
        "significant": read_json(
            workdir / "outputs" / f"significant_earthquakes_{KAHRAMANMARAS_DATE}.geojson"
        ),
    }


@pytest.mark.vcr
def test_sentinel2_copernicus_job_list(workdir: Path, golden: Callable[[str, Any], None]) -> None:
    with time_machine.travel(NOW, tick=False):
        api.main_historic(
            start_date=KAHRAMANMARAS_DATE,
            job_list=True,
            resolution=30,
            sensor="sentinel-2",
            optical_backend="copernicus",
            optical_level="toa",
        )
    outputs = collect_outputs(workdir)

    # Only strike-slip events pass the optical rake filter
    kept = sorted(Path(f["properties"]["url"]).name for f in outputs["significant"]["features"])
    assert kept == ["us6000jllz", "us6000jlqa"]
    # Both events have a finite-fault model; for optical the FFM is buffered 0.15 deg and boxed
    assert len(outputs["aois"]) == 2
    for aoi in outputs["aois"].values():
        polygon = shape(aoi)
        assert polygon.equals(polygon.envelope)
    golden("historic_optical/kahramanmaras_sentinel2_copernicus", outputs)
