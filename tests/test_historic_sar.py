"""Historic SAR mode end to end (USGS + ASF replayed): job lists, AOIs and earthquake info."""

import glob
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict

import pytest
import time_machine

import coseis_api as api
from conftest import read_json

# Output file names contain the current time; freeze it.
NOW = datetime(2026, 10, 5, 12, 0, 0, tzinfo=timezone.utc)

EVENT_DATES = {
    "tibet_2025": "2025-01-07",  # M7.1 Southern Tibetan Plateau (us6000pi9w), has a finite-fault model
    "ende_2026": "2026-08-14",  # M7.7 NNW of Ende, Indonesia (us6000tkt2)
}


def collect_outputs(workdir: Path) -> Dict[str, Any]:
    def only(pattern: str) -> Any:
        matches = glob.glob(str(workdir / pattern))
        assert len(matches) == 1, f"expected one {pattern}, found {matches}"
        return read_json(Path(matches[0]))

    return {
        "jobs": only("jobs_list_*.json"),
        "earthquake_info": only("earthquake_info_*.json"),
        "aois": {
            Path(p).name: read_json(Path(p))
            for p in sorted(glob.glob(str(workdir / "*_AOI.geojson")))
        },
    }


@pytest.mark.vcr
@pytest.mark.parametrize("event", sorted(EVENT_DATES))
def test_hyp3_job_list(event: str, workdir: Path, golden: Callable[[str, Any], None]) -> None:
    with time_machine.travel(NOW, tick=False):
        api.main_historic(
            start_date=EVENT_DATES[event], pairing_mode="coseismic", job_list=True, resolution=30
        )
    outputs = collect_outputs(workdir)
    assert outputs["jobs"], "expected at least one HyP3 job"
    golden(f"historic_sar/{event}_job_list", outputs)


@pytest.mark.vcr
def test_without_job_list_writes_pair_json_and_runs_nothing(
    workdir: Path, golden: Callable[[str, Any], None], topsapp: Dict[str, Any]
) -> None:
    with time_machine.travel(NOW, tick=False):
        api.main_historic(
            start_date=EVENT_DATES["tibet_2025"],
            pairing_mode="coseismic",
            job_list=False,
            resolution=30,
        )
    outputs = collect_outputs(workdir)
    # Without --job_list, historic mode writes topsApp pair JSON plus frame map files; it does not run topsApp
    assert topsapp["calls"] == []
    assert glob.glob(str(workdir / "*_frames.geojson"))
    golden("historic_sar/tibet_2025_pair_json", outputs)
