"""
Forward mode: discovery (GitHub Actions run), processing (local cron run) and the tracker state machine.

Inputs and expected states come from real production commits on `main` (tests/fixtures/production/),
so these tests also check that this branch reproduces what `main` did.
"""

import glob
import os
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, List

import pytest
import time_machine

import coseis_api as api
from conftest import read_json

PRODUCTION = Path(__file__).parent / "fixtures" / "production"
LOCK_FILE = Path("/tmp/coseis_processing.lock")

# When production ran: the GitHub Action that found Tamarindo, and the cron run that processed Ende D61
TAMARINDO_DISCOVERY = datetime(2026, 10, 1, 0, 5, 0, tzinfo=timezone.utc)
TAMARINDO_NO_POST_DATA = datetime(2026, 10, 1, 4, 0, 0, tzinfo=timezone.utc)
ENDE_PROCESSING = datetime(2026, 10, 2, 3, 55, 3, tzinfo=timezone.utc)


@pytest.fixture
def no_lock() -> None:
    """main_forward uses a fixed lock file; don't run while a real forward run holds it."""
    if LOCK_FILE.exists():
        pytest.skip(f"{LOCK_FILE} exists (a forward run may be active)")
    yield
    assert not LOCK_FILE.exists(), "main_forward left its lock file behind"


def seed(workdir: Path, fixture: str) -> None:
    """Copy a production tracker state into the working directory (tracker JSONs + partial job files)."""
    (workdir / "active_jobs" / "partials").mkdir(parents=True, exist_ok=True)
    for path in (PRODUCTION / fixture).glob("*.json"):
        target = (
            workdir
            / ("active_jobs/partials" if path.name.startswith("job_") else "active_jobs")
            / path.name
        )
        shutil.copy(path, target)


def tracker_files(workdir: Path) -> Dict[str, Any]:
    return {
        Path(p).name: read_json(Path(p))
        for p in sorted(glob.glob(str(workdir / "active_jobs" / "*.json")))
    }


def partial_files(workdir: Path) -> Dict[str, Any]:
    return {
        Path(p).name: read_json(Path(p))
        for p in sorted(
            glob.glob(str(workdir / "active_jobs" / "partials" / "job_*_partial.json"))
        )
    }


def without_event_id(jobs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """`main` predates the event_id field on job JSON (added on develop)."""
    return [{k: v for k, v in job.items() if k != "event_id"} for job in jobs]


@pytest.mark.vcr
def test_discovery_matches_production(
    workdir: Path,
    no_lock: None,
    fake_next_pass: Dict[str, Any],
    sent_emails: List[Dict[str, Any]],
    golden: Callable[[str, Any], None],
) -> None:
    """The GitHub Actions run: `coseis.py --forward --pairing coseismic --send_email`."""
    with time_machine.travel(TAMARINDO_DISCOVERY, tick=False):
        api.main_forward(pairing_mode="coseismic", send_email_flag=True)

    trackers = tracker_files(workdir)
    partials = partial_files(workdir)
    golden("forward/discovery_2026-10-01T0005Z", {"trackers": trackers, "partials": partials})

    # Same tracker entry and partial jobs as production commit ff08a73
    assert trackers["us6000tymj.json"] == read_json(
        PRODUCTION / "tamarindo_ff08a73" / "us6000tymj.json"
    )
    for name in (
        "job_m_56_94_km_sw_of_tamarindo_costa_rica_ASCENDING_165_partial.json",
        "job_m_56_94_km_sw_of_tamarindo_costa_rica_DESCENDING_157_partial.json",
    ):
        assert without_event_id(partials[name]) == read_json(
            PRODUCTION / "tamarindo_ff08a73" / name
        )
        assert partials[name][0]["event_id"] == "us6000tymj"

    # New-event email to primary recipients, linking the overpass map published under docs/maps/
    tamarindo = [m for m in sent_emails if "Tamarindo" in m["subject"]]
    assert len(tamarindo) == 1
    assert tamarindo[0]["subject"] == "New Event: M 5.6 - 94 km SW of Tamarindo, Costa Rica"
    assert tamarindo[0]["bcc"] == ["primary@example.com"]
    map_name = "m_56_94_km_sw_of_tamarindo_costa_rica_us6000tymj_overpass_map.html"
    assert (workdir / "docs" / "maps" / map_name).exists()
    assert f"https://cmspeed.github.io/coseis-sar/maps/{map_name}" in tamarindo[0]["contents"][0]


@pytest.mark.vcr
def test_processing_waits_without_post_seismic_data(
    workdir: Path, no_lock: None, topsapp: Dict[str, Any]
) -> None:
    """The local cron run before any post-event SLC exists: nothing changes."""
    seed(workdir, "tamarindo_ff08a73")
    before = (tracker_files(workdir), partial_files(workdir))
    with time_machine.travel(TAMARINDO_NO_POST_DATA, tick=False):
        api.main_forward(
            pairing_mode="coseismic", resolution=30, do_processing=True, process_only=True
        )
    assert (tracker_files(workdir), partial_files(workdir)) == before
    assert topsapp["calls"] == []


def run_processing(workdir: Path, topsapp: Dict[str, Any], outcome: str) -> Dict[str, Any]:
    """The local cron run (`--forward --pairing coseismic --resolution 30 --do_processing --process_only`)."""
    seed(workdir, "ende_98ea2b5")
    topsapp["outcome"] = outcome
    with time_machine.travel(ENDE_PROCESSING, tick=False):
        api.main_forward(
            pairing_mode="coseismic", resolution=30, do_processing=True, process_only=True
        )
    return tracker_files(workdir)["us6000tkt2.json"]["tracks"]["DESCENDING_61"]


@pytest.mark.vcr
def test_processing_success_matches_production(
    workdir: Path, no_lock: None, topsapp: Dict[str, Any]
) -> None:
    track = run_processing(workdir, topsapp, "success")

    # Same state as production commit 92edeb3, apart from the machine-specific processing location
    expected = read_json(PRODUCTION / "ende_92edeb3" / "us6000tkt2.json")["tracks"][
        "DESCENDING_61"
    ]
    pair_dir = (
        "m_77_68_km_nnw_of_ende_indonesia/DESCENDING061/coseismic/DESCENDING061_20260802_20260820"
    )
    assert expected["location"].endswith(pair_dir)
    assert track["location"] == str(Path(api.settings.root_dir) / pair_dir)
    assert {k: v for k, v in track.items() if k != "location"} == {
        k: v for k, v in expected.items() if k != "location"
    }

    # topsApp ran once on the pre/post pair; its resolution comes from the partial job file
    (call,) = topsapp["calls"]
    assert call["cwd"] == track["location"]
    cmd = call["cmd"]
    assert cmd[cmd.index("--reference-scenes") + 1] == (
        "S1C_IW_SLC__1SDV_20260820T213528_20260820T213558_009084_01208A_0DB0 "
        "S1C_IW_SLC__1SDV_20260820T213556_20260820T213631_009084_01208A_BC33"
    )
    assert cmd[cmd.index("--output-resolution") + 1] == "90"

    # Product kept, raw SLCs and intermediate dirs cleaned up, partial file consumed
    location = Path(track["location"])
    assert (location / "S1-GUNW-test-product.nc").exists()
    assert not (location / "S1A_IW_SLC__fake.zip").exists()
    assert not (location / "fine_interferogram").exists()
    assert (
        location / "job_m_77_68_km_nnw_of_ende_indonesia_DESCENDING_61_COMPLETED.json"
    ).exists()
    assert partial_files(workdir) == {}


@pytest.mark.vcr
@pytest.mark.parametrize(
    "outcome, status_prefix",
    [
        (
            "nonzero_exit",
            "Failed: Command '['conda', 'run', '-n', 'topsapp_env_trappist_python11'",
        ),
        ("missing_conda", "Failed: conda"),
    ],
)
def test_processing_failure_is_flagged(
    outcome: str,
    status_prefix: str,
    workdir: Path,
    no_lock: None,
    topsapp: Dict[str, Any],
    sent_emails: List[Dict[str, Any]],
) -> None:
    track = run_processing(workdir, topsapp, outcome)
    assert track["status"] == "READY_FOR_EMAIL"
    assert track["processing_status"].startswith(status_prefix)
    assert partial_files(workdir) == {}

    # The next GitHub Actions run emails the result and parks the track for a human
    api.check_tracker_for_updates(do_processing=False, send_email_flag=True)
    track = tracker_files(workdir)["us6000tkt2.json"]["tracks"]["DESCENDING_61"]
    assert track["status"] == "FAILED_NEEDS_ATTENTION"
    (email,) = sent_emails
    assert email["subject"] == "PROCESSING COMPLETED (WITH FAILURES): 1 Jobs Processed"
    assert email["bcc"] == ["secondary@example.com"]


def test_email_run_reports_and_removes_finished_event(
    workdir: Path, sent_emails: List[Dict[str, Any]]
) -> None:
    """The GitHub Actions run after processing (production commit 595ed24 removed the tracker file)."""
    seed(workdir, "ende_92edeb3")
    api.check_tracker_for_updates(do_processing=False, send_email_flag=True)
    assert tracker_files(workdir) == {}
    (email,) = sent_emails
    assert email["subject"] == "PROCESSING COMPLETED (SUCCESS): 1 Jobs Processed"
    assert email["bcc"] == ["secondary@example.com"]
    assert "DESCENDING_61" in email["contents"][0]


def test_tracker_states_left_alone_by_the_other_runner(
    workdir: Path, sent_emails: List[Dict[str, Any]], topsapp: Dict[str, Any]
) -> None:
    """The email runner ignores AWAITING tracks; the processing runner ignores READY_FOR_EMAIL tracks."""
    seed(workdir, "tamarindo_ff08a73")
    before = tracker_files(workdir)
    api.check_tracker_for_updates(do_processing=False, send_email_flag=True)
    assert tracker_files(workdir) == before and sent_emails == []

    shutil.rmtree(workdir / "active_jobs")
    seed(workdir, "ende_92edeb3")
    before = tracker_files(workdir)
    api.check_tracker_for_updates(do_processing=True, send_email_flag=False)
    assert tracker_files(workdir) == before and topsapp["calls"] == [] and sent_emails == []


def test_lock_file_blocks_a_second_run(workdir: Path, no_lock: None, http_log: List[str]) -> None:
    LOCK_FILE.write_text("running")
    try:
        api.main_forward(pairing_mode="coseismic", send_email_flag=True)
        assert http_log == []
        assert not (workdir / "active_jobs").exists()
    finally:
        LOCK_FILE.unlink()


def dead_pid() -> int:
    """Process ID of a process that has already exited."""
    import subprocess
    import sys

    # Popen, not subprocess.run: run() is replaced by the topsApp fake
    proc = subprocess.Popen([sys.executable, "-c", "pass"])
    proc.wait()
    return proc.pid


@pytest.mark.parametrize(
    "lock_content, age_hours, removed",
    [
        (f"{dead_pid()} 2026-10-07T00:00:00+00:00\n", 0, True),  # owner process is gone: stale
        (f"{os.getpid()} 2026-10-07T00:00:00+00:00\n", 0, False),  # owner still running: held
        ("running", 13, True),  # legacy lock, older than 12 h: stale
        ("running", 1, False),  # legacy lock, recent: held
    ],
)
def test_stale_lock_is_removed_and_live_lock_blocks(
    lock_content: str,
    age_hours: float,
    removed: bool,
    workdir: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture,
) -> None:
    lock = workdir / "forward.lock"
    lock.write_text(lock_content)
    old = time.time() - age_hours * 3600
    os.utime(lock, (old, old))
    monkeypatch.setattr(api.settings, "LOCK_FILE", str(lock))

    api.main_forward(pairing_mode="coseismic", process_only=True)  # no tracker: no network needed

    out = capsys.readouterr().out
    if removed:
        assert "Removing stale lock file" in out
        assert not lock.exists()  # removed, then the run's own lock released at the end
    else:
        assert "Previous processing run still active" in out
        assert lock.read_text() == lock_content
