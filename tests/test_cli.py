"""Command-line routing and validation of scripts/coseis.py."""
from typing import Any, Dict, List

import pytest

import coseis_api as api


def single_call(argv: List[str]) -> Dict[str, Any]:
    calls = api.run_cli(argv)
    assert len(calls) == 1, calls
    return calls[0]


def test_github_actions_command() -> None:
    """.github/workflows/coseis-cron.yml"""
    call = single_call(["--forward", "--pairing", "coseismic", "--send_email"])
    assert call == {"func": "main_forward", "args": [], "kwargs": {
        "pairing_mode": "coseismic", "resolution": 30, "do_processing": False,
        "send_email_flag": True, "process_only": False,
    }}


def test_local_cron_command() -> None:
    """scripts/run_coseis_forward.sh"""
    call = single_call(["--forward", "--pairing", "coseismic", "--resolution", "30", "--do_processing", "--process_only"])
    assert call == {"func": "main_forward", "args": [], "kwargs": {
        "pairing_mode": "coseismic", "resolution": 30, "do_processing": True,
        "send_email_flag": False, "process_only": True,
    }}


def test_historic_sar_job_list() -> None:
    call = single_call(["--historic", "--dates", "2014-10-01", "2026-07-31", "--pairing", "coseismic", "--job_list"])
    assert call == {"func": "main_historic", "args": [], "kwargs": {
        "start_date": "2014-10-01", "end_date": "2026-07-31", "aoi": None, "pairing_mode": "coseismic",
        "job_list": True, "resolution": 30, "sensor": "sar", "optical_backend": "copernicus", "optical_level": "toa",
    }}


def test_historic_optical_gee() -> None:
    call = single_call(["--historic", "--dates", "2023-02-06", "--sensor", "landsat",
                        "--optical_backend", "gee", "--optical_level", "sr"])
    assert call["func"] == "main_historic"
    assert call["kwargs"]["start_date"] == "2023-02-06" and call["kwargs"]["end_date"] is None
    assert (call["kwargs"]["sensor"], call["kwargs"]["optical_backend"], call["kwargs"]["optical_level"]) == \
        ("landsat", "gee", "sr")
    assert call["kwargs"]["pairing_mode"] is None


def test_eq_list() -> None:
    call = single_call(["--eq_list", "events.json", "--sensor", "sentinel-2", "--optical_backend", "gee"])
    assert call["func"] == "main_historic"
    assert call["kwargs"]["eq_list_path"] == "events.json"
    assert "start_date" not in call["kwargs"]


@pytest.mark.parametrize("argv", [
    ["--forward", "--pairing", "coseismic", "--sensor", "sentinel-2"],          # optical forward is Phase 3
    ["--historic", "--dates", "2025-01-07", "--pairing", "coseismic", "--optical_backend", "gee"],  # optical flag, SAR
    ["--historic", "--dates", "2025-01-07"],                                       # SAR needs --pairing
    ["--historic", "--pairing", "coseismic"],                                      # missing --dates
    ["--historic", "--dates", "2025-01-01", "2025-01-02", "2025-01-03", "--pairing", "coseismic"],
    ["--historic", "--dates", "2025-01-07", "--pairing", "coseismic", "--send_email"],
    ["--forward", "--pairing", "coseismic", "--job_list"],
    ["--forward", "--pairing", "coseismic", "--dates", "2025-01-07"],
    ["--eq_list", "events.json", "--dates", "2025-01-07", "--pairing", "coseismic"],
])
def test_invalid_combinations_exit_with_error(argv: List[str]) -> None:
    with pytest.raises(SystemExit) as exit_info:
        api.run_cli(argv)
    assert exit_info.value.code == 1


@pytest.mark.parametrize("argv", [
    [],                                                                    # a run mode is required
    ["--historic", "--forward", "--pairing", "coseismic"],                 # run modes are exclusive
    ["--forward", "--pairing", "interferogram"],                           # invalid choice
])
def test_argparse_rejects(argv: List[str]) -> None:
    with pytest.raises(SystemExit) as exit_info:
        api.run_cli(argv)
    assert exit_info.value.code == 2


def test_entry_point_script_runs_from_scripts_dir() -> None:
    """Cron and GitHub Actions run `cd scripts && python coseis.py ...`."""
    import subprocess
    import sys

    scripts_dir = api.SRC_DIR.parent / "scripts"
    # Popen, not subprocess.run: run() is replaced by the topsApp fake
    with subprocess.Popen([sys.executable, "coseis.py", "--help"], cwd=scripts_dir,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as proc:
        out, err = proc.communicate(timeout=120)
    assert proc.returncode == 0, err
    assert "--forward" in out and "--process_only" in out
