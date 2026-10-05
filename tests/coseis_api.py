"""
The single place the tests import code under test from.

When a function moves between modules of the aria_coseis package, update its import here
and leave the tests unchanged.
"""

import sys
from pathlib import Path
from typing import Any, Callable, Dict, List

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC_DIR))

# Code under test. When a function moves to another module, update its import here only.
# `settings` is the module holding root_dir, TRACKING_DIR, LOCK_FILE and the recipients.
from aria_coseis import cli  # noqa: E402
from aria_coseis import config as settings  # noqa: E402
from aria_coseis.aoi import make_aoi  # noqa: E402
from aria_coseis.modes import main_forward, main_historic  # noqa: E402
from aria_coseis.notify import ascii_table_to_html  # noqa: E402
from aria_coseis.optical.jobs import make_optical_job_json  # noqa: E402
from aria_coseis.pipeline import process_earthquake  # noqa: E402
from aria_coseis.sar.pairing import generate_pairs, make_job_json  # noqa: E402
from aria_coseis.significance import check_significance  # noqa: E402
from aria_coseis.tracking import (  # noqa: E402
    add_to_tracker,
    check_tracker_for_updates,
    load_tracker,
)
from aria_coseis.usgs import parse_custom_eq_list  # noqa: E402
from aria_coseis.utils import convert_time, to_snake_case  # noqa: E402

__all__ = [
    "SRC_DIR",
    "add_to_tracker",
    "ascii_table_to_html",
    "check_significance",
    "check_tracker_for_updates",
    "cli",
    "convert_time",
    "generate_pairs",
    "load_tracker",
    "main_forward",
    "main_historic",
    "make_aoi",
    "make_job_json",
    "make_optical_job_json",
    "parse_custom_eq_list",
    "process_earthquake",
    "run_cli",
    "settings",
    "to_snake_case",
]


def run_cli(argv: List[str]) -> List[Dict[str, Any]]:
    """
    Run the command-line entry point with `argv`, recording calls to
    main_forward/main_historic instead of executing them.
    :return: list of {"func": name, "args": [...], "kwargs": {...}} in call order
    """
    calls: List[Dict[str, Any]] = []

    def recorder(name: str) -> Callable[..., None]:
        def record(*args: Any, **kwargs: Any) -> None:
            calls.append({"func": name, "args": list(args), "kwargs": kwargs})

        return record

    real = (cli.main_forward, cli.main_historic, sys.argv)
    cli.main_forward, cli.main_historic = recorder("main_forward"), recorder("main_historic")
    sys.argv = ["coseis.py", *argv]
    try:
        cli.main()
    finally:
        cli.main_forward, cli.main_historic, sys.argv = real
    return calls
