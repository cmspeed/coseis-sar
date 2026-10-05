"""
The single place the tests import code under test from.

As functions move between modules of the aria_coseis package, update the imports here
and leave the tests unchanged.
"""
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SRC_DIR))

from aria_coseis import cli  # noqa: E402

# --- code under test (generated; keep sorted by module) ---
from aria_coseis import config as settings  # noqa: E402  (holds root_dir, TRACKING_DIR, recipients)
from aria_coseis.legacy import (  # noqa: E402
    add_to_tracker,
    check_significance,
    check_tracker_for_updates,
    generate_pairs,
    load_tracker,
    main_forward,
    main_historic,
    make_aoi,
    make_job_json,
    make_optical_job_json,
    process_earthquake,
)
from aria_coseis.notify import (  # noqa: E402
    ascii_table_to_html,
)
from aria_coseis.usgs import (  # noqa: E402
    parse_custom_eq_list,
)
from aria_coseis.utils import (  # noqa: E402
    convert_time,
    to_snake_case,
)
# --- end code under test ---


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
