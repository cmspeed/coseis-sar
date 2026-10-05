"""
The single place the tests import code under test from.

Phase 2b moves functions out of scripts/coseis.py into the aria_coseis package.
When a function moves, update its import here and leave the tests unchanged.
"""
import ast
import sys
from pathlib import Path
from types import ModuleType
from typing import Any, Callable, Dict, List

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import coseis  # noqa: E402

# Module whose attributes hold runtime settings (root_dir, TRACKING_DIR, recipients).
settings: ModuleType = coseis

to_snake_case = coseis.to_snake_case
convert_time = coseis.convert_time
make_aoi = coseis.make_aoi
generate_pairs = coseis.generate_pairs
make_job_json = coseis.make_job_json
make_optical_job_json = coseis.make_optical_job_json
ascii_table_to_html = coseis.ascii_table_to_html
parse_custom_eq_list = coseis.parse_custom_eq_list
check_significance = coseis.check_significance
process_earthquake = coseis.process_earthquake
load_tracker = coseis.load_tracker
add_to_tracker = coseis.add_to_tracker
check_tracker_for_updates = coseis.check_tracker_for_updates
main_historic = coseis.main_historic
main_forward = coseis.main_forward


def run_cli(argv: List[str]) -> List[Dict[str, Any]]:
    """
    Run the coseis.py command-line entry point with `argv`, recording calls to
    main_forward/main_historic instead of executing them.
    :return: list of {"func": name, "args": [...], "kwargs": {...}} in call order
    """
    calls: List[Dict[str, Any]] = []

    def recorder(name: str) -> Callable[..., None]:
        def record(*args: Any, **kwargs: Any) -> None:
            calls.append({"func": name, "args": list(args), "kwargs": kwargs})
        return record

    # Execute only the `if __name__ == "__main__":` block, in a copy of the module namespace
    tree = ast.parse(Path(coseis.__file__).read_text())
    main_block = next(
        node for node in tree.body
        if isinstance(node, ast.If) and ast.unparse(node.test) == "__name__ == '__main__'"
    )
    code = compile(ast.Module(body=main_block.body, type_ignores=[]), coseis.__file__, "exec")
    namespace = dict(vars(coseis))
    namespace["main_forward"] = recorder("main_forward")
    namespace["main_historic"] = recorder("main_historic")

    old_argv = sys.argv
    sys.argv = ["coseis.py", *argv]
    try:
        exec(code, namespace)
    finally:
        sys.argv = old_argv
    return calls
