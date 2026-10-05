"""
Shared fixtures for the coseis characterization tests.

External boundaries are faked so tests never send email, run topsApp or touch the network:
- HTTP is replayed from vcrpy cassettes (tests/cassettes/); unrecorded requests fail.
- yagmail.SMTP, subprocess.run and the next_pass module are replaced with fakes.
"""
import json
import os
import subprocess
import sys
import types
from pathlib import Path
from typing import Any, Callable, Dict, List

import pytest

# Recipients and credentials are read when coseis is imported: set test values first.
os.environ["COSEIS_PRIMARY_RECIPIENTS"] = "primary@example.com"
os.environ["COSEIS_SECONDARY_RECIPIENTS"] = "secondary@example.com"
os.environ["GMAIL_USER"] = "sender@example.com"
os.environ["GMAIL_APP_PSWD"] = "not-a-password"

import coseis_api as api  # noqa: E402

TESTS_DIR = Path(__file__).resolve().parent
GOLDEN_DIR = TESTS_DIR / "golden"
SHARED_CASSETTES = TESTS_DIR / "cassettes" / "shared"
COASTLINE_CASSETTE = str(SHARED_CASSETTES / "coastline.yaml")
UPDATE_GOLDEN = os.environ.get("UPDATE_GOLDEN") == "1"


@pytest.fixture(scope="module")
def vcr_config() -> Dict[str, Any]:
    # Repeats allowed: one cassette per module, and e.g. the coastline is fetched on every call
    return {"decode_compressed_response": True, "allow_playback_repeats": True}


@pytest.fixture
def default_cassette_name(request: pytest.FixtureRequest) -> str:
    """One cassette per test module (tests/cassettes/<module>/<module>.yaml)."""
    return request.module.__name__.rsplit(".", 1)[-1]


@pytest.fixture
def http_log(monkeypatch: pytest.MonkeyPatch) -> List[str]:
    """URLs of every requests.get call, logged before vcrpy replays (or blocks) them."""
    import requests

    urls: List[str] = []
    real_get = requests.get

    def logging_get(url: str, params: Any = None, **kwargs: Any) -> Any:
        urls.append(requests.Request("GET", url, params=params).prepare().url)
        return real_get(url, params=params, **kwargs)

    monkeypatch.setattr(requests, "get", logging_get)
    return urls


@pytest.fixture
def workdir(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Run in an empty <tmp>/scripts directory, mirroring production's cwd of scripts/."""
    scripts_dir = tmp_path / "scripts"
    scripts_dir.mkdir()
    monkeypatch.chdir(scripts_dir)
    monkeypatch.setattr(api.settings, "root_dir", str(scripts_dir / "data"))
    return scripts_dir


@pytest.fixture(autouse=True)
def sent_emails(monkeypatch: pytest.MonkeyPatch) -> List[Dict[str, Any]]:
    """Replace yagmail.SMTP; every send is recorded here instead of being delivered."""
    import yagmail

    sent: List[Dict[str, Any]] = []

    class FakeSMTP:
        def __init__(self, user: str, password: str) -> None:
            self.user = user

        def send(self, **kwargs: Any) -> None:
            sent.append(kwargs)

    monkeypatch.setattr(yagmail, "SMTP", FakeSMTP)
    return sent


@pytest.fixture(autouse=True)
def topsapp(monkeypatch: pytest.MonkeyPatch) -> Dict[str, Any]:
    """
    Intercept topsApp runs (subprocess.run commands containing "isce2_topsapp"); other commands,
    e.g. matplotlib's fc-list when building its font cache, run normally.
    Set topsapp["outcome"] to choose the behavior: "success" writes a fake .nc product,
    "nonzero_exit" raises CalledProcessError, "missing_conda" raises FileNotFoundError.
    Calls are recorded in topsapp["calls"].
    """
    state: Dict[str, Any] = {"outcome": "unexpected", "calls": []}
    real_run = subprocess.run

    def fake_run(cmd: List[str], cwd: str = None, **kwargs: Any) -> subprocess.CompletedProcess:
        if "isce2_topsapp" not in cmd:
            return real_run(cmd, cwd=cwd, **kwargs)
        state["calls"].append({"cmd": cmd, "cwd": cwd})
        outcome = state["outcome"]
        if outcome == "success":
            Path(cwd, "S1-GUNW-test-product.nc").write_text("fake")
            Path(cwd, "S1A_IW_SLC__fake.zip").write_text("fake")
            Path(cwd, "fine_interferogram").mkdir()
            return subprocess.CompletedProcess(cmd, 0, stdout="ok", stderr="")
        if outcome == "nonzero_exit":
            raise subprocess.CalledProcessError(1, cmd, output="", stderr="topsApp failed")
        if outcome == "missing_conda":
            raise FileNotFoundError("conda")
        raise AssertionError(f"Unexpected subprocess.run call: {cmd}")

    monkeypatch.setattr(subprocess, "run", fake_run)
    return state


@pytest.fixture
def fake_next_pass(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Dict[str, Any]:
    """Install a fake next_pass module that returns fixed overpass info and writes a stub map."""
    calls: Dict[str, Any] = {"bboxes": []}

    def find_next_overpass(args: Any, timestamp_dir: Path) -> Dict[str, Any]:
        calls["bboxes"].append(list(args.bbox))
        return {
            "sentinel-1": {"next_collect_info": "+-----+\n| S1 |\n+-----+\n| next S1 pass |\n+-----+"},
            "nisar": {"next_collect_info": "No NISAR passes."},
        }

    def make_overpasses_map(*args: Any) -> None:
        timestamp_dir = Path(args[-1])
        (timestamp_dir / "satellite_overpasses_map.html").write_text("<html>map</html>")

    package_dir = tmp_path / "fake_next_pass"
    package_dir.mkdir()
    module = types.ModuleType("next_pass")
    module.__file__ = str(package_dir / "__init__.py")
    module.find_next_overpass = find_next_overpass
    plot_maps = types.ModuleType("next_pass.plot_maps")
    plot_maps.make_overpasses_map = make_overpasses_map
    module.plot_maps = plot_maps
    monkeypatch.setitem(sys.modules, "next_pass", module)
    monkeypatch.setitem(sys.modules, "next_pass.plot_maps", plot_maps)
    return calls


@pytest.fixture
def golden() -> Callable[[str, Any], None]:
    """
    Compare JSON-serializable data with tests/golden/<name>.json.
    Run with UPDATE_GOLDEN=1 to (re)write the golden files, then review the diff.
    """
    def check(name: str, data: Any) -> None:
        path = GOLDEN_DIR / f"{name}.json"
        rendered = json.dumps(data, indent=2, sort_keys=True, default=str) + "\n"
        if UPDATE_GOLDEN or not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(rendered)
            if not UPDATE_GOLDEN:
                pytest.fail(f"Golden file {path} did not exist; it has been written. Review it and re-run.")
            return
        assert json.loads(rendered) == json.loads(path.read_text()), f"Output differs from {path}"

    return check


def read_json(path: Path) -> Any:
    with open(path) as f:
        return json.load(f)
