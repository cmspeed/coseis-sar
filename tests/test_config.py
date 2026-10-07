"""Settings defaults and environment overrides (aria_coseis.config is read at import)."""

import importlib
import os
from pathlib import Path
from typing import Dict, Iterator

import pytest

import coseis_api as api

OVERRIDES = (
    "COSEIS_DATA_DIR",
    "COSEIS_TRACKING_DIR",
    "COSEIS_PARTIALS_DIR",
    "COSEIS_MAPS_DIR",
    "COSEIS_OUTPUT_DIR",
    "COSEIS_LOCK_FILE",
)


@pytest.fixture
def reload_config() -> Iterator:
    """Reload config under a given environment, then restore the defaults."""
    saved = {name: os.environ.pop(name, None) for name in OVERRIDES}

    def reload(env: Dict[str, str]):
        os.environ.update(env)
        return importlib.reload(api.settings)

    yield reload
    for name in OVERRIDES:
        os.environ.pop(name, None)
        if saved[name] is not None:
            os.environ[name] = saved[name]
    importlib.reload(api.settings)


def test_defaults(reload_config, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    config = reload_config({})
    assert config.root_dir == os.path.join(str(tmp_path), "data")
    assert config.TRACKING_DIR == "active_jobs"
    assert config.PARTIALS_DIR == os.path.join("active_jobs", "partials")
    assert config.MAPS_DIR == os.path.join("docs", "maps")
    assert config.OUTPUT_DIR == "outputs"
    assert config.LOCK_FILE == "/tmp/coseis_processing.lock"


def test_environment_overrides(reload_config, tmp_path: Path) -> None:
    config = reload_config(
        {
            "COSEIS_DATA_DIR": str(tmp_path / "shadow_data"),
            "COSEIS_TRACKING_DIR": str(tmp_path / "shadow_jobs"),
            "COSEIS_LOCK_FILE": str(tmp_path / "shadow.lock"),
            "COSEIS_MAPS_DIR": str(tmp_path / "shadow_maps"),
            "COSEIS_OUTPUT_DIR": str(tmp_path / "shadow_outputs"),
        }
    )
    assert config.root_dir == str(tmp_path / "shadow_data")
    assert config.TRACKING_DIR == str(tmp_path / "shadow_jobs")
    assert config.PARTIALS_DIR == str(
        tmp_path / "shadow_jobs" / "partials"
    )  # follows TRACKING_DIR
    assert config.LOCK_FILE == str(tmp_path / "shadow.lock")
    assert config.MAPS_DIR == str(tmp_path / "shadow_maps")
    assert config.output_path("x.json") == str(tmp_path / "shadow_outputs" / "x.json")


def test_forward_run_uses_configured_lock_file(
    workdir: Path, monkeypatch: pytest.MonkeyPatch, http_log: list
) -> None:
    lock = workdir / "shadow.lock"
    lock.write_text("running")
    monkeypatch.setattr(api.settings, "LOCK_FILE", str(lock))
    api.main_forward(pairing_mode="coseismic", send_email_flag=True)
    assert http_log == [], "a held lock must stop the run before any work"
