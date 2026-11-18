"""
Precise-orbit reruns (aria_coseis.precise): the queue of urgent-response products, the production
products made from them, and the email outbox the GitHub Action sends.

The queued product is the Ende D61 pair of tests/fixtures/production/ende_98ea2b5, placed as forward
processing places it. ASF's orbit service and the USGS significance check are faked.
"""

import json
import os
from pathlib import Path
from typing import Any, Dict, List

import pytest

import coseis_api as api
from conftest import FAKE_ORBITS, FAKE_PRODUCT, read_json

PRODUCTION = Path(__file__).parent / "fixtures" / "production"
EVENT_ID = "us6000tkt2"
EVENT_FOLDER = "2026/m_77_68_km_nnw_of_ende_indonesia_us6000tkt2"
URGENT = FAKE_PRODUCT.replace("S1-COSEIS_SAR-", "S1-COSEIS_SAR_UR-")
POST_SLCS = [
    "S1C_IW_SLC__1SDV_20260820T213528_20260820T213558_009084_01208A_0DB0",
    "S1C_IW_SLC__1SDV_20260820T213556_20260820T213631_009084_01208A_BC33",
]


def archive(*parts: str) -> Path:
    return Path(api.settings.ARCHIVE_DIR, *parts)


def queue_file() -> Path:
    return archive("MANIFESTS", "precise_queue", f"{URGENT}.json")


def staged_notes() -> List[Dict[str, Any]]:
    return [read_json(p) for p in sorted(Path(api.settings.OUTBOX_STAGING_DIR).glob("*.json"))]


def queue_product(workdir: Path, orbits: str = "RESORB") -> Dict[str, Any]:
    """Place and queue the Ende D61 urgent product as forward processing does."""
    processing_dir = workdir / "data" / "ende" / "DESCENDING061_20260802_20260820"
    product_dir = processing_dir / FAKE_PRODUCT
    product_dir.mkdir(parents=True)
    for ext in (".nc", ".png", ".json"):
        (product_dir / f"{FAKE_PRODUCT}{ext}").write_text(f"urgent {ext}")
    loc = {
        "reference_orbits": FAKE_ORBITS[orbits][:1],
        "secondary_orbits": FAKE_ORBITS[orbits][1:],
    }
    (processing_dir / "loc_data.json").write_text(json.dumps(loc))

    event = read_json(PRODUCTION / "ende_98ea2b5" / f"{EVENT_ID}.json")
    (job,) = read_json(
        PRODUCTION
        / "ende_98ea2b5"
        / "job_m_77_68_km_nnw_of_ende_indonesia_DESCENDING_61_partial.json"
    )
    job["job_parameters"]["granules"] = POST_SLCS
    api.precise.queue_urgent_product(
        str(product_dir / f"{FAKE_PRODUCT}.nc"),
        str(processing_dir),
        EVENT_ID,
        event,
        "DESCENDING_61",
        job,
        "2026-08-02 (Pre) - 2026-08-20 (Post)",
    )
    return read_json(queue_file())


@pytest.fixture
def services(monkeypatch: pytest.MonkeyPatch) -> Dict[str, Any]:
    """Fake ASF orbit service and USGS significance check; set "orbit" and "eligible"."""
    state: Dict[str, Any] = {"orbit": "POEORB", "eligible": True, "scenes": [], "events": []}

    def orbit_type(scene: str) -> Any:
        state["scenes"].append(scene)
        return state["orbit"]

    def meets_criteria(event_id: str) -> Any:
        state["events"].append(event_id)
        return state["eligible"]

    monkeypatch.setattr(api.precise, "orbit_type", orbit_type)
    monkeypatch.setattr(api.precise, "meets_criteria", meets_criteria)
    return state


def test_names() -> None:
    assert api.sar_products.urgent_name(FAKE_PRODUCT) == URGENT
    assert api.sar_products.production_name(URGENT) == FAKE_PRODUCT
    with pytest.raises(ValueError):
        api.sar_products.urgent_name(URGENT)
    with pytest.raises(ValueError):
        api.sar_products.production_name(FAKE_PRODUCT)
    assert api.sar_products.event_year(1786744701564) == "2026"


def test_queue_entry(workdir: Path, topsapp: Dict[str, Any]) -> None:
    entry = queue_product(workdir)
    assert entry["status"] == "pending"
    assert entry["urgent_product"] == f"URGENT_RESPONSE/{EVENT_FOLDER}/{URGENT}"
    assert entry["production_product"] == f"HISTORIC_EVENTS/{EVENT_FOLDER}/{FAKE_PRODUCT}"
    assert entry["urgent_orbits"] == sorted(FAKE_ORBITS["RESORB"])
    assert sorted(p.name for p in archive(entry["urgent_product"]).iterdir()) == [
        f"{URGENT}.nc",
        f"{URGENT}.png",
    ]
    assert archive(entry["urgent_product"], f"{URGENT}.nc").read_text() == "urgent .nc"


def test_waits_for_precise_orbits(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    before = queue_product(workdir)
    services["orbit"] = "RESORB"
    api.precise.check_precise_queue()
    assert read_json(queue_file()) == before
    assert sorted(services["scenes"]) == sorted(
        before["job"]["job_parameters"]["granules"]
        + before["job"]["job_parameters"]["secondary_granules"]
    )
    assert services["events"] == [] and topsapp["calls"] == [] and staged_notes() == []


def test_gives_up_without_precise_orbits(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    entry = queue_product(workdir)
    entry["log"][0] = "2026-08-21T00:00:00Z queued"
    queue_file().write_text(json.dumps(entry))
    services["orbit"] = None
    api.precise.check_precise_queue()
    entry = read_json(queue_file())
    assert entry["status"] == "failed"
    assert "no precise orbits after" in entry["log"][-1]
    (note,) = staged_notes()
    assert note["kind"] == "failed" and note["track"] == "DESCENDING_61"
    assert topsapp["calls"] == []


@pytest.mark.parametrize("eligible, status", [(False, "not_eligible"), (None, "pending")])
def test_significance_decides(
    eligible: Any,
    status: str,
    workdir: Path,
    topsapp: Dict[str, Any],
    services: Dict[str, Any],
) -> None:
    """Not significant: never rerun. Undecidable (USGS down): try again next run."""
    queue_product(workdir)
    services["eligible"] = eligible
    api.precise.check_precise_queue()
    assert read_json(queue_file())["status"] == status
    assert services["events"] == [EVENT_ID]
    assert topsapp["calls"] == [] and staged_notes() == []


def test_production_product(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    queue_product(workdir)
    topsapp.update(outcome="success", orbits="POEORB")
    api.precise.check_precise_queue()

    # The same job, rerun on the precise runner's cores in its own scratch folder
    (call,) = topsapp["calls"]
    cmd = call["cmd"]
    assert cmd[cmd.index("taskset") + 2] == "17-33"
    assert cmd[cmd.index("--reference-scenes") + 1] == " ".join(POST_SLCS)
    work_dir = Path(call["cwd"])
    assert work_dir == Path(
        api.settings.root_dir,
        "m_77_68_km_nnw_of_ende_indonesia",
        "DESCENDING061",
        "precise",
        "DESCENDING061_20260802_20260820",
    )
    assert not (work_dir / "fine_interferogram").exists()

    # Placed in HISTORIC_EVENTS under the production name (.nc and .png only)
    production = archive("HISTORIC_EVENTS", EVENT_FOLDER, FAKE_PRODUCT)
    assert sorted(p.name for p in production.iterdir()) == [
        f"{FAKE_PRODUCT}.nc",
        f"{FAKE_PRODUCT}.png",
    ]
    assert (production / f"{FAKE_PRODUCT}.nc").read_text() == "fake .nc"
    entry = read_json(queue_file())
    assert entry["status"] == "done"
    (note,) = staged_notes()
    assert note["kind"] == "production"
    assert note["location"] == str(production)
    assert note["dates"] == "2026-08-02 (Pre) - 2026-08-20 (Post)"

    # Nothing left to do on the next run
    topsapp["calls"].clear()
    api.precise.check_precise_queue()
    assert topsapp["calls"] == [] and len(staged_notes()) == 1


@pytest.mark.parametrize(
    "change, message",
    [
        ({"product": FAKE_PRODUCT.replace("5d2a", "0000")}, "is not the expected"),
        ({"orbits": "RESORB"}, "did not use precise orbits"),
        ({"outcome": "no_product"}, "wrote no product"),
        ({"outcome": "nonzero_exit"}, "returned non-zero exit status"),
    ],
)
def test_rerun_problems_are_flagged(
    change: Dict[str, str],
    message: str,
    workdir: Path,
    topsapp: Dict[str, Any],
    services: Dict[str, Any],
) -> None:
    queue_product(workdir)
    topsapp.update(outcome="success", orbits="POEORB")
    topsapp.update(change)
    api.precise.check_precise_queue()
    entry = read_json(queue_file())
    assert entry["status"] == "failed"
    assert message in entry["log"][-1]
    assert not archive("HISTORIC_EVENTS").exists()
    (note,) = staged_notes()
    assert note["kind"] == "failed" and message in note["message"]


def test_urgent_product_with_precise_orbits_is_copied(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    """Processed late enough to have precise orbits: the urgent product is the production one."""
    queue_product(workdir, orbits="POEORB")
    api.precise.check_precise_queue()
    production = archive("HISTORIC_EVENTS", EVENT_FOLDER, FAKE_PRODUCT)
    assert (production / f"{FAKE_PRODUCT}.nc").read_text() == "urgent .nc"
    assert read_json(queue_file())["status"] == "done"
    assert topsapp["calls"] == [] and services["scenes"] == []
    assert services["events"] == [EVENT_ID]


def test_existing_production_product_is_left_alone(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    queue_product(workdir)
    production = archive("HISTORIC_EVENTS", EVENT_FOLDER, FAKE_PRODUCT)
    production.mkdir(parents=True)
    api.precise.check_precise_queue()
    assert read_json(queue_file())["status"] == "done"
    assert list(production.iterdir()) == [] and topsapp["calls"] == []


def test_lock_blocks_a_second_run(
    workdir: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    topsapp: Dict[str, Any],
    services: Dict[str, Any],
) -> None:
    queue_product(workdir)
    lock = tmp_path / "precise.lock"
    monkeypatch.setattr(api.settings, "PRECISE_LOCK_FILE", str(lock))
    lock.write_text(f"{os.getpid()} now\n")  # held by a live process
    api.main_precise()
    assert services["scenes"] == [] and lock.exists()

    lock.unlink()
    services["orbit"] = "RESORB"
    api.main_precise()
    assert services["scenes"] and not lock.exists()


def test_outbox_is_emailed_once(
    workdir: Path,
    topsapp: Dict[str, Any],
    services: Dict[str, Any],
    sent_emails: List[Dict[str, Any]],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    queue_product(workdir)
    topsapp.update(outcome="success", orbits="POEORB")
    api.precise.check_precise_queue()

    # run_coseis_precise.sh moves the staged notes into active_jobs/outbox and pushes
    outbox = Path(api.settings.OUTBOX_DIR)
    outbox.mkdir(parents=True)
    for path in Path(api.settings.OUTBOX_STAGING_DIR).glob("*.json"):
        path.rename(outbox / path.name)

    # Without email credentials the notes stay for the next run
    monkeypatch.delenv("GMAIL_USER")
    api.precise.send_outbox()
    assert sent_emails == [] and len(list(outbox.glob("*.json"))) == 1
    monkeypatch.setenv("GMAIL_USER", "sender@example.com")

    api.precise.send_outbox()
    (email,) = sent_emails
    assert email["subject"] == "PRODUCTION PRODUCTS (PRECISE ORBITS): 1 placed"
    assert email["bcc"] == ["secondary@example.com"]
    assert FAKE_PRODUCT in email["contents"][0]
    assert list(outbox.glob("*.json")) == []

    api.precise.send_outbox()
    assert len(sent_emails) == 1


def test_orbit_service_errors_are_retried(
    workdir: Path, topsapp: Dict[str, Any], monkeypatch: pytest.MonkeyPatch
) -> None:
    import requests

    before = queue_product(workdir)

    def unreachable(scene: str) -> None:
        raise requests.ConnectionError("no route to host")

    monkeypatch.setattr(api.precise, "orbit_type", unreachable)
    api.precise.check_precise_queue()
    assert read_json(queue_file()) == before and topsapp["calls"] == []
