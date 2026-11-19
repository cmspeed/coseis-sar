"""
Precise-orbit reruns (aria_coseis.precise): the queue of urgent-response products, shared pairs, the
production products made from them with historic mode's selection, and the email outbox.

The queued product is the Ende D61 pair of tests/fixtures/production/ende_98ea2b5, placed as forward
processing places it. USGS (conftest's usgs_events), ASF's orbit service and historic mode's pair
selection are faked.
"""

import copy
import json
import os
from pathlib import Path
from typing import Any, Dict, List

import pytest

import coseis_api as api
from conftest import FAKE_ORBITS, FAKE_PRODUCT, read_json

PRODUCTION = Path(__file__).parent / "fixtures" / "production"
MAIN = "us6000tkt2"  # M7.8 Ende, 2026-08-14 21:58 UTC
AFTERSHOCK = "us6000tkuc"  # M6.1, 30 min later
SMALL = "us6000tku7"  # M5.6, not significant for historic mode
URGENT = FAKE_PRODUCT.replace("S1-COSEIS_SAR-", "S1-COSEIS_SAR_UR-")
POST_SLCS = [
    "S1C_IW_SLC__1SDV_20260820T213528_20260820T213558_009084_01208A_0DB0",
    "S1C_IW_SLC__1SDV_20260820T213556_20260820T213631_009084_01208A_BC33",
]
TIME = {MAIN: 1786744701564, AFTERSHOCK: 1786746517000, SMALL: 1786745000000}
TITLE = {
    MAIN: "m_77_68_km_nnw_of_ende_indonesia",
    AFTERSHOCK: "m_61_58_km_n_of_ende_indonesia",
    SMALL: "m_56_61_km_nnw_of_ende_indonesia",
}
USGS = {
    MAIN: {"title": "M 7.8 - 66 km NNW of Ende, Indonesia", "mag": 7.8, "historic": True},
    AFTERSHOCK: {"title": "M 6.1 - 62 km N of Ende, Indonesia", "mag": 6.1, "historic": True},
    SMALL: {"title": "M 5.6 - 49 km NNW of Ende, Indonesia", "mag": 5.6, "historic": False},
}
MAIN_FOLDER = "2026/m_78_66_km_nnw_of_ende_indonesia_us6000tkt2"  # current USGS title


def archive(*parts: str) -> Path:
    return Path(api.settings.ARCHIVE_DIR, *parts)


def queue_file() -> Path:
    return archive("MANIFESTS", "precise_queue", f"{URGENT}.json")


def staged_notes() -> List[Dict[str, Any]]:
    return [read_json(p) for p in sorted(Path(api.settings.OUTBOX_STAGING_DIR).glob("*.json"))]


def urgent_job() -> Dict[str, Any]:
    (job,) = read_json(
        PRODUCTION
        / "ende_98ea2b5"
        / "job_m_77_68_km_nnw_of_ende_indonesia_DESCENDING_61_partial.json"
    )
    job["job_parameters"]["granules"] = POST_SLCS
    return job


def queue_product(workdir: Path, event_id: str = MAIN, orbits: str = "RESORB") -> Dict[str, Any]:
    """Place and queue the Ende D61 urgent product for event_id, as forward processing does."""
    processing_dir = workdir / "data" / TITLE[event_id] / "DESCENDING061_20260802_20260820"
    product_dir = processing_dir / FAKE_PRODUCT
    product_dir.mkdir(parents=True, exist_ok=True)
    for ext in (".nc", ".png", ".json"):
        (product_dir / f"{FAKE_PRODUCT}{ext}").write_text(f"urgent {ext}")
    loc = {
        "reference_orbits": FAKE_ORBITS[orbits][:1],
        "secondary_orbits": FAKE_ORBITS[orbits][1:],
    }
    (processing_dir / "loc_data.json").write_text(json.dumps(loc))
    api.precise.queue_urgent_product(
        str(product_dir / f"{FAKE_PRODUCT}.nc"),
        str(processing_dir),
        event_id,
        {"title": TITLE[event_id], "time": TIME[event_id]},
        "DESCENDING_61",
        urgent_job(),
        "2026-08-02 (Pre) - 2026-08-20 (Post)",
    )
    return read_json(queue_file())


@pytest.fixture
def services(monkeypatch: pytest.MonkeyPatch, usgs_events: Dict[str, Any]) -> Dict[str, Any]:
    """
    Fake ASF orbit service and historic-mode selection; USGS knows the three Ende events. Set
    "orbit" (what ASF serves for every scene) and "job" (historic mode's pair; default: the urgent
    product's own pair).
    """
    for eid, info in USGS.items():
        usgs_events[eid] = {"id": eid, "time": TIME[eid], **info}
    state: Dict[str, Any] = {"orbit": "POEORB", "job": "urgent", "scenes": [], "owners": []}

    def orbit_type(scene: str) -> Any:
        state["scenes"].append(scene)
        return state["orbit"]

    def historic_job(owner: Dict[str, Any], track_key: str) -> Any:
        state["owners"].append((owner["id"], track_key))
        if state["job"] == "urgent":
            job = urgent_job()
            job["name"] = "m_78_66_km_nnw_of_ende_indonesia-D061"
            return job
        if isinstance(state["job"], dict) and "job_parameters" not in state["job"]:
            return copy.deepcopy(state["job"][owner["id"]])  # a pair per owner
        return copy.deepcopy(state["job"])

    monkeypatch.setattr(api.precise, "orbit_type", orbit_type)
    monkeypatch.setattr(api.precise, "historic_job", historic_job)
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
    assert entry["urgent_product"] == (
        f"URGENT_RESPONSE/2026/m_77_68_km_nnw_of_ende_indonesia_us6000tkt2/{URGENT}"
    )
    assert [e["event_id"] for e in entry["events"]] == [MAIN]
    assert entry["urgent_orbits"] == sorted(FAKE_ORBITS["RESORB"])
    assert sorted(p.name for p in archive(entry["urgent_product"]).iterdir()) == [
        f"{URGENT}.nc",
        f"{URGENT}.png",
    ]
    assert archive(entry["urgent_product"], f"{URGENT}.nc").read_text() == "urgent .nc"


def test_shared_pair_is_placed_once_under_the_largest_event(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    """Forward mode processes the same pair for each event of a sequence: one urgent copy."""
    entry = queue_product(workdir, SMALL)
    assert TITLE[SMALL] in entry["urgent_product"]
    entry = queue_product(workdir, AFTERSHOCK)  # larger: the copy moves to its folder
    assert TITLE[AFTERSHOCK] in entry["urgent_product"]
    assert not archive("URGENT_RESPONSE", "2026", f"{TITLE[SMALL]}_{SMALL}").exists()
    entry = queue_product(workdir, MAIN)
    entry = queue_product(workdir, AFTERSHOCK)  # processed again: nothing changes
    assert [e["event_id"] for e in entry["events"]] == [SMALL, AFTERSHOCK, MAIN]
    assert entry["urgent_product"] == (
        f"URGENT_RESPONSE/2026/m_77_68_km_nnw_of_ende_indonesia_us6000tkt2/{URGENT}"
    )
    assert len(topsapp["placed"]) == 2  # .nc and .png, once
    assert [p.name for p in archive("URGENT_RESPONSE", "2026").iterdir()] == [
        "m_77_68_km_nnw_of_ende_indonesia_us6000tkt2"
    ]
    assert len(list(archive("MANIFESTS", "precise_queue").glob("*.json"))) == 1


def test_waits_for_precise_orbits(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    queue_product(workdir)
    services["orbit"] = "RESORB"
    api.precise.check_precise_queue()
    entry = read_json(queue_file())
    assert entry["status"] == "pending"
    assert entry["production"]["owner"] == MAIN  # the selection is kept for the next run
    p = entry["production"]["job"]["job_parameters"]
    assert sorted(services["scenes"]) == sorted(p["granules"] + p["secondary_granules"])
    assert topsapp["calls"] == [] and staged_notes() == []

    services["scenes"].clear()
    api.precise.check_precise_queue()
    assert services["owners"] == [(MAIN, "DESCENDING_61")]  # selected once
    assert services["scenes"]


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
    assert "no precise orbits for 4 of 4 SLCs after" in entry["log"][-1]
    (note,) = staged_notes()
    assert note["kind"] == "failed" and note["track"] == "DESCENDING_61"
    assert topsapp["calls"] == []


def test_waits_for_a_historic_pair(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    queue_product(workdir)
    services["job"] = None
    api.precise.check_precise_queue()
    entry = read_json(queue_file())
    assert entry["status"] == "pending" and "job" not in entry["production"]
    assert services["scenes"] == [] and topsapp["calls"] == []


@pytest.mark.parametrize("unreachable, status", [(False, "not_eligible"), (True, "pending")])
def test_significance_decides(
    unreachable: bool,
    status: str,
    workdir: Path,
    topsapp: Dict[str, Any],
    services: Dict[str, Any],
    usgs_events: Dict[str, Any],
) -> None:
    """Only small events used it: never a production product. USGS down: try again next run."""
    queue_product(workdir, SMALL)
    if unreachable:
        usgs_events[SMALL] = None
    usgs_events["calls"].clear()
    api.precise.check_precise_queue()
    assert read_json(queue_file())["status"] == status
    assert usgs_events["calls"] == [SMALL]
    assert services["owners"] == [] and topsapp["calls"] == [] and staged_notes() == []


def test_production_product(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    """Owner: the largest significant event; placed in its folder (current USGS title)."""
    queue_product(workdir, SMALL)
    queue_product(workdir, AFTERSHOCK)
    queue_product(workdir, MAIN)
    topsapp.update(outcome="success", orbits="POEORB")
    api.precise.check_precise_queue()
    assert services["owners"] == [(MAIN, "DESCENDING_61")]

    # Historic mode's job, on the precise runner's cores, in its own scratch folder
    (call,) = topsapp["calls"]
    cmd = call["cmd"]
    assert cmd[cmd.index("taskset") + 2] == "17-33"
    assert cmd[cmd.index("--reference-scenes") + 1] == " ".join(POST_SLCS)
    work_dir = Path(call["cwd"])
    assert work_dir == Path(
        api.settings.root_dir,
        "m_78_66_km_nnw_of_ende_indonesia",
        "D061",
        "precise",
        "D061_20260802_20260820",
    )
    assert not (work_dir / "fine_interferogram").exists()

    # In HISTORIC_EVENTS, in the owner's folder, .nc and .png only
    production = archive("HISTORIC_EVENTS", MAIN_FOLDER, FAKE_PRODUCT)
    assert sorted(p.name for p in production.iterdir()) == [
        f"{FAKE_PRODUCT}.nc",
        f"{FAKE_PRODUCT}.png",
    ]
    assert (production / f"{FAKE_PRODUCT}.nc").read_text() == "fake .nc"
    entry = read_json(queue_file())
    assert entry["status"] == "done"
    assert entry["production"]["product"] == f"HISTORIC_EVENTS/{MAIN_FOLDER}/{FAKE_PRODUCT}"
    (note,) = staged_notes()
    assert note["kind"] == "production"
    assert note["location"] == str(production)
    assert note["title"] == "m_78_66_km_nnw_of_ende_indonesia"
    assert note["dates"] == "20260802 (Pre) - 20260820 (Post)"

    # Nothing left to do on the next run
    topsapp["calls"].clear()
    api.precise.check_precise_queue()
    assert topsapp["calls"] == [] and len(staged_notes()) == 1


def test_production_uses_historic_selection(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    """Historic mode picks another pre-event date than forward mode did: that pair is run."""
    queue_product(workdir)
    job = urgent_job()
    job["name"] = "m_78_66_km_nnw_of_ende_indonesia-D061"
    job["job_parameters"]["secondary_granules"] = [
        "S1D_IW_SLC__1SDV_20260814T213613_20260814T213641_004127_007B2C_A1B2"
    ]
    services["job"] = job
    other = FAKE_PRODUCT.replace("5d2a", "77aa")
    topsapp.update(outcome="success", orbits="POEORB", product=other)
    api.precise.check_precise_queue()
    (call,) = topsapp["calls"]
    cmd = call["cmd"]
    assert (
        cmd[cmd.index("--secondary-scenes") + 1] == job["job_parameters"]["secondary_granules"][0]
    )
    assert call["cwd"].endswith("D061_20260814_20260820")
    assert archive("HISTORIC_EVENTS", MAIN_FOLDER, other).exists()
    assert read_json(queue_file())["status"] == "done"


def test_a_production_pair_runs_once(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    """Two urgent products (different frames) of the same owner and track: one production run."""
    queue_product(workdir, AFTERSHOCK)
    other = FAKE_PRODUCT.replace("5d2a", "beef")
    processing_dir = workdir / "data" / "other"
    (processing_dir / other).mkdir(parents=True)
    (processing_dir / other / f"{other}.nc").write_text("urgent")
    job = urgent_job()
    job["job_parameters"]["granules"] = POST_SLCS[:1]
    api.precise.queue_urgent_product(
        str(processing_dir / other / f"{other}.nc"),
        str(processing_dir),
        MAIN,
        {"title": TITLE[MAIN], "time": TIME[MAIN]},
        "DESCENDING_61",
        job,
        "2026-08-02 (Pre) - 2026-08-20 (Post)",
    )
    topsapp.update(outcome="success", orbits="POEORB")
    api.precise.check_precise_queue()
    assert len(topsapp["calls"]) == 1
    entries = [e for _, e in api.precise.load_queue()]
    assert [e["status"] for e in entries] == ["done", "done"]
    assert {e["production"]["product"] for e in entries} == {
        f"HISTORIC_EVENTS/{MAIN_FOLDER}/{FAKE_PRODUCT}"
    }
    assert len(staged_notes()) == 1


def test_smaller_event_keeps_its_own_pair(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    """
    The aftershock's urgent window contains the main shock, but the pair historic mode selects for
    the aftershock starts after the main shock: each event gets its own production product (as for
    Vanuatu's M6.1 and M6.0 on A154, 2026-10/11).
    """
    queue_product(workdir, MAIN)
    other = FAKE_PRODUCT.replace("5d2a", "beef")
    processing_dir = workdir / "data" / "other"
    (processing_dir / other).mkdir(parents=True)
    (processing_dir / other / f"{other}.nc").write_text("urgent")
    job = urgent_job()
    job["job_parameters"]["granules"] = POST_SLCS[:1]
    api.precise.queue_urgent_product(
        str(processing_dir / other / f"{other}.nc"),
        str(processing_dir),
        AFTERSHOCK,
        {"title": TITLE[AFTERSHOCK], "time": TIME[AFTERSHOCK]},
        "DESCENDING_61",
        job,
        "2026-08-02 (Pre) - 2026-08-20 (Post)",
    )
    between = "S1D_IW_SLC__1SDV_20260814T221500_20260814T221530_004127_007B2C_A1B2"  # 22:15
    main_job, after_job = urgent_job(), urgent_job()
    main_job["name"] = "m_78_66_km_nnw_of_ende_indonesia-D061"
    main_job["job_parameters"]["granules"] = [between]  # 08-02 -> 08-14 22:15: main shock only
    after_job["name"] = "m_61_62_km_n_of_ende_indonesia-D061"
    after_job["job_parameters"]["secondary_granules"] = [between]  # 22:15 -> 08-20: aftershock
    services["job"] = {MAIN: main_job, AFTERSHOCK: after_job}
    topsapp.update(outcome="success", orbits="POEORB")
    api.precise.check_precise_queue()
    assert len(topsapp["calls"]) == 2
    entries = {e["events"][0]["event_id"]: e for _, e in api.precise.load_queue()}
    assert [entries[e]["status"] for e in (MAIN, AFTERSHOCK)] == ["done", "done"]
    assert entries[AFTERSHOCK]["production"]["owner"] == AFTERSHOCK
    assert entries[AFTERSHOCK]["production"]["product"] == (
        f"HISTORIC_EVENTS/2026/m_61_62_km_n_of_ende_indonesia_{AFTERSHOCK}/{FAKE_PRODUCT}"
    )
    assert (
        entries[MAIN]["production"]["product"] == f"HISTORIC_EVENTS/{MAIN_FOLDER}/{FAKE_PRODUCT}"
    )


@pytest.mark.parametrize(
    "change, message",
    [
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
    """Processed late enough for precise orbits, and the pair historic mode selects."""
    queue_product(workdir, orbits="POEORB")
    api.precise.check_precise_queue()
    production = archive("HISTORIC_EVENTS", MAIN_FOLDER, FAKE_PRODUCT)
    assert (production / f"{FAKE_PRODUCT}.nc").read_text() == "urgent .nc"
    assert read_json(queue_file())["status"] == "done"
    assert topsapp["calls"] == [] and services["scenes"] == []


def test_existing_production_product_is_never_overwritten(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    queue_product(workdir)
    production = archive("HISTORIC_EVENTS", MAIN_FOLDER, FAKE_PRODUCT)
    production.mkdir(parents=True)
    topsapp.update(outcome="success", orbits="POEORB")
    api.precise.check_precise_queue()
    entry = read_json(queue_file())
    assert entry["status"] == "failed" and "already exists" in entry["log"][-1]
    assert list(production.iterdir()) == []


def test_events_added_during_a_run_are_kept(
    workdir: Path,
    topsapp: Dict[str, Any],
    services: Dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The urgent runner adds an event while the precise runner works on the entry."""
    queue_product(workdir, AFTERSHOCK)
    real = api.precise.process_entry

    def slow_process(entry: Dict[str, Any], others: List[Dict[str, Any]]) -> None:
        real(entry, others)
        queue_product(workdir, MAIN)  # meanwhile, in the urgent runner

    monkeypatch.setattr(api.precise, "process_entry", slow_process)
    services["orbit"] = "RESORB"
    api.precise.check_precise_queue()
    entry = read_json(queue_file())
    assert [e["event_id"] for e in entry["events"]] == [AFTERSHOCK, MAIN]
    assert entry["status"] == "pending"
    assert MAIN in entry["urgent_product"]


def test_other_statuses_are_left_alone(
    workdir: Path, topsapp: Dict[str, Any], services: Dict[str, Any]
) -> None:
    """A backfilled entry (production made by a historic run) stays as it is, even if an event joins."""
    entry = queue_product(workdir, AFTERSHOCK)
    entry["status"] = "backfill"
    queue_file().write_text(json.dumps(entry))
    queue_product(workdir, MAIN)
    api.precise.check_precise_queue()
    entry = read_json(queue_file())
    assert entry["status"] == "backfill" and services["owners"] == [] and topsapp["calls"] == []


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
    assert services["owners"] == [] and lock.exists()

    lock.unlink()
    services["orbit"] = "RESORB"
    api.main_precise()
    assert services["owners"] and not lock.exists()


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
    workdir: Path,
    topsapp: Dict[str, Any],
    services: Dict[str, Any],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import requests

    queue_product(workdir)

    def unreachable(scene: str) -> None:
        raise requests.ConnectionError("no route to host")

    monkeypatch.setattr(api.precise, "orbit_type", unreachable)
    api.precise.check_precise_queue()
    entry = read_json(queue_file())
    assert entry["status"] == "pending" and topsapp["calls"] == []
