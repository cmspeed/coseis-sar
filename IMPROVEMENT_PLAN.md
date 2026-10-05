# COSEIS Improvement Plan

**Goal:** a single `main` that does SAR and optical, is modular and tested, and supports optical in forward mode. Get there without touching `main` until the full set of changes has been tested.

## Ground rules
1. **`main` is frozen for code** until the cutover (Phase 5). Only bot commits (`active_jobs/`, maps) and urgent SAR hotfixes go to `main`.
2. **`develop` is the integration branch.** Every phase is one or more GitHub issues. Each issue gets a branch named by its issue number, cut from `develop`, with a PR back into `develop` (merge commit, not squash/rebase). No PR targets `main` until Phase 5.
3. **Hotfix flow:** fix on `main` in a small PR, then bring it to `develop` right away.
   - **Since Phase 2b, git can't carry `coseis.py` edits across.** `main` still has the single `scripts/coseis.py`, while `develop`'s code lives in `src/aria_coseis/`.
   - So merge `main` into `develop` (keeping `develop`'s `scripts/coseis.py` shim), then **port the fix by hand** into the matching `aria_coseis` module, with a test that covers it.
   - Keep `main` hotfixes rare and small until cutover.
4. **Sync `main` into `develop` regularly** (at least before starting each phase) to pick up hotfixes. Bot data files will conflict; resolve them by taking `main`'s version.
5. **Commits within a PR are small and single-purpose.** In refactor PRs, a commit that moves code never also changes behavior.
6. **For SAR and shared code, `main` wins.** `main` is the stable reference for SAR and for logic both modes share. `develop` may add optical-only behavior but must not change SAR results. *(decided 2026-10-01)*
7. **SAR/forward output changes must be deliberate.** The Phase 2a tests check output against `main`'s production states (`tests/fixtures/production/`) and golden files. A change that alters them must update them in the same PR, with the reason in the commit message.
8. **Recorded HTTP is a snapshot.** Tests replay USGS/ASF/Copernicus as recorded. Re-recording a cassette (e.g. after USGS revises an event) can change golden files; review those diffs as data updates, not regressions.
9. **Keep the external interface stable.** Cron and GitHub Actions call `cd scripts && python coseis.py --forward ...`. That command, its flags, the `active_jobs/` format and the email env vars must keep working through every phase.

## Status
| Phase | State | Branch / PR |
|---|---|---|
| 1. Reconcile `develop` with `main` | **done** 2026-10-01 | issue #27 → PR #28 |
| 2. Tests + modularize | 2a, 2b **done** 2026-10-05; 2c in progress | 2a: issue #29 → PR #31; 2b: issue #30 → PR #32; 2c: branch `2c-cleanup` |
| 3. Optical forward mode | not started | — |
| 4. Test suite + CI | mostly covered by 2a; gaps remain | — |
| 5. Validation + cutover | not started | — |

Completed tasks are checked off in place and tagged `DONE <date> (<commit/PR>)`.

### Why `main` stays safe
- Scheduled GitHub Actions only run on the default branch (`main`), so PRs and pushes to `develop` never trigger the email workflow.
- The local cron runs `git checkout main && git pull` in its own clone, so branch work done in another clone (e.g. this laptop) can't affect it.
- If optical work must happen on the processing machine, use a separate worktree (`git worktree add ../coseis-dev develop`) so cron's `git checkout main` doesn't switch away from it.
- CI (`.github/workflows/tests.yml`) runs only on pull requests and pushes to `develop`.

## Phase 1: Reconcile `develop` with `main`
Make `develop` a strict superset of `main` (all of `main`'s SAR and forward fixes plus the optical work) before doing anything else. Refactoring two diverged codebases is much harder than refactoring one.

- [x] Merge `main` into `develop`. For `coseis.py`, keep all of `main`'s SAR and forward logic and layer `develop`'s optical additions on top. **DONE 2026-10-01 (c167d73)**
  - `coseis.py` was rebuilt from `main`'s version, not by hand-resolving markers. Every function is byte-identical to `main` or `develop` except seven deliberately merged ones: `add_to_tracker`, `check_significance`, `find_reference_and_secondary_pairs`, `process_earthquake`, `main_forward`, CLI, imports.
  - FFM AOI buffer (0.15°, from `develop` e30e7e9) is applied to **optical only**, so SAR frame selection is unchanged from `main`. *(confirmed 2026-10-01)*
  - `--forward` rejects `--sensor` other than `sar` until Phase 3.
- [x] Adopt `main`'s versions of: the `scripts/active_jobs/` tracker directory, `--process_only`, `run_coseis_forward.sh`, `coseis-cron.yml`, `test-email.yml`, and the 30 m forward default. Remove the stale `scripts/active_job_tracking.json`. **DONE 2026-10-01 (c167d73)**
  - The `--resolution` default is now 30 m for all modes (`develop` had 90). Revisit if optical products can't be generated at 30 m.
- [x] Email tiers: primary and secondary only. Remove TERTIARY (`COSEIS_TERTIARY_RECIPIENTS`) everywhere. *(decided 2026-10-01)* **DONE 2026-10-01 (c167d73)**
- [x] Significance criteria *(decided 2026-10-01)*. **DONE 2026-10-01 (c167d73)**
  - historic: M≥6.0, depth ≤40 km, ≤0.5° from coast. The USGS alert-level requirement was removed *(decided 2026-10-01; 7af333d)*
  - forward: keep the M5.5 rule for now, i.e. (M≥5.5 and ≤15 km) or (M≥6.0 and ≤40 km)
  - optical (historic): additionally requires a strike-slip rake (within 45° of 0°/180°). The rake is fetched for optical only, so SAR makes no extra USGS requests (with a 30 s request timeout).
- [x] **Lazy-import optical dependencies** (`ee`, `pystac_client`, `google.cloud.storage`, `osgeo.gdal`). Verified `import coseis` works in the SAR-only `coseis-sar` env. **DONE 2026-10-01 (d4bd3a7)**
- [x] `batch_autorift.py`: hardcoded test path replaced by `--data_dir` (default `scripts/data`, where `coseis.py` writes when run from `scripts/`). **DONE 2026-10-01 (adbb0a5)**
- [x] `environment.yml`: added `earthengine-api`, `google-cloud-storage`, `pystac-client`, `gdal`, and `next_pass` (pip). **DONE 2026-10-01 (63112c0)**
  - `hyp3_autorift` is still undocumented. Decide whether `batch_autorift.py` gets its own environment.
- [x] Equivalence check: historic SAR `--job_list` for 2025-01-07 (Tibet, M7.1) on `main` vs this branch. Job list, AOI, significance CSV/GeoJSON and earthquake info are identical; the only difference is the new `event_id` field on each job (kept for SAR and optical, decided 2026-10-01). Re-verified after the rake change. **DONE 2026-10-01**
- [x] Forward-mode equivalence: the Phase 2a forward and historic SAR tests pass against `main`'s `coseis.py` (apart from `event_id` and the AOI file name). **DONE 2026-10-05 (Phase 2a)**
  - Expected difference: `_partial.json` jobs now include `event_id`.
- [x] Open the Phase 1 issue (#27), rename the branch to `27`, and push. **DONE 2026-10-01**
- [x] Open the PR from `27` into `develop`. **DONE 2026-10-01 (PR #28, merge commit 2cbf4ea)**

## Phase 2: Safety net, then modularize `coseis.py`
Two issues, two PRs into `develop`. 2a must merge before 2b starts.

**2a. Characterization tests (no changes to `coseis.py`).** Issue #29 → PR #31. Pin down current behavior so the refactor can be checked against it. This is the minimum needed to refactor safely, not the full test suite (that's Phase 4). **DONE 2026-10-05 (PR #31)**: 41 tests, about 3 s offline.
- [x] `pytest` scaffold; `.github/workflows/tests.yml` runs on `pull_request` and on pushes to `develop`, never on a schedule. It uses Python 3.11/3.12 and the same packages as `coseis-cron.yml`.
- [x] HTTP recorded once and replayed (`pytest-recording` / vcrpy, one cassette per test module in `tests/cassettes/`, ~6 MB). Requests missing from a cassette fail.
  - Record new tests with `pytest --record-mode=new_episodes`.
  - Re-record a module by deleting its cassette.
  - Refresh golden files with `UPDATE_GOLDEN=1 pytest`.
- [x] Fakes at the outer boundaries only: `subprocess.run` (topsApp), `yagmail.SMTP`, the `next_pass` module. Optical deps aren't needed.
- [x] Tests import code only through `tests/coseis_api.py`; 2b updates that adapter, not the tests.
- [x] What's covered:
  - pure helpers
  - `check_significance` (SAR fails the test if any rake request is made)
  - historic SAR job lists and AOIs (Tibet 2025, Ende 2026, and the non-job-list path)
  - historic optical: Sentinel-2/Copernicus for Kahramanmaraş 2023 (rake filter, buffered FFM AOI, pairing)
  - forward mode
  - CLI routing, including the two exact cron command lines
- [x] Forward tests are seeded from and checked against **real production states** (`tests/fixtures/production/`):
  - The Tamarindo discovery reproduces production commit ff08a73's tracker entry and partial jobs.
  - Ende D61 processing reproduces 92edeb3's READY_FOR_EMAIL state.
  - Also covered: no post-event data yet, topsApp failure → FAILED_NEEDS_ATTENTION, the email run removing the event (as in 595ed24), the runners ignoring each other's states, and the lock file.
- [x] One-off check (not committed): `main`'s `coseis.py` passes the same forward and historic SAR tests once `event_id` and the AOI file name are normalized.
- Not covered in 2a (moved to Phase 4): the Element84 backend, the GEE export/download/manifest path, `batch_autorift.py`.

**2b. Refactor into a package.** **DONE 2026-10-05 (issue #30 → PR #32)**: one PR, one commit per moved module, all tests green after every commit. Uses the standard **src layout** with import name `aria_coseis`, independent of any GitHub repo rename (see backlog).
```
<repo root>/
  pyproject.toml            # name = "aria-coseis"
  src/aria_coseis/
    config.py               # URLs, thresholds, paths, lock file, recipients, GitHub Pages base URL
    utils.py                # to_snake_case, convert_time
    usgs.py  significance.py  aoi.py  notify.py
    sar/      search.py  pairing.py  topsapp.py
    optical/  copernicus.py  element84.py  gee.py  jobs.py
    tracking.py  pipeline.py  modes.py  cli.py
  tests/
  scripts/coseis.py         # shim: adds src/ to sys.path, calls aria_coseis.cli.main()
```
- [x] `cd scripts && python coseis.py ...` is unchanged: cron, GitHub Actions, `run_coseis_forward.sh` and the `coseis-sar` env need no edits, and nothing is installed in production. A test runs the shim from `scripts/` as a subprocess; a live historic run through it matches the golden output.
- [x] Moved in this order: skeleton (everything in `legacy.py`) → config → utils → notify → usgs → aoi → significance → sar → optical → tracking → pipeline → modes → delete `legacy.py`.
- [x] Moves are verbatim. A script checked all 49 functions against the original `coseis.py`; the only edit is reading `root_dir`, `TRACKING_DIR` and the recipients as `config.X`.
- [x] `GITHUB_PAGES_BASE_URL` moved to `config.py`.
- [x] `COSEIS_DATA_DIR`, `COSEIS_TRACKING_DIR`, `COSEIS_LOCK_FILE` overrides, with defaults unchanged. This enables the Phase 5 shadow run without sharing production's lock file.

**2c. Cleanup (behavior-neutral).** Kept out of 2b so its diff stays pure relocation. Branch `2c-cleanup`.
- [ ] `ruff format` + `ruff check` (PEP 8) across `src/`, `tests/` and `scripts/coseis.py`, in one formatting-only commit; add `ruff` to CI. Other scripts (`batch_autorift.py`, the job-list utilities, `s3_upload_aria_share/`) are out of scope.
- [ ] Type hints on all public functions, one commit per module.
- [ ] Fix stale docstrings (depth limits, "60-day windows", "AOI.geojson").
- [ ] Remove dead code: `scripts/coseis_sar.py` (old CLI, still described in `README.md`), `check_for_new_data` (never called), commented-out blocks. Decide separately about `create_directories_from_json` (see backlog).
- [ ] Rewrite `README.md` for the package layout and current CLI.

## Phase 3: Optical in forward mode
Open design questions to settle in the issue before writing code:
- **Tracker schema:** per-sensor entries under each event (e.g. `tracks.sar[...]`, `tracks.optical[...]`) with independent states. The existing SAR entries must stay readable.
- **Trigger:** what counts as "post-event data ready"? First cloud-free S2/Landsat scene over the AOI, or N days of composite window? GEE composites need a post-event window, so near-real-time optical will be partial.
- **Filtering:** apply the strike-slip rake filter in forward mode? (Rake is often unavailable in the first hours after an event.)
- **Where it runs:** GEE export plus autoRIFT runs locally (like topsApp) or elsewhere, and what credentials that needs.
- **Notifications:** optical results in the existing "processing complete" email or in a separate one.
- [ ] Implement behind an opt-in flag (e.g. `--sensors sar,optical`) that defaults to SAR only, so cutover behavior is unchanged unless enabled.

## Phase 4: Test suite and CI
Most of the planned Phase 4 tests were written in 2a:
- pure helpers, pairing modes, the rake filter
- the full tracker state machine, including failure → FAILED_NEEDS_ATTENTION
- email recipients and subjects
- CLI routing
- CI on pull requests

Remaining gaps:
- [ ] Element84 backend (needs `pystac-client` in CI).
- [ ] GEE path with `ee`/GCS faked: composite export → download → merge → manifest, plus Landsat mission/band selection.
- [ ] `batch_autorift.py` on a small raster pair.
- [ ] S1C/S1D acquisition-date limits in `get_SLCs` (a small synthetic ASF response).
- [ ] Email body rendering (beyond the subject/recipient checks that exist).
- [ ] An opt-in integration test (marked, skipped in CI) that hits live USGS/ASF for one event.
- Lint (`ruff`) in CI moved to 2c.

## Phase 5: Validation and cutover to `main`
- [ ] **Shadow-run** `develop`'s forward mode on the processing machine for 1–2 weeks:
  - separate worktree
  - separate tracking dir
  - `--do_processing` against a scratch data dir
  - emails to the test-recipient workflow only
  - no git push
- [ ] Compare its tracker and outputs against `main`'s production run for the same events.
- [ ] Run historic SAR and optical on known events and compare against existing products.
- [ ] Cutover:
  - tag `main` (`pre-integration`)
  - pick a quiet window (no lock file, no events in `READY_FOR_EMAIL`)
  - merge `develop` into `main` via PR
  - update the GitHub Action's pip dependencies if needed
  - watch the next few cron cycles
- [ ] Rollback: revert the merge commit on `main`. The tracker format is unchanged, so no data migration is needed.
- [ ] Update `CLAUDE.md` and `README.md` to describe the single unified branch.

## Backlog: longer-term improvements
Proposed along the way and not yet scheduled. Move items into a phase when picked up.

**Operations / reliability**
- [ ] The forward lock file is configurable in Python (`COSEIS_LOCK_FILE`, Phase 2b), but `run_coseis_forward.sh` still checks the hardcoded default. Have the script read the same variable.
- [ ] **Pin GitHub Actions dependencies.** The email workflow installs unpinned packages and `next_pass` from git HEAD; an upstream change already broke it once (`main` e71d3d0 "Fix next_pass import path"). Use a `requirements-actions.txt` with versions and a `next_pass` commit SHA.
- [ ] **Stale-lock detection.** If the local run is killed (SIGKILL, reboot), `/tmp/coseis_processing.lock` survives and both the bash script and Python silently skip every later run. Store the PID and timestamp in the lock and clear it if the process is gone or the lock is older than N hours.
- [ ] **Alert on `FAILED_NEEDS_ATTENTION`** and on repeated cron failures (e.g. an email to secondary recipients), instead of relying on someone reading `log_tracking.txt`.
- [ ] **Push-race retry.** The GitHub Action and the local cron both `pull --rebase && push` to `main`, so a simultaneous push fails that run's commit. Add a retry loop.
- [ ] Rotate or trim `scripts/log_tracking.txt`.
- [ ] Fix the misleading cron comment: `*/50` runs at :00 and :50, not every 50 minutes.

**Code quality**
- [ ] **Near-land filter buffers twice.** `get_coastline` buffers the land polygons by 0.5°, then `withinCoastline` buffers that again by 0.5° for every earthquake. So the effective distance is about 1°, not the documented 0.5°, and the repeated buffer is slow. Decide on the intended distance, buffer once, and update the significance tests.
- [ ] Forward discovery fetches the same ASF frame search twice per event (`main_forward`, then `add_to_tracker`). Pass the result through.
- [ ] Historic mode without `--job_list` only writes pair JSON and frame maps; `create_directories_from_json` is never called (dead code), and topsApp runs only in forward mode. Remove it, or wire up local historic processing if that's wanted.
- [x] Untrack `scripts/__pycache__/` (a tracked `.pyc` made every local import dirty the tree). **DONE 2026-10-05 (Phase 2a branch)**
- [ ] Replace `print` with `logging` (levels, timestamps), especially for cron logs.
- [ ] Move tunables (magnitude/depth thresholds, coastline buffer, date windows, rake tolerance) into `aria_coseis.config` (endpoints, paths, recipients and the cloud threshold are there since 2b). Today they're scattered literals, and docstrings already disagree with the code.
- [ ] CLI cleanup: the forward branch repeats `--job_list`, `--dates` and `--aoi` checks twice, and `--pairing` is required for forward even though only `coseismic` is used.
- [ ] SAR historic AOI output is named `<title>_sar_toa_AOI.geojson` (an optical level in a SAR filename). Name it by sensor only.
- [ ] Console entry point (`aria-coseis` command). `pyproject.toml` and `pip install -e .` exist since 2b; production still runs through the `scripts/coseis.py` shim.

**Repo**
- [ ] Rename the GitHub repo to `aria-coseis` (from `cmspeed/coseis-sar`).
  - Pages URLs **don't redirect**, so the same day: hotfix `GITHUB_PAGES_BASE_URL` on `main` (it's in `config.py` after 2b), merge `main` into `develop`, and update `origin` on both machines.
  - Links in already-sent emails break unless a stub `coseis-sar` Pages site redirects.
  - Keep the `coseis-sar` mamba env name and the processing machine's clone folder, which cron and `run_coseis_forward.sh` reference.
- [ ] `.gitignore` lists `scripts/run_coseis_forward.sh`, but the file is tracked. Remove the stale entry or decide whether it should be per-machine.

**Science / products**
- [ ] Record per-product provenance (`coseis.py` commit SHA, parameters, scene IDs) in the outputs and HyP3 job metadata.
- [ ] Bring `s3_upload_aria_share/` into the pipeline (or document it) so finished products are uploaded automatically.
