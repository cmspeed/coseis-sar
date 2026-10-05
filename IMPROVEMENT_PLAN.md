# COSEIS Improvement Plan

**Goal:** a single `main` that does SAR and optical, is modular and tested, and supports optical in forward mode. Get there without touching `main` until the full set of changes has been tested.

## Ground rules
1. **`main` is frozen for code** until the cutover (Phase 5). Only bot commits (`active_jobs/`, maps) and urgent SAR hotfixes go to `main`.
2. **`develop` is the integration branch.** Every phase is one or more GitHub issues. Each issue gets a branch named by its issue number, cut from `develop`, with a PR back into `develop`. No PR targets `main` until Phase 5.
3. **Hotfix flow:** fix on `main` in a small PR, then merge `main` into `develop` right away so the fix isn't lost.
4. **Sync `main` into `develop` regularly** (at least before starting each phase) to pick up hotfixes. Bot data files will conflict; resolve them by taking `main`'s version.
5. **Commits within a PR are small and single-purpose.** In refactor PRs, a commit that moves code never also changes behavior.
6. **For SAR and shared code, `main` wins.** `main` is the stable reference for SAR and for logic both modes share. `develop` may add optical-only behavior but must not change SAR results. *(decided 2026-10-01)*
7. **Keep the external interface stable.** Cron and GitHub Actions call `cd scripts && python coseis.py --forward ...`. That command, its flags, the `active_jobs/` format and the email env vars must keep working through every phase.

## Status
| Phase | State | Branch / PR |
|---|---|---|
| 1. Reconcile `develop` with `main` | **done** 2026-10-01 | issue #27 → PR #28 |
| 2. Tests + modularize | in progress (2a) | 2a: branch `phase2a` (rename to issue #); 2b: not started |
| 3. Optical forward mode | not started | — |
| 4. Test suite + CI | not started | — |
| 5. Validation + cutover | not started | — |

Completed tasks are checked off in place and tagged `DONE <date> (<commit/PR>)`.

### Why `main` stays safe
- Scheduled GitHub Actions only run on the default branch (`main`), so PRs and pushes to `develop` never trigger the email workflow.
- The local cron runs `git checkout main && git pull` in its own clone, so branch work done in another clone (e.g. this laptop) can't affect it.
- If optical work must happen on the processing machine, use a separate worktree (`git worktree add ../coseis-dev develop`) so cron's `git checkout main` doesn't switch away from it.

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
- [ ] Forward-mode equivalence is still unverified offline. **Moved to Phase 2a**, which checks it with the recorded fixtures.
  - Expected difference: `_partial.json` jobs now include `event_id`.
- [x] Open the Phase 1 issue (#27), rename the branch to `27`, and push. **DONE 2026-10-01**
- [x] Open the PR from `27` into `develop`. **DONE 2026-10-01 (PR #28, merge commit 2cbf4ea)**

## Phase 2: Safety net, then modularize `coseis.py`
Two issues, two PRs into `develop`. 2a must merge before 2b starts.

**2a. Characterization tests (no changes to `coseis.py`).** Pin down current behavior so the refactor can be checked against it. This is the minimum needed to refactor safely, not the full test suite (that's Phase 4).
- [ ] `pytest` scaffold and a GitHub workflow that runs it on `pull_request` and on pushes to `develop` (never on a schedule).
- [ ] **Record HTTP once, replay forever** (`pytest-recording` / vcrpy cassettes in `tests/cassettes/`). Almost every external call is `requests.get` (USGS, FFM/rake, coastline, ASF, Copernicus, Element84). Mocking at the HTTP layer means the tests survive functions moving between modules in 2b.
- [ ] Stub only at the outer boundaries:
  - `subprocess.run` (topsApp)
  - `yagmail.SMTP`
  - the `next_pass` module
  - GEE/GCS
- [ ] Tests reach `coseis` code only through one adapter (`tests/coseis_api.py`), so 2b updates the adapter, not the tests.
- [ ] Events:
  - Tibet 2025 (`us6000pi9w`)
  - Ende M7.7 (`us6000tkt2`)
  - Tamarindo M5.6 (`us6000tymj`, the forward M5.5 rule)
  - a strike-slip optical event
  - synthetic events that should be rejected
- [ ] Golden-output tests:
  - pure helpers
  - `check_significance` (SAR never fetches rake)
  - historic SAR job list and AOI
  - forward tracker and `_partial.json`
  - the tracker state machine (success, failure, email)
  - CLI routing, including the two exact cron command lines
- [ ] One-off check (not committed) that `main`'s `coseis.py` produces the same forward tracker output from the same fixtures, apart from `event_id`. This closes the open Phase 1 item.

**2b. Refactor into a package** (one PR, one commit per moved module, 2a tests green after every commit). Use the standard **src layout** with import name `aria_coseis`. This is independent of any GitHub repo rename (see backlog). Don't put a `coseis/` package next to `scripts/coseis.py`: the two would collide on `import coseis`.
```
<repo root>/
  pyproject.toml            # name = "aria-coseis"
  src/aria_coseis/
    config.py               # URLs, thresholds, paths, env/recipients, GitHub Pages base URL
    usgs.py  significance.py  aoi.py
    sar/      search.py  pairing.py  topsapp.py
    optical/  copernicus.py  element84.py  gee.py  jobs.py
    tracking.py  notify.py  pipeline.py  modes.py  cli.py
  tests/
  scripts/coseis.py         # ~10-line shim: add src/ to sys.path, call aria_coseis.cli.main()
```
- `cd scripts && python coseis.py ...` stays unchanged, so cron, GitHub Actions and the `coseis-sar` env need no edits and nothing is installed in production. Tests and CI use `pip install -e .`.
- [ ] Order, from fewest to most dependencies: skeleton + shim (everything in `legacy.py`) → config → notify → usgs/aoi → significance → sar → optical → tracking → pipeline/modes → delete `legacy.py`.
- [ ] Move commits are verbatim relocations: explicit imports only, read settings as `config.X` at call time. Type hints, PEP 8 / `ruff format` and docstring fixes come in separate commits at the end.
- [ ] Make the tracking dir and data root configurable (`COSEIS_TRACKING_DIR`, `COSEIS_DATA_DIR`; defaults unchanged). This enables the Phase 5 shadow runs.
- [ ] Remove dead code (`coseis_sar.py`, commented-out blocks) and fix stale docstrings, in separate commits.

## Phase 3: Optical in forward mode
Open design questions to settle in the issue before writing code:
- **Tracker schema:** per-sensor entries under each event (e.g. `tracks.sar[...]`, `tracks.optical[...]`) with independent states. The existing SAR entries must stay readable.
- **Trigger:** what counts as "post-event data ready"? First cloud-free S2/Landsat scene over the AOI, or N days of composite window? GEE composites need a post-event window, so near-real-time optical will be partial.
- **Filtering:** apply the strike-slip rake filter in forward mode? (Rake is often unavailable in the first hours after an event.)
- **Where it runs:** GEE export plus autoRIFT runs locally (like topsApp) or elsewhere, and what credentials that needs.
- **Notifications:** optical results in the existing "processing complete" email or in a separate one.
- [ ] Implement behind an opt-in flag (e.g. `--sensors sar,optical`) that defaults to SAR only, so cutover behavior is unchanged unless enabled.

## Phase 4: Test suite and CI
- [ ] Unit tests for pure logic:
  - date windows
  - Landsat mission and band selection
  - rake filter
  - `to_snake_case`
  - pairing modes (`all` / `sequential` / `coseismic`)
  - S1C/S1D date limits
  - UTM EPSG calculation
- [ ] Tracker state-machine tests:
  - `AWAITING_POST_SEISMIC → READY_FOR_EMAIL → removed`
  - the failure path to `FAILED_NEEDS_ATTENTION`
  - orphaned partial cleanup
- [ ] Email rendering tests (HTML builds, no sending; `yagmail` mocked).
- [ ] An opt-in integration test (marked, skipped in CI) that hits live USGS/ASF for one event.
- [ ] Lint (`ruff`) in CI.

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
- [ ] **topsApp failures are recorded as "Success".** `run_dockerized_topsApp` catches `CalledProcessError` and doesn't re-raise it, so `check_tracker_for_updates` sets `processing_status = "Success"`. A failed run then never reaches `FAILED_NEEDS_ATTENTION` and is emailed as a success; only the missing `.nc` warning hints at it. Fix after 2b, updating the characterization test that currently pins this behavior.
- [ ] The forward lock file path is hardcoded (`/tmp/coseis_processing.lock`) in both `main_forward` and `run_coseis_forward.sh`. Make it configurable alongside the tracking dir.
- [ ] **Pin GitHub Actions dependencies.** The email workflow installs unpinned packages and `next_pass` from git HEAD; an upstream change already broke it once (`main` e71d3d0 "Fix next_pass import path"). Use a `requirements-actions.txt` with versions and a `next_pass` commit SHA.
- [ ] **Stale-lock detection.** If the local run is killed (SIGKILL, reboot), `/tmp/coseis_processing.lock` survives and both the bash script and Python silently skip every later run. Store the PID and timestamp in the lock and clear it if the process is gone or the lock is older than N hours.
- [ ] **Alert on `FAILED_NEEDS_ATTENTION`** and on repeated cron failures (e.g. an email to secondary recipients), instead of relying on someone reading `log_tracking.txt`.
- [ ] **Push-race retry.** The GitHub Action and the local cron both `pull --rebase && push` to `main`, so a simultaneous push fails that run's commit. Add a retry loop.
- [ ] Rotate or trim `scripts/log_tracking.txt`.
- [ ] Fix the misleading cron comment: `*/50` runs at :00 and :50, not every 50 minutes.

**Code quality**
- [x] Untrack `scripts/__pycache__/` (a tracked `.pyc` made every local import dirty the tree). **DONE 2026-10-05 (Phase 2a branch)**
- [ ] Replace `print` with `logging` (levels, timestamps), especially for cron logs.
- [ ] Move tunables (magnitude/depth thresholds, coastline buffer, date windows, rake tolerance, cloud threshold) into one config module or file. Today they're scattered literals, and docstrings already disagree with the code.
- [ ] CLI cleanup: the forward branch repeats `--job_list`, `--dates` and `--aoi` checks twice, and `--pairing` is required for forward even though only `coseismic` is used.
- [ ] SAR historic AOI output is named `<title>_sar_toa_AOI.geojson` (an optical level in a SAR filename). Name it by sensor only.
- [ ] Package with `pyproject.toml` (installable `coseis`, console entry point) once Phase 2 lands.

**Repo**
- [ ] Rename the GitHub repo to `aria-coseis` (from `cmspeed/coseis-sar`).
  - Pages URLs **don't redirect**, so the same day: hotfix `GITHUB_PAGES_BASE_URL` on `main` (it's in `config.py` after 2b), merge `main` into `develop`, and update `origin` on both machines.
  - Links in already-sent emails break unless a stub `coseis-sar` Pages site redirects.
  - Keep the `coseis-sar` mamba env name and the processing machine's clone folder, which cron and `run_coseis_forward.sh` reference.
- [ ] `.gitignore` lists `scripts/run_coseis_forward.sh`, but the file is tracked. Remove the stale entry or decide whether it should be per-machine.

**Science / products**
- [ ] Record per-product provenance (`coseis.py` commit SHA, parameters, scene IDs) in the outputs and HyP3 job metadata.
- [ ] Bring `s3_upload_aria_share/` into the pipeline (or document it) so finished products are uploaded automatically.
