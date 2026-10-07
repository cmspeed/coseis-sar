# COSEIS Improvement Plan

**Goal:** a single `main` that does SAR and optical, is modular and tested, and supports optical in forward mode. Get there without touching `main` until the full set of changes has been tested.

## Ground rules
1. **`main` is frozen for code** until the cutover (Phase 5). Only bot commits (`active_jobs/`, maps) and urgent SAR hotfixes go to `main`.
2. **`develop` is the integration branch.** Every phase is one or more GitHub issues. Each issue gets a branch named by its issue number, cut from `develop`, with a PR back into `develop` (merge commit, not squash/rebase). No PR targets `main` until Phase 5.
3. **Hotfix flow:** fix on `main` in a small PR, then bring it to `develop` right away.
   - **Since Phase 2b, git can't carry `coseis.py` edits across.** `main` still has the single `scripts/coseis.py`, while `develop`'s code lives in `src/aria_coseis/`.
   - So **port the fix by hand** into the matching `aria_coseis` module on a `develop` branch, with a test that covers it. Don't merge `main` into `develop` (see rule 4).
   - Keep `main` hotfixes rare and small until cutover.
4. **No routine `main` → `develop` syncs** *(since Phase 2d, decided 2026-10-06)*. `main`'s new commits are bot tracker/map data, which `develop` doesn't carry: its `active_jobs/` holds only placeholders. Hotfixes are ported by hand (rule 3). The live tracker data moves into the new layout at cutover (Phase 5).
5. **Commits within a PR are small and single-purpose.** In refactor PRs, a commit that moves code never also changes behavior.
6. **For SAR and shared code, `main` wins.** `main` is the stable reference for SAR and for logic both modes share. `develop` may add optical-only behavior but must not change SAR results. *(decided 2026-10-01)*
7. **SAR/forward output changes must be deliberate.** The Phase 2a tests check output against `main`'s production states (`tests/fixtures/production/`) and golden files. A change that alters them must update them in the same PR, with the reason in the commit message.
8. **Recorded HTTP is a snapshot.** Tests replay USGS/ASF/Copernicus as recorded. Re-recording a cassette (e.g. after USGS revises an event) can change golden files; review those diffs as data updates, not regressions.
9. **Keep the external interface stable.** The CLI flags, the tracker JSON format and the email env vars must keep working through every phase. Since Phase 2d, `develop` runs as `python -m aria_coseis ...` from the repository root (`ops/run_coseis_forward.sh`, `coseis-cron.yml`); production keeps `cd scripts && python coseis.py ...` on `main` until the cutover switches both runners together.

## Status
| Phase | State | Branch / PR |
|---|---|---|
| 1. Reconcile `develop` with `main` | **done** 2026-10-01 | issue #27 → PR #28 |
| 2. Tests, package, layout | **done** (2a–2c 2026-10-05, 2d 2026-10-06) | 2a: #29 → PR #31; 2b: #30 → PR #32; 2c: #33 → PR #34; 2d: #35 → PR #36 |
| 3. Optical forward mode | not started | — |
| 4. Test suite + CI | mostly covered by 2a; gaps remain | — |
| 5. Validation + cutover | cutover in progress 2026-10-07 | branch `cutover` → PR into `main` |

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

**2c. Cleanup (behavior-neutral).** Kept out of 2b so its diff stays pure relocation. **DONE 2026-10-05 (issue #33)**
- [x] `ruff format` (PEP 8, 99 columns) over `src/`, `tests/` and `scripts/coseis.py`, in a formatting-only commit; the AST of every file is unchanged apart from docstring whitespace.
- [x] `ruff check` (pycodestyle, pyflakes, import sorting) with auto-fixes only, and CI runs `ruff check` and `ruff format --check`.
  - E501 (line length) is ignored: code fits in 99 columns, while long strings, docstrings and HTML templates are left as is.
  - Other scripts (`batch_autorift.py`, the job-list utilities, `s3_upload_aria_share/`) are out of scope.
- [x] Type hints on all top-level functions, one commit per module, with `from __future__ import annotations`. Verified that the code is identical once annotations are removed.
- [x] Stale docstrings, comments and one log message fixed (near-land ~1°, CLI examples, GEE 90-day windows, ASF ±90-day search).
- [x] Dead code removed: `scripts/coseis_sar.py`, `check_for_new_data`, a commented-out debug print. `create_directories_from_json` is kept pending a decision (see backlog).
- [x] `README.md` rewritten for the package layout and current CLI.

**2d. Repository layout.** **DONE 2026-10-06 (issue #35 → PR #36)**: GitHub holds only what historic and forward mode need. Decisions (2026-10-06): partial jobs in `active_jobs/partials/`, untracked outputs in `outputs/`, cron wrapper tracked in `ops/`, archive tag plus a local copy, and no routine `main` → `develop` syncs.
```
.github/workflows/   coseis-cron.yml, tests.yml, test-email.yml
src/aria_coseis/     package; python -m aria_coseis (aria-coseis); optical/autorift.py (was scripts/batch_autorift.py)
tests/
active_jobs/         tracker JSONs; partials/ holds partial HyP3 jobs (bot-committed)
docs/maps/           GitHub Pages overpass maps (bot-committed)
ops/                 run_coseis_forward.sh (processing-machine cron wrapper)
pyproject.toml  environment.yml  README.md  IMPROVEMENT_PLAN.md
```
- [x] Dependencies declared in `pyproject.toml`, with optional `forward` (next_pass), `optical` and `test` groups; `requirements-test.txt` removed; `environment.yml` installs the package.
- [x] All file locations come from `config` relative to the repo root (`TRACKING_DIR`, `PARTIALS_DIR`, `MAPS_DIR`, `OUTPUT_DIR`, `root_dir`), each overridable by `COSEIS_*`. Trackers keep bare partial-job file names, so live entries still work after the move.
- [x] `scripts/batch_autorift.py` moved to `aria_coseis.optical.autorift` (`git mv`; default `--data_dir` is now the data dir) and formatted with ruff.
- [x] Cron wrapper moved to `ops/` and updated: repo root, `python -m aria_coseis`, `logs/forward.log`, `COSEIS_LOCK_FILE`, `git add -A active_jobs/`.
- [x] `coseis-cron.yml` installs `.[forward]`, runs from the root, and stages `active_jobs/` and `docs/maps/`. This takes effect only on `main`, at cutover.
- [x] Untracked, but kept locally and gitignored: `scripts/`, `hyp3/`, `job_lists/`, `s3_upload_aria_share/`, `docs/*.txt`, `requirements.txt`. Recoverable from tag `archive/pre-layout` and a local copy in `../coseis-archive/`.
- [x] Verified from a fresh clone: install, ruff, and all tests pass; a live historic run from the root writes only to `outputs/` and matches the golden job list.
- [x] The processing machine doesn't use `hyp3/`, `job_lists/` or `s3_upload_aria_share/` (confirmed 2026-10-06), so nothing there needs a backup at cutover.

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
**Decision (2026-10-06): cut over directly, without a week-long shadow run.** Forward discovery and processing are already checked against real production runs (Phase 2a), and they also passed against `main`'s own code. A smoke test on the processing machine covers its environment instead, and a brief cron hiccup is acceptable. The significance criteria are deliberately left as they are, to be revisited later together with the manuscript.

- [x] Health check: no job running, no lock file; the next Tamarindo acquisition (2026-10-08) is after the cutover. **DONE 2026-10-07**
- [x] Smoke test on the processing machine: temporary clone of `develop` run with `PYTHONPATH=src` and the `COSEIS_*` overrides on a scratch copy of the live tracker, `--forward --process_only --do_processing`. Ran cleanly: no post-event data, tracker unchanged. **DONE 2026-10-07**
- [x] Paused the local cron and disabled the `coseis-cron.yml` workflow. **DONE 2026-10-07**
- [x] Tagged `main` as `pre-cutover` (595ed24). **DONE 2026-10-07**
- [x] Cutover branch from `main`:
  - `develop` merged with an explicit merge commit, so the cutover can be reverted as a single commit
  - the merge removes `scripts/`, so the live tracker (`us6000tymj.json`) and its two partial jobs were restored from `pre-cutover`, byte-identical, and moved into `active_jobs/` and `active_jobs/partials/`
  - verified that each tracker entry's partial file resolves, and that a `--process_only` run through the new code on a copy leaves the tracker unchanged

  **DONE 2026-10-07**
- [ ] Cutover PR into `main`: CI green, then merge with a merge commit.
- [ ] Processing machine:
  - `git pull` on `main`
  - `pip install -e . --no-deps` in `coseis-sar` (`--no-deps` so pip doesn't replace conda-installed packages)
  - check `python -m aria_coseis --help`
  - crontab: `COSEIS_DATA_DIR=<repo>/scripts/data <repo>/ops/run_coseis_forward.sh`, keeping the existing topsApp products in place
  - run the wrapper once by hand; check `logs/forward.log` and `git log`
- [ ] Re-enable the `coseis-cron.yml` workflow and trigger it once by hand (exercises discovery and the `.[forward]` install); watch the next few cycles.
- [ ] First real processing on the new code: Tamarindo (A165/D157) once the 2026-10-08 acquisition is available.
- [ ] Rollback if needed: revert the cutover merge commit on `main` (this also reverses the tracker moves) and point the crontab back at `scripts/run_coseis_forward.sh`.
- [ ] After cutover: decide whether to keep `develop` as the integration branch or merge issue branches directly into `main`; update `CLAUDE.md` and the ground rules accordingly.
- Deferred: historic SAR and optical comparisons against existing products (covered for SAR by the golden tests; optical with Phase 4).

## Backlog: longer-term improvements
Proposed along the way and not yet scheduled. Move items into a phase when picked up.

**Operations / reliability**
- [x] The forward lock file is configurable in both Python and `ops/run_coseis_forward.sh` (`COSEIS_LOCK_FILE`). **DONE 2026-10-06 (Phase 2d)**
- [ ] **Pin forward-mode dependencies.** The email workflow installs `.[forward]`, which pulls `next_pass` from git HEAD and unpinned core packages; an upstream change already broke it once (`main` e71d3d0 "Fix next_pass import path"). Pin `next_pass` to a commit or release in `pyproject.toml` and add a constraints file for the workflow.
- [ ] **Stale-lock detection.** If the local run is killed (SIGKILL, reboot), `/tmp/coseis_processing.lock` survives and both the bash script and Python silently skip every later run. Store the PID and timestamp in the lock and clear it if the process is gone or the lock is older than N hours.
- [ ] **Alert on `FAILED_NEEDS_ATTENTION`** and on repeated cron failures (e.g. an email to secondary recipients), instead of relying on someone reading `log_tracking.txt`.
- [ ] **Push-race retry.** The GitHub Action and the local cron both `pull --rebase && push` to `main`, so a simultaneous push fails that run's commit. Add a retry loop.
- [ ] Rotate or trim `logs/forward.log` (formerly `scripts/log_tracking.txt`).
- [x] Fix the misleading cron comment: `*/50` runs at :00 and :50, not every 50 minutes. **DONE 2026-10-06 (Phase 2d)**

**Code quality**
- [ ] `--eq_list` runs skip `check_significance`, so custom event lists are processed as given. Decide whether lists should be filtered too.
- [ ] Single-date historic queries build `starttime` as `YYYY-MM-DD00:00:00` (no separator). USGS accepts it, but it's malformed; use `YYYY-MM-DDT00:00:00`.
- [ ] `optical/element84.py` has a bare `except:` (marked `noqa: E722`), which also swallows `KeyboardInterrupt`/`SystemExit`. Narrow it to `except Exception:` with a test.
- [ ] **Near-land filter buffers twice.** `get_coastline` buffers the land polygons by 0.5°, then `withinCoastline` buffers that again by 0.5° for every earthquake. So the effective distance is about 1°, not the documented 0.5°, and the repeated buffer is slow. Decide on the intended distance, buffer once, and update the significance tests.
- [ ] Forward discovery fetches the same ASF frame search twice per event (`main_forward`, then `add_to_tracker`). Pass the result through.
- [ ] Historic mode without `--job_list` only writes pair JSON and frame maps; `create_directories_from_json` is never called (dead code), and topsApp runs only in forward mode. Remove it, or wire up local historic processing if that's wanted.
- [x] Untrack `scripts/__pycache__/` (a tracked `.pyc` made every local import dirty the tree). **DONE 2026-10-05 (Phase 2a branch)**
- [ ] Replace `print` with `logging` (levels, timestamps), especially for cron logs.
- [ ] Move tunables (magnitude/depth thresholds, coastline buffer, date windows, rake tolerance) into `aria_coseis.config` (endpoints, paths, recipients and the cloud threshold are there since 2b). Today they're scattered literals, and docstrings already disagree with the code.
- [ ] CLI cleanup: the forward branch repeats `--job_list`, `--dates` and `--aoi` checks twice, and `--pairing` is required for forward even though only `coseismic` is used.
- [ ] SAR historic AOI output is named `<title>_sar_toa_AOI.geojson` (an optical level in a SAR filename). Name it by sensor only.
- [x] Console entry point: `aria-coseis` and `python -m aria_coseis`. **DONE 2026-10-06 (Phase 2d)**

**Repo**
- [ ] Rename the GitHub repo to `aria-coseis` (from `cmspeed/coseis-sar`).
  - Pages URLs **don't redirect**, so the same day: hotfix `GITHUB_PAGES_BASE_URL` on `main` (it's in `config.py` after 2b), merge `main` into `develop`, and update `origin` on both machines.
  - Links in already-sent emails break unless a stub `coseis-sar` Pages site redirects.
  - Keep the `coseis-sar` mamba env name and the processing machine's clone folder, which cron and `run_coseis_forward.sh` reference.
- [x] Stale `.gitignore` entry for the tracked cron wrapper removed; the wrapper now lives in `ops/`. **DONE 2026-10-06 (Phase 2d)**

**Science / products**
- [ ] Record per-product provenance (`coseis.py` commit SHA, parameters, scene IDs) in the outputs and HyP3 job metadata.
- [ ] Bring `s3_upload_aria_share/` into the pipeline (or document it) so finished products are uploaded automatically.
