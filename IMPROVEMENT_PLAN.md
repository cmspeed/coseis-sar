# COSEIS Improvement Plan

**Goal:** a single `main` that does SAR and optical, is modular and tested, and supports optical in forward mode. Get there without touching `main` until the full set of changes has been tested.

## Ground rules
1. **`main` is frozen for code** until the cutover (Phase 5). Only bot commits (`active_jobs/`, maps) and urgent SAR hotfixes go to `main`.
2. **`develop` is the integration branch.** Every phase is one or more GitHub issues. Each issue gets a branch named by its issue number, cut from `develop`, with a PR back into `develop`. No PR targets `main` until Phase 5.
3. **Hotfix flow:** fix on `main` in a small PR, then merge `main` into `develop` right away so the fix isn't lost.
4. **Sync `main` into `develop` regularly** (at least before starting each phase) to pick up hotfixes. Bot data files will conflict; resolve them by taking `main`'s version.
5. **Commits within a PR are small and single-purpose.** In refactor PRs, a commit that moves code never also changes behavior.
6. **Keep the external interface stable.** Cron and GitHub Actions call `cd scripts && python coseis.py --forward ...`. That command, its flags, the `active_jobs/` format and the email env vars must keep working through every phase.

## Status
| Phase | State | Branch / PR |
|---|---|---|
| 1. Reconcile `develop` with `main` | in progress | — |
| 2. Tests + modularize | not started | — |
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

- [ ] Merge `main` into `develop`. For `coseis.py`, keep all of `main`'s SAR and forward logic and layer `develop`'s optical additions on top.
- [ ] Adopt `main`'s versions of: the `scripts/active_jobs/` tracker directory, `--process_only`, `run_coseis_forward.sh`, `coseis-cron.yml`, `test-email.yml`, and the 30 m forward default. Remove the stale `scripts/active_job_tracking.json`.
- [ ] Email tiers: primary and secondary only. Remove TERTIARY (`COSEIS_TERTIARY_RECIPIENTS`) everywhere. *(decided 2026-10-01)*
- [ ] Significance criteria *(decided 2026-10-01)*:
  - historic: M≥6.0, depth ≤40 km, ≤0.5° from coast
  - forward: keep the M5.5 rule for now, i.e. (M≥5.5 and ≤15 km) or (M≥6.0 and ≤40 km)
  - open question: should historic still require a USGS alert level? (`main` requires it; `develop` commented it out)
- [ ] **Lazy-import optical dependencies.** Today `develop` imports `ee` and `google.cloud.storage` at module level, but the GitHub Action doesn't install them. Merging `develop` into `main` as it stands would break the email workflow.
- [ ] Finish optical stabilization: remove `batch_autorift.py`'s hardcoded test path, make `coseis.py` and `batch_autorift.py` agree on where `data/` lives, and add the optical dependencies to `environment.yml`.
- [ ] Check: `--forward --pairing coseismic` (no flags) runs on `develop` and produces the same tracker JSONs as `main` for the same event.

## Phase 2: Safety net, then modularize `coseis.py`
**2a. Characterization tests (before any refactor).** Pin down current behavior so the refactor can be checked against it. This is the minimum needed to refactor safely; it is not the full test suite (that's Phase 4).
- [ ] Add a `pytest` scaffold and a GitHub workflow that runs it on `pull_request` (it never runs on a schedule, so it can't interfere with cron).
- [ ] Record USGS/ASF/coastline responses for 2–3 known events (e.g. Ierapetra, Khovd, one Japan event) as JSON fixtures, with the network mocked.
- [ ] Add golden-output tests:
  - `check_significance`
  - AOI construction
  - SAR frame selection and pairing (job JSON for HyP3, `_partial.json` for forward)
  - optical pair and manifest generation
  - CLI argument validation

**2b. Refactor into a package** (one PR per extraction, each passing the 2a tests):
```
scripts/
  coseis.py              # thin CLI shim: same flags, same invocation
  coseis/
    config.py            # constants, env vars, paths (tracking dir, data root)
    events.py            # USGS queries, coastline, significance, rake
    aoi.py               # AOI / FFM geometry
    sar/                 # ASF search, frame selection, pairing, job JSON, topsApp
    optical/             # copernicus.py, element84.py, gee.py, manifest
    tracking.py          # active_jobs state machine
    notify.py            # email, HTML, maps, next_pass
    modes/               # historic.py, forward.py
```
- [ ] Suggested order, from fewest to most dependencies: config → notify → events/aoi → sar → optical → tracking → modes/CLI.
- [ ] Make the tracking dir and data root configurable (env var or flag). This enables the Phase 5 shadow runs.
- [ ] Remove dead code (`coseis_sar.py`, commented-out blocks) and fix stale docstrings in separate commits.

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
