# Forward-mode operations

How the automated (forward) mode runs in production. Most users only need the main
[README](../README.md); this page is for whoever operates the scheduled runners.

## Two runners, one tracker
Both runners work on the `main` branch and share the tracker in `active_jobs/`:

| Runner | Where | Command | Does |
|---|---|---|---|
| Discovery and alerts | GitHub Actions, `.github/workflows/coseis-cron.yml` (minutes 0 and 50 of every hour) | `python -m aria_coseis --forward --pairing coseismic --send_email` | Queries USGS for events in the last 48 h; starts tracking significant ones; publishes overpass maps to `docs/maps/` (GitHub Pages); emails new-event alerts and processing results; commits with `[skip ci]` |
| Processing | Processing machine, cron → `ops/run_coseis_forward.sh` | `python -m aria_coseis --forward --pairing coseismic --resolution 30 --do_processing --process_only` | Checks ASF for post-event Sentinel-1 data for tracked events; runs topsApp; updates the tracker; commits with `[skip ci]` |

Tracker states, per event and track:
`AWAITING_POST_SEISMIC` → (post-event SLCs found, topsApp run) → `READY_FOR_EMAIL` → result emailed and track removed, or `FAILED_NEEDS_ATTENTION` if processing failed.

`active_jobs/<event_id>.json` holds one event; `active_jobs/partials/` holds its HyP3-style job files waiting for post-event granules.

## Processing machine setup
```bash
git clone git@github.com:cmspeed/coseis-sar.git && cd coseis-sar   # SSH: cron pushes with a key
mamba env create -f environment.yml          # or, in an existing env: pip install -e . --no-deps
crontab -e
```
Crontab entry (runs at minutes 5 and 15):
```
15,05 * * * * COSEIS_DATA_DIR=<repo>/scripts/data <repo>/ops/run_coseis_forward.sh
```
- `COSEIS_DATA_DIR` sets where topsApp products go. The current machine keeps them in `scripts/data/` from before the 2026 layout change; without the variable they go to `<repo>/data/`.
- The wrapper loads conda through `~/.bashrc`, activates `coseis-sar`, pulls `main`, runs processing, and pushes tracker changes. Its log is `logs/forward.log`.
- **Pushing from cron needs SSH.** `origin` must be the SSH URL, with a key that GitHub accepts for writing (currently a deploy key on `coseis-sar` with write access) and no passphrase, since cron has no SSH agent. An HTTPS remote can appear to work in an interactive terminal (for example through an editor's login helper) and still fail under cron. To test it the way cron runs:
  ```bash
  env -i HOME="$HOME" PATH=/usr/bin:/bin GIT_TERMINAL_PROMPT=0 \
    bash -c 'cd <repo> && git fetch origin && git push --dry-run origin main'
  ```
  "Everything up-to-date" means it works.
- topsApp runs through `conda run -n topsapp_env_trappist_python11 isce2_topsapp ...`, so that environment must exist. It also needs NASA Earthdata credentials (e.g. `~/.netrc`) to download SLCs.

## GitHub Actions setup
Repository secrets used by `coseis-cron.yml`:
- `GMAIL_USER`, `GMAIL_APP_PSWD`: sender account and app password
- `COSEIS_PRIMARY_RECIPIENTS`: new-event alerts (comma-separated)
- `COSEIS_SECONDARY_RECIPIENTS`: processing results (comma-separated)

`test-email.yml` (run manually) sends a test message to `TEST_RECIPIENTS`.

## Day-to-day
- **Pause:** comment out the crontab line, and disable the workflow under Actions → *COSEIS Forward Email Notifications* → *Disable workflow*.
- **Resume:** reverse both. Trigger the workflow once by hand (*Run workflow*) to check it.
- **Check health:**
  - `tail logs/forward.log`
  - recent runs on the Actions tab
  - recent `[skip ci]` commits on `main`
- **Stuck runs:** the lock file `/tmp/coseis_processing.lock` (override with `COSEIS_LOCK_FILE`) prevents overlapping runs. If a run was killed, the lock stays behind and every later run exits immediately. Delete it once you're sure nothing is running.
- **Don't hand-edit `active_jobs/` while the runners are active.** Both commit to it.
- **Push failures:** git's own output, including errors, goes to `logs/forward.log` just above the "Push attempt N failed" lines. "Authentication failed", "could not read Username" or "Permission denied (publickey)" means a credentials problem (see setup above), not a conflict.
- **Push conflicts:** both runners retry their push, rebasing onto each other's commits. If the processing machine can't push because of a real conflict (both edited the same tracker lines), `logs/forward.log` shows "WARNING: local tracker commits are not on GitHub yet". In that case, resolve it by hand in the production clone (`git pull --rebase origin main`, fix the JSON, `git rebase --continue`, `git push`) while the cron is paused. A failed push in the Action fails that workflow run.

## Overrides
Every path is relative to the repository root and can be overridden with an environment variable:
- `COSEIS_DATA_DIR`
- `COSEIS_TRACKING_DIR`
- `COSEIS_PARTIALS_DIR`
- `COSEIS_MAPS_DIR`
- `COSEIS_OUTPUT_DIR`
- `COSEIS_LOCK_FILE`

These are useful for a dry run against a scratch copy of the tracker, without touching production.
