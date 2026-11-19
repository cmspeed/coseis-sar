# Forward-mode operations

How the automated (forward) mode runs in production. Most users only need the main
[README](../README.md); this page is for whoever operates the scheduled runners.

## Three runners, one tracker
All runners work on the `main` branch and share the tracker in `active_jobs/`:

| Runner | Where | Command | Does |
|---|---|---|---|
| Discovery and alerts | GitHub Actions, `.github/workflows/coseis-cron.yml` (minutes 0 and 50 of every hour) | `python -m aria_coseis --forward --pairing coseismic --send_email` | Queries USGS for events in the last 48 h; starts tracking significant ones; publishes overpass maps to `docs/maps/` (GitHub Pages); emails new-event alerts, processing results and the precise runner's notifications; commits with `[skip ci]` |
| Urgent processing | Processing machine, cron → `ops/run_coseis_forward.sh` | `python -m aria_coseis --forward --pairing coseismic --resolution 30 --do_processing --process_only` | Checks ASF for post-event Sentinel-1 data for tracked events; runs topsApp (restituted orbits); places the product in the archive's `URGENT_RESPONSE/` and queues it for a precise-orbit rerun; updates the tracker; commits with `[skip ci]` |
| Precise-orbit reruns | Processing machine, cron → `ops/run_coseis_precise.sh` | `python -m aria_coseis --precise` | Reruns queued urgent products once precise orbits exist, for events that meet the significance criteria, and places them in `HISTORIC_EVENTS/`; commits its email notifications with `[skip ci]` |

Only one Action run happens at a time (a `concurrency` group). The two cron jobs share the clone, so every git step in their wrappers holds a lock (`.git/coseis-git.lock`, `flock`); processing runs without it, and each runner has its own processing lock and CPU cores, so a long rerun never delays urgent processing.

Tracker states, per event and track:
`AWAITING_POST_SEISMIC` → (post-event SLCs found, topsApp run) → `READY_FOR_EMAIL` → result emailed and track removed, or `FAILED_NEEDS_ATTENTION` if processing failed.

`active_jobs/<event_id>.json` holds one event; `active_jobs/partials/` holds its HyP3-style job files waiting for post-event granules; `active_jobs/outbox/` holds the precise runner's notifications until the Action emails them.

## Urgent-response and production products
Precise orbits (AUX_POEORB) appear about 3 weeks after an acquisition, so the urgent runner always works with restituted orbits. Each event therefore gets two records:

```
<archive>/URGENT_RESPONSE/<year>/<title>_<event id>/S1-COSEIS_SAR_UR-<...>/S1-COSEIS_SAR_UR-<...>.{nc,png}
<archive>/HISTORIC_EVENTS/<year>/<title>_<event id>/S1-COSEIS_SAR-<...>/S1-COSEIS_SAR-<...>.{nc,png}
<archive>/MANIFESTS/precise_queue/S1-COSEIS_SAR_UR-<...>.json
```
- **Urgent response:** what forward mode delivered, as it delivered it: its pair, its frames, restituted orbits. Every forward-processed event gets one (forward significance criteria). The folder uses the title the event had when forward mode found it.
- **Production:** made the way historic mode makes it, for events that meet the historic significance criteria with their current USGS parameters: historic mode's pair on that track (its AOI, coverage rule and SLC search), so the dates or frames can differ from the urgent product's, with precise orbits. The folder uses the event's current USGS title, as in the rest of `HISTORIC_EVENTS/`.
- **Shared pairs** (foreshocks, aftershocks, sequences) follow the archive's rules. Before running topsApp, the urgent runner looks for an urgent product made from exactly the same SLCs; if there is one, the event is linked to it and topsApp is not run again. Each urgent product sits in the folder of the largest event that uses it (ties: the later one); it moves there if a larger event uses it later. A production product belongs to the largest significant event its pair contains, and a production pair is run once, however many urgent products lead to it.
- Both runners place products compressed, with the code that built the compressed archive (`src/aria_coseis/sar/repack.py`, run in the topsApp environment): the extra layers DockerizedTopsApp writes uncompressed get gzip, everything else is copied as is, and the result is checked against topsApp's file before it is renamed into place. Nothing in the archive is ever overwritten. topsApp's own folder (under `COSEIS_DATA_DIR`) keeps its copy of the product; SLCs and intermediate files are deleted.
- **Queue entries** (`MANIFESTS/precise_queue/`, one per urgent product) list the events that used the product, the urgent job and orbit files, the production owner and job once chosen, and a status with a log:
  - `pending`: waiting. Each precise run finds the owner, selects historic mode's pair for it (kept once found), asks ASF's orbit service whether precise orbits exist for every SLC of that pair, and then runs topsApp, checks that it used precise orbits, and places the product.
  - `done`: the production product is in `HISTORIC_EVENTS/` (or was made through another entry; or the urgent product was the same pair with precise orbits and was copied).
  - `not_eligible`: no event the pair contains meets the historic criteria; no production product.
  - `failed`: something needs a person (rerun error, a product already in the archive, or no historic pair or precise orbits after 30 days); it is emailed. To retry, fix the cause and set `"status"` back to `"pending"`. A failed rerun keeps its topsApp folder (`<title>/<track>/precise/...` under `COSEIS_DATA_DIR`).
- Both runners rewrite queue entries; they take turns through a lock file on local disk (`/tmp/coseis_precise_queue.lock`, `COSEIS_QUEUE_LOCK_FILE`).
- **Notifications:** the processing machine can't send email. The precise runner writes one JSON per placed or failed product to `outputs/outbox_pending/`; its wrapper moves them to `active_jobs/outbox/` and pushes; the next Action run emails them together to the secondary recipients and deletes them. Urgent products are reported as before, by the tracker's processing-results email.

## Processing machine setup
```bash
git clone git@github.com:cmspeed/coseis-sar.git && cd coseis-sar   # SSH: cron pushes with a key
mamba env create -f environment.yml          # or, in an existing env: pip install -e . --no-deps
crontab -e
```
Crontab entries (urgent processing at minutes 5 and 15; precise reruns every 6 hours):
```
15,05 * * * * COSEIS_DATA_DIR=<repo>/scripts/data COSEIS_ARCHIVE_DIR=<archive> <repo>/ops/run_coseis_forward.sh
30 */6 * * * COSEIS_DATA_DIR=<repo>/scripts/data COSEIS_ARCHIVE_DIR=<archive> <repo>/ops/run_coseis_precise.sh
```
- `COSEIS_ARCHIVE_DIR` is the product archive (the one-stop shop). If the archive is renamed or moved, only these two crontab lines change. Without the variable, products go to `<repo>/data/COSEIS-ONE_STOP_SHOP/`.
- `COSEIS_DATA_DIR` sets where topsApp works (SLCs, intermediate files and its copy of the product). The current machine keeps this in `scripts/data/` from before the 2026 layout change; without the variable it is `<repo>/data/`.
- Cores: urgent processing runs topsApp on cores 0-16 and precise reruns on 17-33 (`COSEIS_TOPSAPP_CORES`, `COSEIS_PRECISE_CORES`).
- The forward wrapper loads conda through `~/.bashrc`, activates `coseis-sar`, pulls `main`, runs processing, and pushes tracker changes. Its log is `logs/forward.log`. The precise wrapper does the same without pulling (the forward wrapper keeps the clone up to date); its log is `logs/precise.log`.
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
- **Pause:** comment out the crontab lines, and disable the workflow under Actions → *COSEIS Forward Email Notifications* → *Disable workflow*.
- **Resume:** reverse both. Trigger the workflow once by hand (*Run workflow*) to check it.
- **Check health:**
  - `tail logs/forward.log logs/precise.log`
  - the queue: `grep -h '"status"' <archive>/MANIFESTS/precise_queue/*.json | sort | uniq -c`
  - recent runs on the Actions tab
  - recent `[skip ci]` commits on `main`
- **Stuck runs:** the lock files `/tmp/coseis_processing.lock` and `/tmp/coseis_precise.lock` (override with `COSEIS_LOCK_FILE`, `COSEIS_PRECISE_LOCK_FILE`) prevent overlapping runs. Each holds its run's process ID; a lock whose process is gone is removed automatically.
- **Don't hand-edit `active_jobs/` while the runners are active.** Both commit to it.
- **Push failures:** git's own output, including errors, goes to `logs/forward.log` just above the "Push attempt N failed" lines. "Authentication failed", "could not read Username" or "Permission denied (publickey)" means a credentials problem (see setup above), not a conflict.
- **Push conflicts:** both runners retry their push, rebasing onto each other's commits. If the processing machine can't push because of a real conflict (both edited the same tracker lines), `logs/forward.log` shows "WARNING: local tracker commits are not on GitHub yet". In that case, resolve it by hand in the production clone (`git pull --rebase origin main`, fix the JSON, `git rebase --continue`, `git push`) while the cron is paused. A failed push in the Action fails that workflow run.

## Overrides
Every path is relative to the repository root and can be overridden with an environment variable:
- `COSEIS_DATA_DIR`
- `COSEIS_ARCHIVE_DIR`
- `COSEIS_TRACKING_DIR`
- `COSEIS_PARTIALS_DIR`
- `COSEIS_MAPS_DIR`
- `COSEIS_OUTPUT_DIR`
- `COSEIS_LOCK_FILE`, `COSEIS_PRECISE_LOCK_FILE`
- `COSEIS_OUTBOX_DIR`, `COSEIS_OUTBOX_STAGING_DIR`

These are useful for a dry run against a scratch copy of the tracker, without touching production.
