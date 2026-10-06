# ARIA COSEIS

Automated coseismic surface displacement for significant earthquakes worldwide.

ARIA COSEIS finds significant shallow earthquakes in the USGS catalog. It then finds the Sentinel-1 SAR and Sentinel-2/Landsat optical imagery that brackets each event, and prepares or runs the processing:

- **SAR**: ISCE2 `topsApp` interferograms and offsets, run locally or on [HyP3](https://hyp3-docs.asf.alaska.edu), giving line-of-sight displacement.
- **Optical**: pre- and post-event composites from Google Earth Engine (GEE), processed with `autoRIFT`, giving 2-D horizontal displacement.

It runs in two modes:
- **Historic** mode processes past events.
- **Forward** mode runs on a schedule. It watches for new events, emails alerts with predicted satellite overpasses, and processes SAR pairs once post-event data arrive.

## Installation
```bash
git clone https://github.com/cmspeed/coseis-sar.git
cd coseis-sar
mamba env create -f environment.yml    # creates `coseis-sar` and installs this package (editable)
mamba activate coseis-sar
```
- Without conda, `pip install -e ".[forward]"` installs the package and its core dependencies, plus `next_pass` for forward-mode alerts.
- SAR-only use needs only the core packages. The optical packages (`earthengine-api`, `google-cloud-storage`, `pystac-client`, `gdal`) are imported only when an optical backend runs.
- `aria_coseis.optical.autorift` additionally needs `hyp3_autorift`.

## Usage
Run from the repository root. `aria-coseis` is equivalent to `python -m aria_coseis`.

```bash
# Historic SAR: HyP3 job list for every significant event in a date range
python -m aria_coseis --historic --dates 2014-10-01 2026-07-31 --pairing coseismic --job_list

# Historic optical: Sentinel-2 composites from Google Earth Engine, then autoRIFT
python -m aria_coseis --historic --dates 2023-02-06 --sensor sentinel-2 --optical_backend gee --optical_level toa
python -m aria_coseis.optical.autorift [--filter] [--data_dir DIR]

# A custom list of events instead of a date range
python -m aria_coseis --eq_list events.json --sensor landsat --optical_backend gee

# Forward mode (see Operations)
python -m aria_coseis --forward --pairing coseismic --send_email
```

Where files go (all relative to the repository root, each overridable with an environment variable):

| Location | Contents | Tracked | Override |
|---|---|---|---|
| `outputs/` | AOIs, frame maps, job lists, significance tables, next-pass results | no | `COSEIS_OUTPUT_DIR` |
| `data/` | topsApp products, Earth Engine downloads and autoRIFT manifests | no | `COSEIS_DATA_DIR` |
| `active_jobs/` | forward-mode tracker, one JSON per event; `partials/` holds jobs awaiting post-event data | yes | `COSEIS_TRACKING_DIR`, `COSEIS_PARTIALS_DIR` |
| `docs/maps/` | overpass maps published on GitHub Pages | yes | `COSEIS_MAPS_DIR` |
| `logs/` | forward-mode cron log | no | — |

| Option | Meaning |
|---|---|
| `--historic` / `--forward` / `--eq_list FILE` | Run mode (choose one) |
| `--dates START [END]` | Historic date or date range (`YYYY-MM-DD`) |
| `--pairing {coseismic,sequential,all}` | SAR pair selection (required for SAR) |
| `--job_list` | Write HyP3 job JSON instead of local pair JSON |
| `--resolution M` | topsApp output resolution in meters (default 30) |
| `--sensor {sar,sentinel-2,landsat}` | Imagery (default `sar`); forward mode is SAR only for now |
| `--optical_backend {copernicus,element84,gee}`, `--optical_level {raw,toa,sr}` | Optical data source and product level |
| `--do_processing`, `--send_email`, `--process_only` | Forward mode: run topsApp, send emails, skip discovery |
| `--aoi FILE` | Use this GeoJSON AOI instead of the automatic one |

**Significance criteria.** Historic mode keeps M ≥ 6.0 events at ≤ 40 km depth. Forward mode keeps (M ≥ 5.5 and ≤ 15 km) or (M ≥ 6.0 and ≤ 40 km). Both modes drop mid-ocean events. Optical runs keep only strike-slip events (rake within 45° of 0°/180°).

**AOI.** The AOI is the USGS finite-fault model when one exists (buffered 0.15° for optical runs), otherwise a 1° box around the epicenter.

## Operations (forward mode)
Two scheduled runners share the tracker files in `active_jobs/` on the `main` branch:

- **GitHub Actions** (`.github/workflows/coseis-cron.yml`) checks USGS for new events and starts tracking them. It publishes overpass maps to `docs/maps/` (GitHub Pages) and emails alerts and processing results.
- **The local cron job** on the processing machine (`ops/run_coseis_forward.sh`) checks ASF for post-event SLCs and runs `topsApp`.

Email settings come from environment variables: `GMAIL_USER`, `GMAIL_APP_PSWD`, `COSEIS_PRIMARY_RECIPIENTS` (new events) and `COSEIS_SECONDARY_RECIPIENTS` (processing results).

`COSEIS_LOCK_FILE` overrides the lock file that prevents overlapping forward runs (default `/tmp/coseis_processing.lock`).

## Code layout
```
src/aria_coseis/   package: config, usgs, significance, aoi, notify,
                   sar/ (search, pairing, topsapp), optical/ (copernicus, element84, gee, jobs, autorift),
                   tracking, pipeline, modes, cli (python -m aria_coseis)
tests/             pytest suite; HTTP is replayed from recorded cassettes
active_jobs/       forward-mode tracker (written by the runners)
docs/maps/         overpass maps on GitHub Pages (written by the runners)
ops/               processing-machine cron wrapper
```

## Development
```bash
pip install -e ".[test]"
pytest                                   # offline: replays tests/cassettes/
ruff check src tests && ruff format --check src tests
```
- To record a cassette for a new test, run `pytest --record-mode=new_episodes`.
- To refresh golden files after an intended output change, run `UPDATE_GOLDEN=1 pytest` and review the diff.
- The development workflow and roadmap (branches, phases, backlog) are in `IMPROVEMENT_PLAN.md`.

## Contact
For questions or issues, open an issue in this repository or contact [cole.speed@jpl.nasa.gov](mailto:cole.speed@jpl.nasa.gov).
