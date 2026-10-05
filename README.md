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
mamba env create -f environment.yml    # creates the `coseis-sar` environment
mamba activate coseis-sar
```
- SAR-only use needs only the core packages. The optical packages (`earthengine-api`, `google-cloud-storage`, `pystac-client`, `gdal`) are imported only when an optical backend runs.
- `scripts/batch_autorift.py` additionally needs `hyp3_autorift`.

## Usage
Run from the `scripts/` directory. Outputs (`data/`, `active_jobs/`, job lists, maps) are written relative to it.

```bash
cd scripts

# Historic SAR: HyP3 job list for every significant event in a date range
python coseis.py --historic --dates 2014-10-01 2026-07-31 --pairing coseismic --job_list

# Historic optical: Sentinel-2 composites from Google Earth Engine, then autoRIFT
python coseis.py --historic --dates 2023-02-06 --sensor sentinel-2 --optical_backend gee --optical_level toa
python batch_autorift.py [--filter] [--data_dir DIR]

# A custom list of events instead of a date range
python coseis.py --eq_list events.json --sensor landsat --optical_backend gee

# Forward mode (see Operations)
python coseis.py --forward --pairing coseismic --send_email
```

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
Two scheduled runners share the tracker files in `scripts/active_jobs/` on the `main` branch:

- **GitHub Actions** (`.github/workflows/coseis-cron.yml`) checks USGS for new events and starts tracking them. It publishes overpass maps to `docs/maps/` (GitHub Pages) and emails alerts and processing results.
- **The local cron job** on the processing machine (`scripts/run_coseis_forward.sh`) checks ASF for post-event SLCs and runs `topsApp`.

Email settings come from environment variables: `GMAIL_USER`, `GMAIL_APP_PSWD`, `COSEIS_PRIMARY_RECIPIENTS` (new events) and `COSEIS_SECONDARY_RECIPIENTS` (processing results).

`COSEIS_DATA_DIR`, `COSEIS_TRACKING_DIR` and `COSEIS_LOCK_FILE` override the output directory, the tracker directory and the lock file.

## Code layout
```
src/aria_coseis/   package: config, usgs, significance, aoi, notify,
                   sar/ (search, pairing, topsapp), optical/ (copernicus, element84, gee, jobs),
                   tracking, pipeline, modes, cli
scripts/coseis.py  command-line entry point (adds src/ to the path and calls aria_coseis.cli)
tests/             pytest suite; HTTP is replayed from recorded cassettes
```

## Development
```bash
pip install -r requirements-test.txt
pytest                                   # offline: replays tests/cassettes/
ruff check src tests scripts/coseis.py && ruff format --check src tests scripts/coseis.py
```
- To record a cassette for a new test, run `pytest --record-mode=new_episodes`.
- To refresh golden files after an intended output change, run `UPDATE_GOLDEN=1 pytest` and review the diff.
- The development workflow and roadmap (branches, phases, backlog) are in `IMPROVEMENT_PLAN.md`.

## Contact
For questions or issues, open an issue in this repository or contact [cole.speed@jpl.nasa.gov](mailto:cole.speed@jpl.nasa.gov).
