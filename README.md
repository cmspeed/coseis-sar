# ARIA COSEIS

Automated coseismic surface displacement for significant earthquakes worldwide.

ARIA COSEIS selects significant earthquakes from the USGS catalog and prepares displacement processing from satellite imagery that brackets each event:

- **SAR (Sentinel-1):** finds coseismic SLC pairs for every intersecting track and writes jobs for ISCE2 `topsApp` (interferometry and pixel offsets) on [HyP3](https://hyp3-docs.asf.alaska.edu), giving line-of-sight displacement.
- **Optical (Sentinel-2 / Landsat):** builds cloud-free pre- and post-event composites in Google Earth Engine and runs `autoRIFT` on them, giving 2-D horizontal displacement.

The examples below cover historic mode (past events). A separate forward mode runs automatically for new earthquakes; see [Forward mode](#forward-mode).

## Installation
```bash
git clone https://github.com/cmspeed/coseis-sar.git
cd coseis-sar
mamba env create -f environment.yml     # creates the `coseis-sar` env and installs this package
mamba activate coseis-sar
```
Without conda, `pip install -e .` installs the package and its core dependencies. Optional extras:
- `.[optical]`: Earth Engine, Cloud Storage, STAC. GDAL is best installed from conda.
- `.[forward]`: `next_pass`, used for overpass predictions.
- `.[test]`: test and lint tools.

Optical autoRIFT processing also needs `hyp3_autorift`.

## Quick start
Run commands from the repository root. `aria-coseis` and `python -m aria_coseis` are the same command.

### SAR: HyP3 job list for past earthquakes
```bash
# One day
python -m aria_coseis --historic --dates 2025-01-07 --pairing coseismic --job_list

# A date range (e.g. the whole Sentinel-1 era)
python -m aria_coseis --historic --dates 2014-10-01 2026-07-31 --pairing coseismic --job_list
```
Writes `outputs/jobs_list_<timestamp>.json`, ready to submit to HyP3, plus a summary of each event (`outputs/earthquake_info_<timestamp>.json`). Without `--job_list`, it writes per-pair topsApp parameters and a map of the intersecting frames instead.

### Optical: Sentinel-2 or Landsat with autoRIFT
```bash
# 1. Build pre/post-event composites in Earth Engine and download them to data/
python -m aria_coseis --historic --dates 2023-02-06 --sensor sentinel-2 --optical_backend gee --optical_level toa

# 2. Run autoRIFT on everything downloaded so far
python -m aria_coseis.optical.autorift --filter
```
- **Requirements:**
  - an Earth Engine account authenticated with `earthengine authenticate`, with access to the Cloud project `coseis-1`
  - Google Cloud credentials for downloading (e.g. `gcloud auth application-default login`)
  - the export bucket set in `COSEIS_GCS_BUCKET`
- **Event selection:** only strike-slip events are processed for optical.
- **Results:** composites, manifests and autoRIFT offsets (in meters) go to `data/GEE_Optical_Downloads/<event>/`.
- **Other backends:** `--optical_backend copernicus` or `element84` only search and pair Sentinel-2 scenes (they don't build composites) and write autoRIFT job JSON. Copernicus writes it only with `--job_list`.

### A custom list of events
```bash
python -m aria_coseis --eq_list events.json --pairing coseismic --job_list
```
`events.json` uses the same format as the `earthquake_info_*.json` files that historic runs write:
```json
[{"title": "M 7.1 - 2025 Southern Tibetan Plateau Earthquake",
  "time": "2025-01-07 01:05:16 UTC",
  "epicenter": {"latitude": 28.639, "longitude": 87.361, "depth_km": 10.0,
                "usgs_event_url": "https://earthquake.usgs.gov/earthquakes/eventpage/us6000pi9w"}}]
```
Listed events are processed as given; the significance criteria below are not applied.

## How events and areas are chosen
- **Significant earthquakes (historic mode):** M ≥ 6.0, depth ≤ 40 km, epicenter on land or near the coast (open-ocean events are excluded). Optical processing additionally requires strike-slip faulting (USGS rake within 45° of 0° or 180°).
- **Area of interest:** the USGS finite-fault model when one exists (buffered by 0.15° for optical), otherwise a 1° × 1° box centered on the epicenter. `--aoi FILE` supplies your own GeoJSON instead.
- **SAR pairs:** per track, the nearest acquisitions before and after the event that cover the area (searched from 90 days before to 30 days after).

## Options
| Option | Meaning |
|---|---|
| `--historic` / `--eq_list FILE` / `--forward` | Run mode (choose one) |
| `--dates START [END]` | Day or date range, `YYYY-MM-DD` (historic) |
| `--pairing {coseismic,sequential,all}` | SAR pair selection; required for SAR |
| `--job_list` | Write HyP3 job JSON |
| `--resolution M` | topsApp output resolution in meters (default 30) |
| `--sensor {sar,sentinel-2,landsat}` | Imagery (default `sar`) |
| `--optical_backend {copernicus,element84,gee}` | Optical data source (default `copernicus`) |
| `--optical_level {toa,sr,raw}` | Optical product level (default `toa`; `raw` is Landsat only) |
| `--aoi FILE` | Custom area of interest (GeoJSON) |

Run `python -m aria_coseis --help` for the full list.

## Where files go
Relative to the repository root; each location can be changed with an environment variable.

| Location | Contents | Override |
|---|---|---|
| `outputs/` | job lists, event summaries, significance tables, AOIs, frame maps | `COSEIS_OUTPUT_DIR` |
| `data/` | Earth Engine downloads, autoRIFT results, topsApp products | `COSEIS_DATA_DIR` |

`outputs/` and `data/` are not tracked by git.

## Forward mode
Forward mode watches USGS for new significant earthquakes. For each one it emails an alert with predicted Sentinel-1 and NISAR overpasses, and processes the SAR pair once post-event data are available. It runs on a schedule (GitHub Actions plus a cron job on a processing machine) and keeps its state in `active_jobs/`. Operator notes are in [ops/README.md](ops/README.md).

## Development
```bash
pip install -e ".[test]"
pytest                                             # offline: replays recorded HTTP in tests/cassettes/
ruff check src tests && ruff format --check src tests
```
- Record a cassette for a new test with `pytest --record-mode=new_episodes`.
- After an intended output change, refresh golden files with `UPDATE_GOLDEN=1 pytest` and review the diff.
- Work happens on feature branches with pull requests into `main`; CI runs the tests on every PR. The roadmap and backlog are in [IMPROVEMENT_PLAN.md](IMPROVEMENT_PLAN.md).

Code layout:
```
src/aria_coseis/   config, usgs, significance, aoi, notify, tracking, pipeline, modes, cli
                   sar/      search, pairing, topsapp
                   optical/  copernicus, element84, gee, jobs, autorift
tests/             pytest suite
ops/               forward-mode cron wrapper and operator notes
active_jobs/       forward-mode tracker (written by the scheduled runners)
docs/maps/         overpass maps published on GitHub Pages (written by the runners)
```

## Acknowledgements
The 2026 restructuring of this codebase (the `aria_coseis` package, test suite and CI) was developed with assistance from [Claude Code](https://claude.com/claude-code) (Anthropic); all changes were reviewed by the maintainers.

## Contact
Open an issue in this repository or contact [cole.speed@jpl.nasa.gov](mailto:cole.speed@jpl.nasa.gov).
