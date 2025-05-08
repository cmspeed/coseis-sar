# coseis-sar

**coseis-sar** is a Python package designed to generate coseismic Sentinel-1 Single Look Complex synthetic aperture radar (SAR) image pairs for interferometric synthetic aperture radar (InSAR) processing using ISCE2's topsApp. It supports both local processing and cloud-based processing via HYP3.

## Features
- Generates Sentinel-1 SLC pairs for InSAR processing based on the user-specified `pairing` mode.
- Supports **historic** and **forward** processing modes
- Produces a `job_list.json` file for cloud processing via HYP3 (if specified).
- Runs ISCE2 topsApp locally when `job_list` is not specified.
- Generates an AOI around an earthquake epicenter dynamically and and automatically identifies intersecting Sentinel-1 SLCs intersecting the AOI and temporally bounding the earthquake event.

## Installation
Clone the repository and install any required dependencies:

```bash
git clone https://github.com/cmspeed/coseis-sar.git
cd coseis-sar
pip install -r requirements.txt
```

## Usage

```bash
python coseis_sar.py --historic --dates <date1> <date2> --pairing <pairing_mode> [--job_list]
python coseis_sar.py --forward --pairing <pairing_mode>
```

#### Arguments
- `--historic` : Used for processing past earthquakes or generating a job list for cloud processing.
- `--forward` : Used for detecting Sentinel-1 data for recent earthquakes (within the last hour), typically run as a cron job.
- `--dates <date1> <date2>` : Specifies the date range for historic processing.
- `--pairing <mode>` : Defines how SLC pairs are selected. Options:
  - `coseismic` : Pairs spanning the earthquake event.
  - `sequential` : Consecutive SLC acquisitions.
  - `all` : All possible pairs.
- `--job_list` (optional) : If specified, creates `job_list.json` for HYP3 cloud processing instead of running topsApp locally.

## Examples

### Historic Mode (Local Processing of 2025 Southern Tibetan Plateau Earthquake)
```bash
python coseis_sar.py --historic --dates 2025-01-07 2025-01-07 --pairing coseismic
```

This will generate `coseis-sar` products for all intersecting Sentinel-1 tracks (ascending and descending) that intersect the automatically defined AOI, which is based on the epicenter coordinates provided by the USGS Earthquake Portal API.

### Historic Mode (HYP3 Cloud Processing)
```bash
python coseis_sar.py --historic --dates 2014-01-01 2025-12-31 --pairing coseismic --job_list
```

This will generate a `job_list.json` with all of the required parameters to deploy InSAR processing in the cloud with HYP3 and generate `coseis-sar` products for each earthquake coseismic pair defined in the `job_list.json`. 


### Forward Mode (Real-time Processing)
```bash
python coseis_sar.py --forward --pairing coseismic
```

This will initiate a workflow that requests all earthquakes that have occurred in the last hour from the USGS Earthquake Portal API and determines the intersecting SLCs needed for InSAR processing. Information about the earthquake, as well as information about the corresponding Sentinel-1 data, is distributed automatically via email to a user-defined list of recipients. Currently, `forward` mode is largely applied via a cronjob to notify researchers of new earthquakes and prepare for subsequent InSAR processing workflows.


## HYP3 Processing
For cloud-based processing, `job_list.json` will be generated when `--job_list` is specified. This file contains the parameters necessary for running ISCE2 topsApp via HYP3 for automated processing. [More details on HYP3 processing here.](http://hyp3-docs.asf.alaska.edu)

## Contact
For questions or issues, open an issue in the GitHub repository or contact [cole.speed@jpl.nasa.gov](mailto:cole.speed@jpl.nasa.gov).
