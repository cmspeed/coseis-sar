"""Settings: API endpoints, thresholds, working paths and email recipients."""
import logging
import os

# Set logging level to WARNING to suppress DEBUG and INFO logs
logging.basicConfig(level=logging.WARNING)

# API endpoints
USGS_api_hourly = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_hour.geojson"  # USGS Earthquake API - Hourly
USGS_api_daily = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_day.geojson"  # USGS Earthquake API - Daily
USGS_api_30day = "https://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/all_month.geojson"  # USGS Earthquake API - Monthly
USGS_api_alltime = "https://earthquake.usgs.gov/fdsnws/event/1/query" # USGS Earthquake API - All Time
coastline_api = "https://raw.githubusercontent.com/OSGeo/PROJ/refs/heads/master/docs/plot/data/coastline.geojson" # Coastline API
ASF_DAAC_API = "https://api.daac.asf.alaska.edu/services/search/param" # ASF DAAC API endpoint
CMR_API_URL = "https://cmr.earthdata.nasa.gov/search/granules.json" # NASA CMR API endpoint

# Published overpass maps (docs/maps/ on GitHub Pages). Pages URLs do not redirect if the repo is renamed.
GITHUB_PAGES_BASE_URL = "https://cmspeed.github.io/coseis-sar"
# Processing outputs. Override with COSEIS_DATA_DIR (e.g. for a shadow run on the processing machine).
root_dir = os.environ.get("COSEIS_DATA_DIR") or os.path.join(os.getcwd(), "data")

# Global variables
OPTICAL_CLOUD_THRESHOLD = 20.0  # Maximum cloud cover percentage for optical data

# Forward-mode tracker directory (relative to the working directory unless absolute)
TRACKING_DIR = os.environ.get("COSEIS_TRACKING_DIR") or "active_jobs"

# Prevents overlapping forward runs; run_coseis_forward.sh checks the default path too
LOCK_FILE = os.environ.get("COSEIS_LOCK_FILE") or "/tmp/coseis_processing.lock"

def get_recipients_from_env(var_name):
    """
    Retrieves a list of emails from an environment variable.
    """
    env_val = os.getenv(var_name, "")
    # Split by comma and strip whitespace
    return [email.strip() for email in env_val.split(',') if email.strip()]

# Load recipients from environment variables
PRIMARY_RECIPIENTS = get_recipients_from_env('COSEIS_PRIMARY_RECIPIENTS')
SECONDARY_RECIPIENTS = get_recipients_from_env('COSEIS_SECONDARY_RECIPIENTS')
