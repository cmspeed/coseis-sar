import re
import unicodedata
import os
import argparse
import asf_search as asf
from dateutil import parser as dateparser
from dateutil.parser import isoparse
from pathlib import Path
import requests
import json
import folium
import geojson
import geopandas as gpd
import glob
import csv
import math
from shapely import wkt
from shapely.geometry import mapping, shape, box, Point, Polygon, LineString, MultiLineString, MultiPolygon
from shapely.ops import unary_union, linemerge
import subprocess
import sys
from datetime import datetime, timedelta, timezone
from collections import defaultdict
from itertools import combinations
import logging
import yagmail
import time
from time import sleep
from types import SimpleNamespace
from urllib.parse import urlparse
from typing import List, Dict, Any, Optional

from aria_coseis import config
from aria_coseis.config import (
    ASF_DAAC_API,
    OPTICAL_CLOUD_THRESHOLD,
    USGS_api_alltime,
    coastline_api,
)
from aria_coseis.utils import convert_time, to_snake_case
from aria_coseis.notify import ascii_table_to_html, get_next_pass, make_interactive_map, send_email
from aria_coseis.usgs import (
    get_event_rake,
    get_ffm_geojson_url,
    get_historic_earthquake_data_date_range,
    get_historic_earthquake_data_single_date,
    parse_custom_eq_list,
    parse_geojson,
)
from aria_coseis.aoi import load_aoi_from_json, make_aoi
from aria_coseis.significance import check_significance
from aria_coseis.sar.search import get_SLCs, get_path_and_frame_numbers
from aria_coseis.sar.pairing import find_reference_and_secondary_pairs, make_job_json
from aria_coseis.sar.topsapp import run_dockerized_topsApp
from aria_coseis.optical.jobs import get_utm_zone, make_optical_job_json
from aria_coseis.optical.copernicus import find_optical_pairs_copernicus, search_copernicus_public
from aria_coseis.optical.element84 import find_optical_pairs_element84, search_element84_stac
from aria_coseis.optical.gee import (
    assign_nodata,
    download_from_gcs,
    export_gee_landsat_composite,
    export_gee_sentinel2_composite,
    merge_and_compress_chips,
    wait_for_gee_tasks,
)
from aria_coseis.tracking import add_to_tracker, check_tracker_for_updates, load_tracker
from aria_coseis.pipeline import process_earthquake
