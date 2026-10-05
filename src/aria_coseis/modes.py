"""Historic and forward processing modes."""

import os
from pathlib import Path
import requests
import json
import geojson
from datetime import datetime, timedelta, timezone

from aria_coseis import config
from aria_coseis.aoi import make_aoi
from aria_coseis.config import GITHUB_PAGES_BASE_URL, USGS_api_alltime
from aria_coseis.notify import ascii_table_to_html, get_next_pass, send_email
from aria_coseis.pipeline import process_earthquake
from aria_coseis.sar.search import get_path_and_frame_numbers
from aria_coseis.significance import check_significance
from aria_coseis.tracking import add_to_tracker, check_tracker_for_updates, load_tracker
from aria_coseis.usgs import (
    get_historic_earthquake_data_date_range,
    get_historic_earthquake_data_single_date,
    parse_custom_eq_list,
    parse_geojson,
)
from aria_coseis.utils import convert_time, to_snake_case


def main_historic(
    start_date=None,
    end_date=None,
    eq_list_path=None,
    aoi=None,
    pairing_mode=None,
    job_list=False,
    resolution=90,
    sensor="sar",
    optical_backend="copernicus",
    optical_level="toa",
):
    """
    Runs the main query and processing workflow in historic processing mode.
    Used to produce 'pre-seismic', 'co-seismic', and 'post-seismic' displacement products for historic earthquakes.
    Pre-seismic and post-seismic data are generated for 90 days before and 30 days after the event.
    A single date or date range can be provided for processing.
    :param start_date: The query start date in YYYY-MM-DD format
    :param end_date: The query end date in YYYY-MM-DD format (Optional)
    :param eq_list_path: The path to a custom JSON file containing a list of earthquakes to process (Optional)
    :param aoi: The path to a JSON file representing the area of interest (AOI) (Optional)
    :param pairing_mode: 'all', 'sequential', or 'coseismic' for specifying desired SLC pairing.
    :param job_list: If True, create a list of jobs in HYP3 format for cloud processing.
    :param resolution: Output resolution for the topsApp processing, default is 90m
    :param sensor: 'sar' for SAR processing, 'optical' for optical processing
    :optical_backend: 'copernicus', 'element84' for source data file nomenclature (only applicable if sensor is 'optical')
    :optical_level: 'toa' or 'sr' for the desired optical product level (only applicable if sensor is 'optical')
    """
    # Generate the list of earthquakes
    geojson_data = None

    # If a custom earthquake list is provided, use that instead of querying the USGS API, else use the provided dates to query the USGS API for earthquakes in that time range
    if eq_list_path:
        eq_sig = parse_custom_eq_list(eq_list_path)
        if not eq_sig:
            print("No valid earthquakes found in the provided list.")
            return
    else:
        if start_date and not end_date:
            print("=========================================")
            print(f"Running historic processing in single-date mode for date: {start_date}")
            print("=========================================")
            geojson_data = get_historic_earthquake_data_single_date(
                USGS_api_alltime, str(start_date)
            )

        elif start_date and end_date:
            print("=========================================")
            print(
                f"Running historic processing in date range mode for dates: {start_date} to {end_date}"
            )
            print("=========================================")
            geojson_data = get_historic_earthquake_data_date_range(
                USGS_api_alltime, str(start_date), str(end_date)
            )

        if geojson_data:
            earthquakes = parse_geojson(geojson_data)
            eq_sig = check_significance(
                earthquakes, start_date, end_date, sensor=sensor, mode="historic"
            )
        else:
            eq_sig = None

    # Process the list of earthquakes
    if eq_sig is not None:
        jobs_dict = []
        master_scene_features = []
        earthquake_infos = []

        for eq in eq_sig:
            try:
                eq_jsons, eq_features = process_earthquake(
                    eq,
                    aoi,
                    pairing_mode,
                    job_list,
                    resolution,
                    sensor,
                    optical_backend,
                    optical_level,
                )

                if eq_features:
                    master_scene_features.extend(eq_features)

                if eq_jsons:
                    event_dt = convert_time(eq["time"])
                    eq_info = {
                        "title": eq.get("title"),
                        "epicenter": {
                            "latitude": eq["coordinates"][1],
                            "longitude": eq["coordinates"][0],
                            "depth_km": eq["coordinates"][2],
                            "usgs_event_url": eq.get("url", ""),
                        },
                        "time": event_dt.strftime("%Y-%m-%d %H:%M:%S UTC"),
                    }
                    earthquake_infos.append(eq_info)

            except Exception as e:
                print(f"Error processing {eq['title']}: {e}")
                continue

            if eq_jsons:
                for i, eq_json in enumerate(eq_jsons):
                    for j, json_data in enumerate(eq_json):
                        jobs_dict.append(json_data)

        if jobs_dict:
            current_time = datetime.now(timezone.utc).strftime("%Y-%m-%d_%H-%M-%S_UTC")

            with open(f"jobs_list_{current_time}.json", "w") as f:
                json.dump(jobs_dict, f, indent=4)

            with open(f"earthquake_info_{current_time}.json", "w", encoding="utf-8") as f:
                json.dump(earthquake_infos, f, indent=4, ensure_ascii=False)

        if master_scene_features:
            current_time = datetime.now(timezone.utc).strftime("%Y-%m-%d_%H-%M-%S")
            feature_filename = (
                f"all_selected_scenes_{sensor}_{optical_backend}_{current_time}.geojson"
            )
            fc = {"type": "FeatureCollection", "features": master_scene_features}
            with open(feature_filename, "w") as f:
                json.dump(fc, f, indent=2)
            print(f"Saved master scene footprints to {feature_filename}")
    else:
        if eq_list_path:
            print("No significant earthquakes found in the provided list.")
        else:
            print(f"No significant earthquakes found between {start_date} and {end_date}.")


def main_forward(
    pairing_mode=None,
    resolution=30,
    do_processing=False,
    send_email_flag=False,
    process_only=False,
):
    """
    Runs the main query and processing workflow in forward processing mode.
    Used to produce co-seismic product for new earthquakes when new SLC data becomes available.
    :param pairing_mode: 'all', 'sequential', or 'coseismic' for specifying desired SLC pairing
    :param resolution: Output resolution for the topsApp processing, default is 30m
    :param do_processing: If True, runs the dockerized topsApp processing workflow after generating the JSONs. Default is False.
    :param send_email_flag: If True, sends an email alert after processing. Default is False.
    :param process_only: If True, only runs the processing workflow without generating new JSONs or sending emails. Default is False.
    """
    import shutil

    # A lock file to prevent overlapping runs
    lock_file = config.LOCK_FILE

    if os.path.exists(lock_file):
        print("Previous processing run still active. Exiting.")
        return
    try:
        with open(lock_file, "w") as f:
            f.write("running")

        print("=========================================")
        print("Running cronjob to check for new earthquakes...")
        print("=========================================")

        # Initialize the tracking directory if it doesn't exist
        if not os.path.exists(config.TRACKING_DIR):
            os.makedirs(config.TRACKING_DIR, exist_ok=True)

        if not process_only:
            # Check for New Earthquakes over 48-hour window to ensre no events are missed due to API delays
            two_days_ago = (datetime.now(timezone.utc) - timedelta(days=2)).strftime(
                "%Y-%m-%dT%H:%M:%S"
            )

            # Define parameters for a custom search on the USGS 'alltime' endpoint
            params = {"format": "geojson", "starttime": two_days_ago, "minmagnitude": 5.5}
            print(f"Checking for earthquakes since {two_days_ago}...")

            # Use the query endpoint instead of the static summary feeds
            response = requests.get(USGS_api_alltime, params=params)
            response.raise_for_status()
            geojson_data = response.json()

            start_date = datetime.now().strftime("%Y-%m-%d")
            current_time = datetime.now(timezone.utc).strftime("%Y-%m-%d at %H:%M:%S UTC")

            if geojson_data:
                # Parse GeoJSON and create variables for each feature's properties
                earthquakes = parse_geojson(geojson_data)
                eq_sig = check_significance(earthquakes, start_date, end_date=None, mode="forward")

                if eq_sig is not None:
                    for eq in eq_sig:
                        # Check for duplicate entry in the pending queue
                        tracker = load_tracker()
                        if eq.get("id") in tracker:
                            print(
                                f"Earthquake with ID {eq.get('id')} is already in the pending queue. Skipping."
                            )
                            continue

                        title = eq.get("title", "")
                        title_snake = to_snake_case(title)
                        print(f"title: {title_snake}")
                        coords = eq.get("coordinates", [])

                        # Initial AOI creation (a 1-degree box)
                        aoi = make_aoi(coords)

                        # Write AOI to a geojson file
                        with open(f"{title}_AOI.geojson", "w") as f:
                            geojson.dump(aoi, f, indent=2)

                        # Get path/frame numbers for the initial AOI
                        path_frame_numbers, frame_dataframe = get_path_and_frame_numbers(
                            aoi, eq.get("time")
                        )

                        # Create a timestamp string
                        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

                        # Create the output directory
                        timestamp_dir = Path(f"nextpass_outputs_{timestamp}")
                        timestamp_dir.mkdir(parents=True, exist_ok=True)

                        # Run next_pass to get the next overpasses
                        s1_info, nisar_info, overpass_map = get_next_pass(aoi, timestamp_dir)

                        # Setup Github pages directory
                        docs_maps_dir = Path(os.getcwd()).parent / "docs" / "maps"
                        docs_maps_dir.mkdir(parents=True, exist_ok=True)

                        # Create a unique ID for the filenames so they aren't overwritten
                        unique_id = eq.get("id", datetime.now().strftime("%Y%m%d%H%M%S"))

                        # Route Joint Map
                        map_url = ""
                        if overpass_map and os.path.exists(overpass_map):
                            new_map_name = f"{title_snake}_{unique_id}_overpass_map.html"
                            shutil.copy(overpass_map, docs_maps_dir / new_map_name)
                            map_url = f"{GITHUB_PAGES_BASE_URL}/maps/{new_map_name}"

                        # Construct and send the initial email alert
                        message_dict = {
                            "title": eq.get("title", ""),
                            "time": convert_time(eq["time"]).strftime("%Y-%m-%d %H:%M:%S"),
                            "coordinates": [
                                round(coord, 3) for coord in eq.get("coordinates", [])
                            ],
                            "magnitude": eq.get("mag", ""),
                            "depth": round(eq.get("coordinates", [])[2], 1),
                            "alert": eq.get("alert", ""),
                            "url": eq.get("url", ""),
                        }

                        # Convert the raw text table to an HTML table
                        html_s1_table = ascii_table_to_html(s1_info)
                        html_nisar_table = ascii_table_to_html(nisar_info)

                        # Construct and send the email
                        header_html = f"""
                        <div style="font-family: Arial, sans-serif; color: #333; max-width: 850px; margin: auto; border: 1px solid #e0e0e0; border-radius: 8px; overflow: hidden;">
                            <div style="background-color: #003366; color: white; padding: 20px;">
                                <h2 style="margin: 0; font-size: 22px;">{message_dict["title"]}</h2>
                                <p style="margin: 5px 0 0; font-size: 14px; color: #b3d4fc;">{message_dict["time"]} UTC</p>
                            </div>
                            <div style="padding: 20px;">
                                <h3 style="margin: 0 0 10px 0; border-bottom: 2px solid #f0f0f0; padding-bottom: 8px; color: #003366;">Event Details</h3>
                                <table style="width: 100%; text-align: left; margin-bottom: 25px; border-collapse: collapse;">
                                    <tr>
                                        <th style="width: 150px; padding: 4px 0;">Epicenter (Lat, Lon):</th>
                                        <td style="padding: 4px 0;">{message_dict["coordinates"][1]}, {message_dict["coordinates"][0]}</td>
                                    </tr>
                                    <tr>
                                        <th style="padding: 4px 0;">Depth:</th>
                                        <td style="padding: 4px 0;">{message_dict["depth"]} km</td>
                                    </tr>
                                </table>
                                <a href="{message_dict["url"]}" style="display: inline-block; padding: 10px 18px; background-color: #0055a4; color: white; text-decoration: none; border-radius: 5px; font-weight: bold; margin-bottom: 30px;">View on USGS Hazard Portal</a>
                        """

                        s1_section = f"""
                                <h3 style="margin: 0 0 10px 0; border-bottom: 2px solid #f0f0f0; padding-bottom: 8px; color: #003366;">Sentinel-1 Acquisitions</h3>
                                {html_s1_table}
                        """

                        nisar_section = f"""
                                <h3 style="margin: 25px 0 10px 0; border-bottom: 2px solid #f0f0f0; padding-bottom: 8px; color: #003366;">NISAR Acquisitions</h3>
                                {html_nisar_table}
                        """

                        footer_html = """
                            </div>
                            <div style="background-color: #f9f9f9; padding: 15px; text-align: center; font-size: 12px; color: #888; border-top: 1px solid #e0e0e0;">
                                This is an automated message. Please do not reply.<br>
                                For product-specific inquiries, contact Dr. Cole Speed (<a href="mailto:cole.speed@jpl.nasa.gov">cole.speed@jpl.nasa.gov</a>) and Dr. Grace Bato (<a href="mailto:bato@jpl.nasa.gov">bato@jpl.nasa.gov</a>).
                            </div>
                        </div>
                        """

                        # Helper function to generate the clickable button to route to HTML map
                        def get_button_html(url):
                            if not url:
                                return ""
                            return f"""
                            <div style="text-align: center; margin: 30px 0;">
                                <a href="{url}" style="background-color: #003366; color: white; padding: 15px 40px; text-decoration: none; border-radius: 50px; display: inline-block; font-family: Arial, sans-serif; font-size: 18px;">
                                    Click here for interactive overpass map
                                </a>
                            </div>
                            """

                        if send_email_flag:
                            subject_text = f"New Event: {message_dict['title']}"

                            # Send joint S1 + NISAR email to PRIMARY_RECIPIENTS
                            if config.PRIMARY_RECIPIENTS:
                                primary_body = (
                                    header_html
                                    + s1_section
                                    + nisar_section
                                    + get_button_html(map_url)
                                    + footer_html
                                ).replace("\n", "")
                                send_email(
                                    subject=subject_text,
                                    body=primary_body,
                                    recipients=config.PRIMARY_RECIPIENTS,
                                )
                                print("=========================================")
                                print("Joint S1 and NISAR email sent to primary recipients.")
                                print("=========================================")
                        else:
                            print("=========================================")
                            print("Email sending is disabled (--send_email not provided).")
                            print("=========================================")

                        # START TRACKING FOR THIS EVENT
                        # Finds pre-seismic SLCs, creates partial job list, and saves to tracking file
                        add_to_tracker(eq, aoi, resolution)

                else:
                    print(f"No new significant earthquakes found as of {current_time}.")
        else:
            print("Running in --process_only mode. Skipping USGS earthquake discovery.")

        # Check ASF DAAC for available SLCs for pending earthquakes
        check_tracker_for_updates(do_processing, send_email_flag)

    finally:
        if os.path.exists(lock_file):
            os.remove(lock_file)
