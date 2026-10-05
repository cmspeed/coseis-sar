"""Forward-mode job tracking: tracker files and the AWAITING -> READY_FOR_EMAIL state machine."""

import os
import json
import glob
from shapely.geometry import mapping, shape
from shapely.ops import unary_union
from datetime import datetime
from collections import defaultdict

from aria_coseis import config
from aria_coseis.notify import send_email
from aria_coseis.sar.pairing import make_job_json
from aria_coseis.sar.search import get_SLCs, get_path_and_frame_numbers
from aria_coseis.sar.topsapp import run_dockerized_topsApp
from aria_coseis.utils import convert_time, to_snake_case


def load_tracker():
    """Loads all active jobs from the tracking directory."""
    tracker = {}
    if not os.path.exists(config.TRACKING_DIR):
        os.makedirs(config.TRACKING_DIR, exist_ok=True)
        return tracker

    for file in glob.glob(os.path.join(config.TRACKING_DIR, "*.json")):
        try:
            with open(file, "r") as f:
                event_data = json.load(f)
                # The filename (minus .json) is the event_id
                event_id = os.path.basename(file).replace(".json", "")
                tracker[event_id] = event_data
        except json.JSONDecodeError:
            continue
    return tracker


def save_tracker(data):
    """Saves the tracking data back to individual files."""
    if not os.path.exists(config.TRACKING_DIR):
        os.makedirs(config.TRACKING_DIR, exist_ok=True)

    # Save current tracker state to individual files
    for event_id, event_data in data.items():
        file_path = os.path.join(config.TRACKING_DIR, f"{event_id}.json")
        with open(file_path, "w") as f:
            json.dump(event_data, f, indent=4)

    # Remove files for events that are no longer in the tracker dictionary
    for file in glob.glob(os.path.join(config.TRACKING_DIR, "*.json")):
        event_id = os.path.basename(file).replace(".json", "")
        if event_id not in data:
            os.remove(file)


def add_to_tracker(eq, aoi, resolution=30):
    """
    Initializes tracking for a new earthquake.
    Identifies intersecting tracks, finds pre-seismic SLCs for each track,
    creates a partial job file (granules=Empty, secondary_granules=Filled), and adds entry to tracking file.
    :param eq: dictionary containing earthquake data
    :param aoi: shapely Polygon object representing the area of interest
    :param resolution: desired output resolution for processing (default=30)
    """
    tracker = load_tracker()
    event_id = eq.get("id")
    title = to_snake_case(eq.get("title"))
    event_time = eq.get("time")

    if event_id in tracker:
        print(f"Event {title} is already being tracked.")
        return

    print("=========================================")
    print(f"Initializing tracking for {title}...")
    print("=========================================")

    # Get intersecting tracks
    path_frame_numbers, _ = get_path_and_frame_numbers(aoi, event_time)

    tracks_info = {}

    for (flight_direction, path_number), frame_numbers_set in path_frame_numbers.items():
        frame_numbers = list(set(fn[0] for fn in frame_numbers_set))

        # Unique key for this track
        track_key = f"{flight_direction}_{path_number}"

        # Fetch SLCs intersecting the AOI (Historical search relative to event time)
        slcs = get_SLCs(
            flight_direction, path_number, aoi.wkt, event_time, processing_mode="historic"
        )

        rupture_dt = convert_time(event_time).replace(tzinfo=None)
        pre_slcs = []
        reference_date = None

        if slcs:
            # Determine the maximum footprint this specific track has over the AOI
            scenes_by_date = defaultdict(list)
            for s in slcs:
                scenes_by_date[s["date"][:10]].append(s)

            max_track_area = 0
            for date_str, scenes in scenes_by_date.items():
                union_geom = unary_union([s["geometry"] for s in scenes])
                intersection_area = union_geom.intersection(aoi).area
                if intersection_area > max_track_area:
                    max_track_area = intersection_area

            # Filter for pre-seismic scenes and select the closest valid date
            valid_pre_scenes = [
                s for s in slcs if datetime.strptime(s["date"], "%Y-%m-%dT%H:%M:%SZ") < rupture_dt
            ]

            if valid_pre_scenes:
                # Sort dates newest to oldest (closest to earthquake first)
                dates = sorted(list(set(s["date"][:10] for s in valid_pre_scenes)), reverse=True)

                for ref_date in dates:
                    candidate_scenes = [s for s in valid_pre_scenes if s["date"][:10] == ref_date]
                    union_geom = unary_union([s["geometry"] for s in candidate_scenes])
                    intersection_area = union_geom.intersection(aoi).area

                    # Ensure coverage is at least 95% of the track's max expected footprint
                    if max_track_area > 0 and (intersection_area / max_track_area) > 0.95:
                        pre_slcs = [s["fileID"].removesuffix("-SLC") for s in candidate_scenes]
                        reference_date = ref_date
                        break  # Successfully found valid coverage

        if not pre_slcs:
            print(
                f"No fully overlapping pre-seismic coverage found for {track_key}. Skipping track."
            )
            continue

        # Create Partial Job List
        # We leave 'granules' (post-seismic) empty for now
        # We fill 'secondary_granules' (pre-seismic)
        job_filename = f"job_{title}_{track_key}_partial.json"

        # Create the standard HYP3 structure
        job_json = make_job_json(
            title, event_id, flight_direction, path_number, [], pre_slcs, resolution
        )

        # Save Partial File
        with open(job_filename, "w") as f:
            json.dump([job_json], f, indent=4)  # List of 1 job

        # Add to tracks info
        tracks_info[track_key] = {
            "flight_direction": flight_direction,
            "path_number": path_number,
            "frame_numbers": frame_numbers,
            "partial_job_file": job_filename,
            "reference_date": reference_date,
            "status": "AWAITING_POST_SEISMIC",
        }
        print(
            f"  Initialized track {track_key}. Pre-seismic date: {reference_date}. Waiting for post-seismic."
        )

    if tracks_info:
        tracker[event_id] = {
            "title": title,
            "time": event_time,
            "aoi": mapping(aoi),
            "tracks": tracks_info,
        }
        save_tracker(tracker)
        print(f"Added {title} to tracking file with {len(tracks_info)} tracks.")
    else:
        print(f"No valid tracks initialized for {title}.")


def check_tracker_for_updates(do_processing=False, send_email_flag=False):
    """
    Iterates through the tracking file using a unified state machine.
    - Local Machine (do_processing=True): Looks for 'AWAITING_POST_SEISMIC', runs topsApp, updates state to 'READY_FOR_EMAIL'.
    - GitHub Action (send_email_flag=True): Looks for 'READY_FOR_EMAIL', sends email, removes event.

    Args:
        do_processing (bool): If True, this machine will perform local processing for tracks awaiting post-seismic data.
        send_email_flag (bool): If True, this machine will send emails for tracks that are ready for email.
    """
    tracker = load_tracker()
    if not tracker:
        return

    events_to_remove = []
    completed_jobs_summary = []

    for event_id, event_data in tracker.items():
        title = event_data["title"]
        event_time = event_data["time"]
        tracks = event_data["tracks"]

        tracks_to_remove = []

        for track_key, track_info in tracks.items():
            # Local Processing
            if track_info["status"] == "AWAITING_POST_SEISMIC":
                if not do_processing:
                    continue

                flight_dir = track_info["flight_direction"]
                path_num = track_info["path_number"]

                # Reconstruct the AOI geometry from the tracker dictionary
                aoi_geom = shape(event_data["aoi"])

                slcs = get_SLCs(
                    flight_dir, path_num, aoi_geom.wkt, event_time, processing_mode="forward"
                )

                rupture_dt = convert_time(event_time).replace(tzinfo=None)
                post_slcs = []
                secondary_date = None

                if slcs:
                    # Determine the maximum footprint this specific track has over the AOI
                    scenes_by_date = defaultdict(list)
                    for s in slcs:
                        scenes_by_date[s["date"][:10]].append(s)

                    max_track_area = 0
                    for date_str, scenes in scenes_by_date.items():
                        union_geom = unary_union([s["geometry"] for s in scenes])
                        intersection_area = union_geom.intersection(aoi_geom).area
                        if intersection_area > max_track_area:
                            max_track_area = intersection_area

                    # Filter for post-seismic scenes and select the closest valid date
                    slcs.sort(key=lambda x: x["date"])
                    valid_post_scenes = [
                        s
                        for s in slcs
                        if datetime.strptime(s["date"], "%Y-%m-%dT%H:%M:%SZ") > rupture_dt
                    ]

                    if valid_post_scenes:
                        # Sort dates oldest to newest (closest to earthquake first)
                        post_dates = sorted(list(set(s["date"][:10] for s in valid_post_scenes)))

                        for sec_date in post_dates:
                            candidate_scenes = [
                                s for s in valid_post_scenes if s["date"][:10] == sec_date
                            ]
                            union_geom = unary_union([s["geometry"] for s in candidate_scenes])
                            intersection_area = union_geom.intersection(aoi_geom).area

                            # Ensure coverage is at least 95% of the track's max expected footprint
                            if max_track_area > 0 and (intersection_area / max_track_area) > 0.95:
                                post_slcs = [
                                    s["fileID"].removesuffix("-SLC") for s in candidate_scenes
                                ]
                                secondary_date = sec_date
                                break  # Successfully found valid coverage

                if post_slcs:
                    pre_seismic_date = track_info["reference_date"]
                    post_seismic_date = secondary_date

                    partial_file = track_info["partial_job_file"]
                    try:
                        with open(partial_file, "r") as f:
                            job_list = json.load(f)
                            job = job_list[0]

                        job["job_parameters"]["granules"] = post_slcs

                        older_date_str = pre_seismic_date.split("T")[0].replace("-", "")
                        newer_date_str = post_seismic_date.split("T")[0].replace("-", "")
                        pair_folder_name = (
                            f"{flight_dir}{int(path_num):03d}_{older_date_str}_{newer_date_str}"
                        )
                        processing_dir = os.path.join(
                            config.root_dir,
                            title,
                            f"{flight_dir}{int(path_num):03d}",
                            "coseismic",
                            pair_folder_name,
                        )

                        print(f"    Starting automatic processing for {pair_folder_name}")
                        try:
                            run_dockerized_topsApp(job, processing_dir)
                            track_info["processing_status"] = "Success"

                            # Delete raw SLCs and intermediate files if final .nc product exists
                            nc_files = glob.glob(
                                os.path.join(processing_dir, "**", "*.nc"), recursive=True
                            )

                            if nc_files:
                                print(
                                    f"    Output .nc file found. Cleaning up SLCs and heavy intermediate files in {pair_folder_name}..."
                                )
                                import shutil

                                # Delete raw .zip and .SAFE files
                                for slc_zip in glob.glob(
                                    os.path.join(processing_dir, "S1[A-D]*.zip")
                                ):
                                    os.remove(slc_zip)
                                for slc_safe in glob.glob(
                                    os.path.join(processing_dir, "S1[A-D]*.SAFE")
                                ):
                                    shutil.rmtree(slc_safe, ignore_errors=True)

                                # Delete intermediate ISCE2 folders
                                intermediate_dirs = [
                                    "geom_reference",
                                    "ion",
                                    "fine_interferogram",
                                    "fine_offsets",
                                    "fine_coreg",
                                    "mask",
                                    "reference",
                                    "secondary",
                                    "PICKLE",
                                    "aux_cal",
                                    "orbits",
                                ]
                                for idir in intermediate_dirs:
                                    dir_path = os.path.join(processing_dir, idir)
                                    if os.path.exists(dir_path):
                                        shutil.rmtree(dir_path, ignore_errors=True)
                            else:
                                print(
                                    f"    WARNING: No final .nc product found in {pair_folder_name}. Retaining raw and intermediate files for debugging."
                                )

                        except Exception as e:
                            print(f"    Processing failed for {pair_folder_name}: {e}")
                            track_info["processing_status"] = f"Failed: {str(e)}"

                        completed_filename = os.path.join(
                            processing_dir, f"job_{title}_{track_key}_COMPLETED.json"
                        )
                        with open(completed_filename, "w") as f:
                            json.dump([job], f, indent=4)

                        # Notify GitHub Actions that the job is done and ready for email
                        track_info["status"] = "READY_FOR_EMAIL"
                        track_info["dates"] = (
                            f"{pre_seismic_date} (Pre) - {post_seismic_date} (Post)"
                        )
                        track_info["location"] = processing_dir

                        if os.path.exists(partial_file):
                            os.remove(partial_file)

                    except Exception as e:
                        print(f"    Error processing partial file {partial_file}: {e}")
                else:
                    print(f"    No post-seismic data yet.")

            # Email with Github Actions
            elif track_info["status"] == "READY_FOR_EMAIL":
                if not send_email_flag:
                    continue

                completed_jobs_summary.append(
                    {
                        "title": title,
                        "track": track_key,
                        "dates": track_info.get("dates", "Unknown"),
                        "status": track_info.get("processing_status", "Unknown"),
                        "location": track_info.get("location", "Unknown"),
                    }
                )

                # Determine whether to delete or quarantine based on success/failure
                if track_info.get("processing_status", "").startswith("Failed"):
                    print(f"    Job {track_key} failed. Moving to FAILED_NEEDS_ATTENTION state.")
                    track_info["status"] = "FAILED_NEEDS_ATTENTION"
                else:
                    # Mark successful track for removal
                    tracks_to_remove.append(track_key)

        # Remove completed tracks
        for tk in tracks_to_remove:
            del tracks[tk]

        # If no tracks left, mark event for removal
        if not tracks:
            events_to_remove.append(event_id)

    # Clean up tracker
    for eid in events_to_remove:
        print(f"Event {tracker[eid]['title']} fully processed. Removing from tracker.")
        del tracker[eid]

    save_tracker(tracker)

    # Send "PROCESSING COMPLETED" email
    if completed_jobs_summary and send_email_flag:
        has_failures = any(item["status"].startswith("Failed") for item in completed_jobs_summary)
        status_tag = "WITH FAILURES" if has_failures else "SUCCESS"
        subject = (
            f"PROCESSING COMPLETED ({status_tag}): {len(completed_jobs_summary)} Jobs Processed"
        )

        # Start HTML body with inline CSS for universal email client support
        body = """
        <html>
            <body style="font-family: Arial, sans-serif; color: #333; line-height: 1.6; margin: 0; padding: 20px;">
            <div style="max-width: 600px; margin: 0 auto; border: 1px solid #e0e0e0; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05);">
                <h2 style="color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; margin-top: 0;">Automated Processing Update</h2>
                <p>The following SAR processing jobs have concluded:</p>
        """

        # Generate a styled card for each completed job
        for item in completed_jobs_summary:
            status_color = "#e74c3c" if item["status"].startswith("Failed") else "#27ae60"

            body += f"""
                <div style="background-color: #f8f9fa; padding: 15px; margin-bottom: 15px; border-left: 5px solid {status_color}; border-radius: 4px;">
                <p style="margin: 0 0 5px;"><strong>Event:</strong> {item["title"]}</p>
                <p style="margin: 0 0 5px;"><strong>Track:</strong> {item["track"]}</p>
                <p style="margin: 0 0 5px;"><strong>Dates:</strong> {item["dates"]}</p>
                <p style="margin: 0 0 5px;"><strong>Status:</strong> <span style="color: {status_color}; font-weight: bold;">{item["status"]}</span></p>
                <p style="margin: 0;"><strong>Location:</strong> <br>
                    <code style="background: #e9ecef; padding: 4px; display: block; margin-top: 5px; word-wrap: break-word; font-size: 12px; color: #c0392b;">
                    {item["location"]}
                    </code>
                </p>
                </div>
            """

        # Close the HTML body
        body += """
                <p style="font-size: 12px; color: #95a5a6; border-top: 1px solid #e0e0e0; padding-top: 15px; margin-top: 20px; text-align: center;">
                This is an automated message from the COSEIS pipeline.
                </p>
            </div>
            </body>
        </html>
        """

        print("Sending completion email to secondary recipients...")
        send_email(subject, body, recipients=config.SECONDARY_RECIPIENTS)
