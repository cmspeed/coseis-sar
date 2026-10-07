"""Local dockerized topsApp execution."""

from __future__ import annotations

import os
import subprocess
from typing import Any


def create_directories_from_json(
    eq_jsons: list[list[dict[str, Any]]], root_dir: str
) -> tuple[list[list[str]], int]:
    """
    Create directories for each group of SLCs based on the JSON data provided. These directories will be used to store the outputs of the dockerized topsApp.
    The directories are created in the root directory specified and will have a name following this format: 'flight_directionpath_number_secondary_date_reference_date'.
    For example 'A012_20220101_20220112' for an ascending path 12 earthquake with secondary date 2022-01-01 and reference date 2022-01-12.
    :param eq_jsons: List of JSON objects containing the parameters for each pair of SLCs
    :param root_dir: Root directory where the directories will be created
    :return: List of directory names
    """
    dirnames = []
    total = 0
    for isce_jsons in eq_jsons:
        sub_dirnames = []  # To hold the directories for each group in `isce_jsons`
        for json_data in isce_jsons:
            title = json_data["title"]
            timing = json_data["timing"]
            flight_direction = (
                "A" if json_data["flight-direction"] == "ASCENDING" else "D"
            )  # Reformat 'fight-direction' to shorten dirname
            path_number = json_data["path-number"]
            secondary_date = json_data["secondary-date"].replace("-", "")
            reference_date = json_data["reference-date"].replace("-", "")

            # Build the full path
            base_path = os.path.join(root_dir, title, flight_direction + path_number, timing)
            sub_path = f"{flight_direction}{path_number}_{secondary_date}_{reference_date}"
            full_path = os.path.join(base_path, sub_path)

            # Create directories, ensuring no overwriting
            os.makedirs(full_path, exist_ok=True)
            sub_dirnames.append(full_path)
            print(f"Created: {full_path}")
            total += 1
        dirnames.append(sub_dirnames)
    return dirnames, total


def run_dockerized_topsApp(json_data: dict[str, Any], working_dir: str) -> None:
    """
    Run dockerized topsApp InSAR processing workflow using the provided JSON data.
    Outputs are added to the root dir + an extension for each pair.
    :param json_data: JSON object containing the parameters for dockerized topsApp
    :param working_dir: Working directory where the outputs will be stored
    """
    # Extract the parameters from the JSON data
    params = json_data["job_parameters"]
    reference_scenes = " ".join(params["granules"])
    secondary_scenes = " ".join(params["secondary_granules"])

    # Define static/dynamic parameters
    frame_id = params.get("frame_id", -1)
    estimate_ionosphere_delay = str(params.get("estimate_ionosphere_delay"))
    esd_coherence_threshold = str(params.get("esd_coherence_threshold"))
    compute_solid_earth_tide = str(params.get("compute_solid_earth_tide"))
    goldstein_filter_power = str(params.get("goldstein_filter_power"))
    output_resolution = str(params.get("output_resolution"))
    unfiltered_coherence = str(params.get("unfiltered_coherence", True))
    dense_offsets = str(params.get("dense_offsets", True))

    # Construct the command
    cmd = [
        "conda",
        "run",
        "-n",
        "topsapp_env_trappist_python11",
        "taskset",
        "-c",
        "0-16",
        "isce2_topsapp",
        "--reference-scenes",
        reference_scenes,
        "--secondary-scenes",
        secondary_scenes,
        "--frame-id",
        str(frame_id),
        "--estimate-ionosphere-delay",
        estimate_ionosphere_delay,
        "--esd-coherence-threshold",
        esd_coherence_threshold,
        "--compute-solid-earth-tide",
        compute_solid_earth_tide,
        "--goldstein-filter-power",
        goldstein_filter_power,
        "--output-resolution",
        output_resolution,
        "--unfiltered-coherence",
        unfiltered_coherence,
        "--dense-offsets",
        dense_offsets,
        "++process",
        "coseis_sar",
    ]

    # Set Environment Variables
    env = os.environ.copy()
    env["XLA_PYTHON_CLIENT_MEM_FRACTION"] = ".20"

    print("=========================================")
    print(f"Running dockerized topsApp in {working_dir}...")
    print(f"Command: {' '.join(cmd)}")
    print("=========================================")

    # Run the command
    try:
        # Check if working dir exists, create if not
        os.makedirs(working_dir, exist_ok=True)

        # Execute
        result = subprocess.run(
            cmd,
            cwd=working_dir,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
        print("Processing completed successfully.")

        # Optional: Save a log file in that directory
        with open(os.path.join(working_dir, "topsapp_stdout.log"), "w") as log:
            log.write(result.stdout)

    except subprocess.CalledProcessError as e:
        print("Error occurred while running topsApp.")
        print("Stderr:\n", e.stderr)
        # Optional: Write error log
        with open(os.path.join(working_dir, "topsapp_error.log"), "w") as log:
            log.write(e.stderr)
        raise e  # Re-raise to handle it in the calling function if needed
