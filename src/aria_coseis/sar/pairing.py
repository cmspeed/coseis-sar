"""Sentinel-1 reference/secondary pairing and topsApp/HyP3 job JSON."""

from __future__ import annotations

from datetime import datetime
from itertools import combinations
from typing import Any

from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union

from aria_coseis.utils import convert_time


def generate_pairs(pairs: list[Any], mode: str) -> list[tuple[Any, Any]]:
    """
    Generate pairs of SLCs based on the selected pairing mode.
    :param pairs: List of SLC pairs sorted by date
    :param: mode: 'sequential' for temporally consecutive pairs, 'all' for all possible pairs, 'conseismic' for pairs bounding the rupture date only
    :return: List of SLC pairs based on the mode
    """
    if mode == "sequential":
        return [(pairs[i], pairs[i + 1]) for i in range(len(pairs) - 1)]
    elif mode == "all":
        all_pairs = list(combinations(pairs, 2))
        return all_pairs
    elif mode == "coseismic":
        return []


def find_reference_and_secondary_pairs(
    SLCs: list[dict[str, Any]],
    time: int,
    flight_direction: str,
    path_number: int,
    title: str,
    aoi: BaseGeometry,
    event_id: str,
    pairing_mode: str = "sequential",
    job_list: bool = False,
    resolution: int = 90,
) -> list[dict[str, Any]]:
    """
    Find the reference and secondary pairs of SLCs necessary to run dockerized topsApp,
    and determine whether each pair is pre-seismic, co-seismic, or post-seismic based on the rupture date and SLC dates.
    :param SLCs: List of dictionaries containing SLC fileIDs and their respective dates
    :param time: Unix timestamp representing the earthquake's origin time
    :param flight_direction: 'ASCENDING' or 'DESCENDING'
    :param path_number: Sentinel-1 path number
    :param title: USGS title of the earthquake event, used for file organization
    :param aoi: Shapely Polygon object representing the Area of Interest
    :param event_id: USGS event ID of the earthquake, used for file organization
    :param pairing_mode: 'sequential' for temporally consecutive pairs, 'all' for all possible pairs, 'coseismic' for pairs bounding the rupture date only
    :param job_list: True if the JSON objects are for HYP3 job submission, False otherwise
    :param resolution: Output resolution for the topsApp processing, default is 90m
    :return: List of JSON objects containing the parameters for each pair of SLCs
    """
    # Get the rupture date in the format YYYY-MM-DD
    rupture_date = convert_time(time)
    rupture_date = rupture_date.strftime("%Y-%m-%d")
    rupture_date_dt = convert_time(time).replace(tzinfo=None)

    # Reformatting for dictionary keys for later use
    flight_direction = "A" if flight_direction == "ASCENDING" else "D"
    path_number = f"{int(path_number):03}"

    # Pair SLCs by date
    slc_by_date = {}
    for slc in SLCs:
        date_obj = datetime.strptime(slc["date"][:10], "%Y-%m-%d")
        key = date_obj
        if key not in slc_by_date:
            slc_by_date[key] = []
        slc_by_date[key].append(slc)

    sorted_dates = sorted(slc_by_date.keys())
    initial_pairs = []

    # Determine the maximum footprint this specific track has over the AOI
    max_track_area = 0
    for date in sorted_dates:
        frames = slc_by_date[date]
        union_geom = unary_union([slc["geometry"] for slc in frames])
        intersection_area = union_geom.intersection(aoi).area
        if intersection_area > max_track_area:
            max_track_area = intersection_area

    # Filter dates based on coverage
    for date in sorted_dates:
        frames = slc_by_date[date]
        union_geom = unary_union([slc["geometry"] for slc in frames])
        intersection_area = union_geom.intersection(aoi).area

        # Only accept dates where the available data covers >95% of the track's max expected footprint
        if max_track_area > 0 and (intersection_area / max_track_area) > 0.95:
            initial_pairs.append((date, frames))

    pre_seismic = []
    post_seismic = []

    for pair in initial_pairs:
        # pair[1] is the list of frames. Take the exact datetime of the first frame.
        exact_time = datetime.strptime(pair[1][0]["date"], "%Y-%m-%dT%H:%M:%SZ")

        if exact_time < rupture_date_dt:
            pre_seismic.append(pair)
        elif exact_time > rupture_date_dt:
            post_seismic.append(pair)

    co_seismic = []

    for i in range(len(initial_pairs) - 1):
        # Extract the exact datetime for the first slice in each track pass
        exact_time_1 = datetime.strptime(initial_pairs[i][1][0]["date"], "%Y-%m-%dT%H:%M:%SZ")
        exact_time_2 = datetime.strptime(initial_pairs[i + 1][1][0]["date"], "%Y-%m-%dT%H:%M:%SZ")

        # Safely evaluate using precise hour/minute/second boundaries
        if exact_time_1 < rupture_date_dt < exact_time_2:
            co_seismic = [(initial_pairs[i], initial_pairs[i + 1])]
            break

    # Pair based on the selected pairing_mode
    if pairing_mode == "coseismic":
        paired_results = {"co-seismic": co_seismic}
    else:
        paired_results = {
            "pre_seismic": generate_pairs(pre_seismic, pairing_mode),
            "co_seismic": co_seismic,
            "post_seismic": generate_pairs(post_seismic, pairing_mode),
        }

    # Create JSON objects for each pair
    isce_jsons = []
    for timing, pairs in paired_results.items():
        for secondary, reference in pairs:
            reference_date, reference_scenes = reference
            secondary_date, secondary_scenes = secondary
            reference_scenes_ids = [slc["fileID"].removesuffix("-SLC") for slc in reference_scenes]
            secondary_scenes_ids = [slc["fileID"].removesuffix("-SLC") for slc in secondary_scenes]

            # Dynamically pull frame numbers for the JSON
            current_frame_numbers = list(set(slc["frameNumber"] for slc in reference_scenes))

            if job_list:
                json_output = make_job_json(
                    title,
                    event_id,
                    flight_direction,
                    path_number,
                    reference_scenes_ids,
                    secondary_scenes_ids,
                    resolution,
                )
            else:
                json_output = make_json(
                    title,
                    timing,
                    flight_direction,
                    path_number,
                    current_frame_numbers,
                    {"date": reference_date.strftime("%Y-%m-%d")},
                    {"date": secondary_date.strftime("%Y-%m-%d")},
                    reference_scenes_ids,
                    secondary_scenes_ids,
                )
            isce_jsons.append(json_output)
    return isce_jsons


def make_json(
    title: str,
    timing: str,
    flight_direction: str,
    path_number: int,
    frame_numbers: list[int],
    reference: str,
    secondary: str,
    reference_scenes: list[str],
    secondary_scenes: list[str],
) -> dict[str, Any]:
    """Create a JSON object containing parameters for dockerized topsApp.
    Note: Not all params here are used in the final dockerized topsApp. Some are used for file organzation.
    Note: Several params are 'hardcoded', as these should not vary between individual products.
    :param title: USGS title of the earthquake event
    :param timing: 'pre-seismic', 'co-seismic', or 'post-seismic'
    :param flight_direction: 'A' or 'D' for ascending or descending
    :param path_number: Sentinel-1 path number
    :param frame_numbers: List of intersecting Sentinel-1 frame numbers
    :param reference: Dictionary containing the reference date
    :param secondary: Dictionary containing the secondary date
    :param reference_scenes: List of reference SLC fileIDs
    :param secondary_scenes: List of secondary SLC fileIDs
    :return: JSON object containing the parameters for dockerized topsApp
    """
    # Reformatting 'fight-direction' for readability in the json
    flight_direction = "ASCENDING" if flight_direction == "A" else "DESCENDING"

    isce_json = {
        "title": title,
        "timing": timing,
        "flight-direction": flight_direction,
        "path-number": path_number,
        "frame-numbers": frame_numbers,
        "reference-date": reference["date"],
        "secondary-date": secondary["date"],
        "reference-scenes": reference_scenes,
        "secondary-scenes": secondary_scenes,
        "frame-id": -1,
        "estimate-ionosphere-delay": True,
        "esd-coherence-threshold": -1,
        "compute-solid-earth-tide": True,
        "goldstein-filter-power": 0.5,
        "output-resolution": 30,
        "unfiltered-coherence": True,
        "dense-offsets": True,
    }
    return isce_json


def make_job_json(
    title: str,
    event_id: str,
    flight_direction: str,
    path_number: int,
    reference_scenes: list[str],
    secondary_scenes: list[str],
    resolution: int,
) -> dict[str, Any]:
    """
    Create a JSON object containing parameters for dockerized topsApp on HYP3.
    :param title: USGS title of the earthquake event
    :param event_id: USGS event ID of the earthquake event
    :param flight_direction: 'A' or 'D' for ascending or descending
    :param path_number: Sentinel-1 path number
    :param reference_scenes: List of reference SLC fileIDs
    :param secondary_scenes: List of secondary SLC fileIDs
    :param resolution: Output resolution for the topsApp processing, default is 90m
    :return: JSON object containing the parameters for dockerized topsApp
    """

    job_json = {
        "name": f"{title}-{flight_direction}{path_number}",
        "job_type": "ARIA_S1_COSEIS",
        "event_id": event_id,
        "job_parameters": {
            "granules": reference_scenes,
            "secondary_granules": secondary_scenes,
            "frame_id": -1,
            "estimate_ionosphere_delay": True,
            "esd_coherence_threshold": -1,
            "compute_solid_earth_tide": True,
            "goldstein_filter_power": 0.5,
            "unfiltered_coherence": True,
            "dense_offsets": True,
            "output_resolution": resolution,
        },
    }
    return job_json
