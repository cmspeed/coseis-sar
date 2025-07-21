import json
import argparse
import re
from pathlib import Path


def extract_magnitude(name):
    """
    Extracts magnitude (XX) from name formatted as 'm_XX_place_earthquake'.
    Returns an integer magnitude.
    """
    match = re.match(r"m_(\d+)_", name)
    return int(match.group(1)) if match else -1


def normalize_granules(granules):
    """
    Returns a sorted tuple of granules for consistent deduplication.
    """
    return tuple(sorted(granules))


def deduplicate_jobs(jobs):
    """
    Deduplicates jobs based on granules + secondary_granules, keeping the one with the highest magnitude.
    """
    deduped = {}

    for job in jobs:
        params = job.get("job_parameters", {})
        granules = normalize_granules(params.get("granules", []))
        secondary_granules = normalize_granules(params.get("secondary_granules", []))
        key = (granules, secondary_granules)

        current_mag = extract_magnitude(job.get("name", ""))

        # Keep highest magnitude
        if key not in deduped or current_mag > extract_magnitude(deduped[key]["name"]):
            deduped[key] = job

    return list(deduped.values())


def main(input_file):
    input_path = Path(input_file)

    with input_path.open("r") as f:
        data = json.load(f)

    if not isinstance(data, list):
        print("Error: Top-level JSON structure must be a list of jobs.")
        return

    deduped_jobs = deduplicate_jobs(data)

    output_path = input_path.with_name(f"{input_path.stem}_deduplicated.json")
    with output_path.open("w") as f:
        json.dump(deduped_jobs, f, indent=2)

    print(f"Deduplicated jobs written to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Deduplicate ARIA_S1_COSEIS jobs by granule content and magnitude.")
    parser.add_argument("input_file", type=str, help="Path to job_list.json")

    args = parser.parse_args()
    main(args.input_file)