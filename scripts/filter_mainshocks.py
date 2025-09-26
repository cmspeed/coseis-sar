import json
import re
import sys
from pathlib import Path

def extract_magnitude(name: str) -> int:
    """Extract the first integer from the job name (magnitude)."""
    match = re.match(r"m_(\d+)", name)
    return int(match.group(1)) if match else -1

def filter_jobs(input_file: str, output_file: str, removed_file: str):
    # Load jobs
    with open(input_file, "r") as f:
        jobs = json.load(f)

    unique = {}
    removed = []

    for job in jobs:
        granules = tuple(job["job_parameters"]["granules"])
        secondary = tuple(job["job_parameters"]["secondary_granules"])
        key = (granules, secondary)

        magnitude = extract_magnitude(job["name"])

        if key not in unique:
            unique[key] = job
        else:
            # Compare with the one we already kept
            existing = unique[key]
            if magnitude > extract_magnitude(existing["name"]):
                removed.append(existing)
                unique[key] = job
            else:
                removed.append(job)

    # Save filtered jobs
    with open(output_file, "w") as f:
        json.dump(list(unique.values()), f, indent=4)

    # Save removed jobs
    with open(removed_file, "w") as f:
        json.dump(removed, f, indent=4)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: python {sys.argv[0]} input.json kept.json removed.json")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    removed_file = sys.argv[3]

    filter_jobs(input_file, output_file, removed_file)