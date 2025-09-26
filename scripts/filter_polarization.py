#!/usr/bin/env python3
"""
filter_polarization.py

Filters jobs from a JSON file whose granules or secondary_granules
contain something other than IW and 1SDV, and writes them to a new JSON file.
"""

import json
import argparse
import sys

def is_non_standard_granule(granule: str) -> bool:
    """Return True if granule does NOT contain both 'IW' and '1SDV'."""
    return not ("IW" in granule and ("1SDV" in granule or "1SSV" in granule))

def main():
    parser = argparse.ArgumentParser(description="Filter jobs by granule polarization.")
    parser.add_argument("input_json", help="Path to the input jobs JSON file")
    parser.add_argument("output_json", help="Path to the output JSON file for non-conforming jobs")
    args = parser.parse_args()

    # Load input JSON file
    try:
        with open(args.input_json) as f:
            jobs = json.load(f)
    except Exception as e:
        print(f"Error reading input JSON file: {e}", file=sys.stderr)
        sys.exit(1)

    # Filter non-conforming jobs
    non_standard_jobs = []
    for job in jobs:
        granules = job.get("job_parameters", {}).get("granules", [])
        secondary = job.get("job_parameters", {}).get("secondary_granules", [])
        all_granules = granules + secondary

        if any(is_non_standard_granule(g) for g in all_granules):
            non_standard_jobs.append(job)

    # Write non-conforming jobs to output JSON file
    try:
        with open(args.output_json, "w") as f:
            json.dump(non_standard_jobs, f, indent=2)
        print(f"Wrote {len(non_standard_jobs)} non-conforming jobs to {args.output_json}")
    except Exception as e:
        print(f"Error writing output JSON file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
