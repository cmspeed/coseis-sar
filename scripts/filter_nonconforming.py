#!/usr/bin/env python3
"""
filter_polarization.py

Filters jobs from a JSON file whose granules or secondary_granules
contain something other than IW and 1SDV/1SSV, writes the non-conforming
jobs to one JSON file, and writes the conforming jobs to another JSON file.
"""

import json
import argparse
import sys

def is_non_standard_granule(granule: str) -> bool:
    """Return True if granule does NOT contain both 'IW' and '1SDV' or '1SSV'."""
    # A standard granule must contain 'IW' AND either '1SDV' or '1SSV'.
    return not ("IW" in granule and ("1SDV" in granule or "1SSV" in granule))

def main():
    parser = argparse.ArgumentParser(description="Filter jobs by granule polarization.")
    parser.add_argument("input_json", help="Path to the input jobs JSON file.")
    parser.add_argument("non_conforming_json", help="Path to the output JSON file for non-conforming jobs.")
    parser.add_argument("conforming_json", help="Path to the output JSON file for conforming (filtered) jobs.")
    args = parser.parse_args()

    # Load input JSON file
    try:
        with open(args.input_json) as f:
            jobs = json.load(f)
    except Exception as e:
        print(f"Error reading input JSON file: {e}", file=sys.stderr)
        sys.exit(1)

    # Filter jobs
    non_conforming_jobs = []
    conforming_jobs = []
    
    for job in jobs:
        granules = job.get("job_parameters", {}).get("granules", [])
        secondary = job.get("job_parameters", {}).get("secondary_granules", [])
        all_granules = granules + secondary

        # Check if ANY granule is non-standard
        if any(is_non_standard_granule(g) for g in all_granules):
            non_conforming_jobs.append(job)
        else:
            # If all granules are standard, add to the conforming list
            conforming_jobs.append(job)

    # Write non-conforming jobs to a separate JSON file
    try:
        with open(args.non_conforming_json, "w") as f:
            json.dump(non_conforming_jobs, f, indent=2)
        print(f"Wrote {len(non_conforming_jobs)} non-conforming jobs to {args.non_conforming_json}")
    except Exception as e:
        print(f"Error writing non-conforming JSON file: {e}", file=sys.stderr)
        sys.exit(1)

    # Write conforming jobs to another JSON file
    try:
        with open(args.conforming_json, "w") as f:
            json.dump(conforming_jobs, f, indent=2)
        print(f"Wrote {len(conforming_jobs)} conforming jobs to {args.conforming_json}")
    except Exception as e:
        print(f"Error writing conforming JSON file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()