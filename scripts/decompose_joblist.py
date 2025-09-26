import json
from pathlib import Path
from collections import defaultdict

def decompose_job_list(json_path, out_dir):
    # Load JSON
    with open(json_path, "r") as f:
        jobs = json.load(f)

    grouped = defaultdict(list)
    mismatched = []

    for job in jobs:
        granules_len = len(job["job_parameters"].get("granules", []))
        secondary_len = len(job["job_parameters"].get("secondary_granules", []))

        if granules_len == secondary_len:
            grouped[granules_len].append(job)
        else:
            mismatched.append(job)

    # Ensure output dir exists
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Write matched groups
    for nframes, group in grouped.items():
        out_file = out_dir / f"job_list_{nframes}frame.json"
        with open(out_file, "w") as f:
            json.dump(group, f, indent=4)
        print(f"Wrote {len(group)} jobs → {out_file}")

    # Write mismatched group if any
    if mismatched:
        out_file = out_dir / "job_list_mismatched.json"
        with open(out_file, "w") as f:
            json.dump(mismatched, f, indent=4)
        print(f"Wrote {len(mismatched)} jobs → {out_file}")

    # Print summary
    print("\nSummary:")
    for nframes, group in grouped.items():
        print(f"  {nframes}-frame: {len(group)} jobs")
    if mismatched:
        print(f"  mismatched: {len(mismatched)} jobs")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description="Decompose job_list.json into subjsons by frame count"
    )
    parser.add_argument("json_file", help="Path to job_list.json")
    parser.add_argument("out_dir", help="Directory to store subjsons")
    args = parser.parse_args()

    decompose_job_list(args.json_file, args.out_dir)
