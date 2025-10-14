import json
from pathlib import Path

def update_output_resolution(json_path, out_path, new_resolution=30):
    # Load JSON
    with open(json_path, "r") as f:
        jobs = json.load(f)

    # Update field
    for job in jobs:
        if "job_parameters" in job:
            job["job_parameters"]["output_resolution"] = new_resolution

    # Write updated JSON
    out_path = Path(out_path)
    with open(out_path, "w") as f:
        json.dump(jobs, f, indent=4)

    print(f"Updated output_resolution to {new_resolution} for {len(jobs)} jobs → {out_path}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Update output_resolution in job_list.json")
    parser.add_argument("json_file", help="Path to input job_list.json")
    parser.add_argument("out_file", help="Path to save updated JSON")
    parser.add_argument(
        "--resolution",
        type=int,
        default=30,
        help="New output_resolution value (default: 30)"
    )
    args = parser.parse_args()

    update_output_resolution(args.json_file, args.out_file, args.resolution)