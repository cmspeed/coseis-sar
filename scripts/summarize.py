import json
import csv
import re
from pathlib import Path

def summarize_jobs(json_path, csv_path):
    # Load JSON
    with open(json_path, "r") as f:
        jobs = json.load(f)

    rows = []
    for job in jobs:
        job_name = job["name"]

        # Split into name and suffix (-Axxx / -Dxxx)
        base, suffix = job_name.rsplit("-", 1)

        # Direction + Path
        direction = "Ascending" if suffix.startswith("A") else "Descending"
        path = suffix[1:]  # drop A/D

        rows.append({
            "Name": base,
            "Direction": direction,
            "Path": path
        })

    # Count metadata
    unique_earthquakes = len(set(r["Name"] for r in rows))
    num_jobs = len(rows)

    # Write CSV
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)

        # Header rows
        writer.writerow([f"Number of unique earthquakes: {unique_earthquakes}"])
        writer.writerow([f"Number of unique jobs: {num_jobs}"])
        
        # Column header
        writer.writerow(["Name", "Direction", "Path"])
        
        # Data rows
        for r in rows:
            writer.writerow([r["Name"], r["Direction"], r["Path"]])

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Summarize job_list.json into a CSV")
    parser.add_argument("json_file", help="Path to job_list.json")
    parser.add_argument("csv_file", help="Path to output CSV")
    args = parser.parse_args()

    summarize_jobs(args.json_file, args.csv_file)
