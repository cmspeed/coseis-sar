#!/usr/bin/env python3
import json
import sys
import csv
from collections import defaultdict
from pathlib import Path

def parse_name(name: str):
    """
    Parse job 'name' field into (earthquake_id, track_id).
    Example: "m_60_37_km_wsw_of_asadabad_afghanistan-D078"
    Returns:
        earthquake_id = "m_60_37_km_wsw_of_asadabad_afghanistan"
        track_id      = "D078"
    """
    if "-" not in name:
        raise ValueError(f"Unexpected name format: {name}")
    earthquake_id, track_id = name.rsplit("-", 1)
    return earthquake_id, track_id

def analyze_jobs(json_path: Path):
    with open(json_path) as f:
        jobs = json.load(f)

    total_jobs = len(jobs)
    eq_to_tracks = defaultdict(set)

    for job in jobs:
        name = job.get("name", "")
        earthquake_id, track_id = parse_name(name)
        eq_to_tracks[earthquake_id].add(track_id)

    num_unique_eqs = len(eq_to_tracks)

    stats = {
        "unique_earthquakes": num_unique_eqs,
        "jobs_total": total_jobs,
        "earthquakes": {
            eq: {
                "num_tracks": len(tracks),
                "tracks": sorted(tracks),
            }
            for eq, tracks in eq_to_tracks.items()
        }
    }
    return stats

def write_csv(stats, csv_path: Path):
    with open(csv_path, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(["Earthquake ID", "Num Tracks", "Tracks"])
        # Write rows
        for eq, details in stats["earthquakes"].items():
            writer.writerow([eq, details["num_tracks"], " ".join(details["tracks"])])
        # Write summary rows at the bottom
        writer.writerow([])
        writer.writerow(["Total unique earthquakes", stats["unique_earthquakes"]])
        writer.writerow(["Total jobs", stats["jobs_total"]])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <jobs.json> <output.csv>")
        sys.exit(1)

    json_path = Path(sys.argv[1])
    csv_path = Path(sys.argv[2])

    if not json_path.exists():
        print(f"Error: {json_path} does not exist")
        sys.exit(1)

    stats = analyze_jobs(json_path)
    write_csv(stats, csv_path)

    print(f"Wrote summary statistics to {csv_path}")