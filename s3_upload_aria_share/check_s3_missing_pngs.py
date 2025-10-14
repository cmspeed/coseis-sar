#!/usr/bin/env python3

import json
import os
import subprocess
import sys

# ===== CONFIGURATION =====
JSON_FILE = "all_jobs_with_year.json"
DEST_BUCKET = "aria-share"
DEST_PREFIX = "COSEIS-ONE_STOP_SHOP/HISTORIC_EVENTS"
AWS_PROFILE = "saml-pub"
LOCAL_ROOT = "S3_local_copy" # still used for reconstructing paths
# =========================

# Counters
present_count = 0
missing_png_dirs = set() # Use a set to store unique missing directories

def check_s3_dir_for_png(dest_bucket, dest_dir_prefix):
    """Check if directory in S3 bucket contains a .png file and prints the full directory if missing."""
    global present_count, missing_png_dirs

    # The S3 directory path to check for contents
    s3_path = f"s3://{dest_bucket}/{dest_dir_prefix}/"

    cmd_check = [
        "aws", "s3", "ls",
        s3_path, # list contents of the directory/prefix
        "--recursive", # Check subdirectories as well, in case
        "--human-readable",
        "--profile", AWS_PROFILE
    ]

    result = subprocess.run(cmd_check, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Check for errors in the AWS command execution
    if result.returncode != 0:
        # If the directory doesn't exist or there's an error, report it as missing a PNG
        sys.stderr.write(f"Error or non-existent S3 path {s3_path}: {result.stderr.strip()}\n")
        print(f"✗ MISSING PNG in: {dest_dir_prefix}") # Print the missing directory
        missing_png_dirs.add(dest_dir_prefix)
        return

    # Check if any line in stdout (the list of files) contains a '.png' file extension
    png_found = any(line.strip().endswith('.png') for line in result.stdout.splitlines())

    if png_found:
        # ONLY update the count, DO NOT print confirmation
        present_count += 1
    else:
        # Print the full directory for missing PNGs
        print(f"✗ MISSING PNG in: {dest_dir_prefix}")
        missing_png_dirs.add(dest_dir_prefix)


def main():
    global present_count, missing_png_dirs

    # Load JSON
    try:
        with open(JSON_FILE) as f:
            jobs = json.load(f)
    except FileNotFoundError:
        print(f"Error: JSON file not found: {JSON_FILE}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON: {e}")
        sys.exit(1)

    # Use a set to store unique product directories to check only once
    product_directories_to_check = set()

    for job in jobs:
        year = job.get('year')
        if not year:
            print(f"Skipping job {job.get('name', 'UNKNOWN')} because year is missing")
            continue

        eq_name = job.get('name', 'UNKNOWN').rsplit('-', 1)[0]

        # Extract product directories from '.nc files' section
        for file_info in job.get('files', []):
            filename = file_info['filename']
            product_name = os.path.splitext(filename)[0]
            # Construct the product directory prefix
            dir_prefix = f"{DEST_PREFIX}/{year}/{eq_name}/{product_name}"
            product_directories_to_check.add(dir_prefix)

        # Extract product directories from 'browse_images' section (for redundancy/completeness)
        for img_url in job.get('browse_images', []):
            filename = os.path.basename(img_url)
            product_name = os.path.splitext(filename)[0]
            # Construct the product directory prefix
            dir_prefix = f"{DEST_PREFIX}/{year}/{eq_name}/{product_name}"
            product_directories_to_check.add(dir_prefix)

    print(f"Found {len(product_directories_to_check)} unique product directories to check. Printing only those missing a .png file.\n")
    print("Starting S3 PNG check...")

    # Iterate and check each unique product directory
    for dir_prefix in sorted(list(product_directories_to_check)):
        check_s3_dir_for_png(DEST_BUCKET, dir_prefix)

    # ---- Summary ----
    print("\n===== S3 PNG CHECK SUMMARY =====")
    print(f"Total directories checked: {len(product_directories_to_check)}")
    print(f"Total directories with PNGs found: {present_count}")
    if missing_png_dirs:
        print(f"Total directories MISSING PNGs: {len(missing_png_dirs)}")
        print("\nFull list of missing directories:")
        for m in sorted(list(missing_png_dirs)):
            print(f" - {m}")
    else:
        print("All product directories have at least one .png file! 🎉")
    print("================================")
    print("Done.")


if __name__ == "__main__":
    main()