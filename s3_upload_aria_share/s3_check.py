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
LOCAL_ROOT = "S3_local_copy"  # still used for reconstructing paths
# =========================

# Counters
present_count = 0
missing = []


def check_s3_file(dest_bucket, dest_key):
    """Check if file exists in S3 bucket."""
    global present_count, missing

    cmd_check = [
        "aws", "s3", "ls",
        f"s3://{dest_bucket}/{dest_key}",
        "--profile", AWS_PROFILE
    ]
    result = subprocess.run(cmd_check, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode == 0 and result.stdout.strip():
        print(f"✓ Exists: {dest_bucket}/{dest_key}")
        present_count += 1
    else:
        print(f"✗ MISSING: {dest_bucket}/{dest_key}")
        missing.append(dest_key)


def main():
    global present_count, missing

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

    for job in jobs:
        year = job.get('year')
        if not year:
            print(f"Skipping job {job.get('name', 'UNKNOWN')} because year is missing")
            continue

        eq_name = job.get('name', 'UNKNOWN').rsplit('-', 1)[0]

        # ---- Check .nc files ----
        for file_info in job.get('files', []):
            filename = file_info['filename']
            product_name = os.path.splitext(filename)[0]

            dest_key = f"{DEST_PREFIX}/{year}/{eq_name}/{product_name}/{filename}"
            check_s3_file(DEST_BUCKET, dest_key)

        # ---- Check .png browse images ----
        for img_url in job.get('browse_images', []):
            filename = os.path.basename(img_url)
            product_name = os.path.splitext(filename)[0]

            dest_key = f"{DEST_PREFIX}/{year}/{eq_name}/{product_name}/{filename}"
            check_s3_file(DEST_BUCKET, dest_key)

    # ---- Summary ----
    print("\n===== S3 CHECK SUMMARY =====")
    print(f"Total files present: {present_count}")
    if missing:
        print(f"Total missing: {len(missing)}")
        for m in missing:
            print(f" - {m}")
    else:
        print("All files present!")
    print("============================")
    print("Done.")


if __name__ == "__main__":
    main()