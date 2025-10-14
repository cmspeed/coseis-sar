#!/usr/bin/env python3

import json
import os
import subprocess
import sys

# ===== CONFIGURATION =====
JSON_FILE = "all_jobs_with_year_remaining_3.json"
DEST_BUCKET = "aria-share"
DEST_PREFIX = "COSEIS-ONE_STOP_SHOP/HISTORIC_EVENTS"
AWS_PROFILE = "saml-pub"
LOCAL_ROOT = "S3_local_copy"  # local root folder where files already exist
# =========================

# Counters

success_count = 0
failures = []


def upload_file(local_path, dest_bucket, dest_key):
    """Upload an existing local file to destination bucket."""
    global success_count, failures

    if not os.path.exists(local_path):
        print(f"Local file missing: {local_path}")
        failures.append(f"Missing local file: {local_path}")
        return

    cmd_upload = [
        "aws", "s3", "cp",
        local_path,
        f"s3://{dest_bucket}/{dest_key}",
        "--profile", AWS_PROFILE
    ]
    try:
        subprocess.run(cmd_upload, check=True)
        print(f"Uploaded: {local_path} → {dest_bucket}/{dest_key}")
        success_count += 1
    except subprocess.CalledProcessError as e:
        print(f"Error uploading {dest_key}: {e}")
        failures.append(f"Upload failed: {dest_key}")


def main():
    global success_count, failures

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

        # ---- Upload .nc files ----
        for file_info in job.get('files', []):
            filename = file_info['filename']
            product_name = os.path.splitext(filename)[0]

            dest_key = f"{DEST_PREFIX}/{year}/{eq_name}/{product_name}/{filename}"
            local_path = os.path.join(LOCAL_ROOT, year, eq_name, product_name, filename)

            upload_file(local_path, DEST_BUCKET, dest_key)

        # ---- Upload .png browse images ----
        for img_url in job.get('browse_images', []):
            filename = os.path.basename(img_url)
            product_name = os.path.splitext(filename)[0]

            dest_key = f"{DEST_PREFIX}/{year}/{eq_name}/{product_name}/{filename}"
            local_path = os.path.join(LOCAL_ROOT, year, eq_name, product_name, filename)

            upload_file(local_path, DEST_BUCKET, dest_key)

    # ---- Summary ----
    print("\n===== UPLOAD SUMMARY =====")
    print(f"Total files successfully uploaded: {success_count}")
    if failures:
        print(f"Total failures: {len(failures)}")
        for f in failures:
            print(f" - {f}")
    else:
        print("No failures!")
    print("==========================")
    print("All done.")


if __name__ == "__main__":
    main()
