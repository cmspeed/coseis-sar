#!/usr/bin/env python3

import json
import os
import subprocess
import sys

# ===== CONFIGURATION =====
JSON_FILE = "all_jobs_with_year_remaining.json"
DEST_BUCKET = "aria-share"
DEST_PREFIX = "COSEIS-ONE_STOP_SHOP/HISTORIC_EVENTS"
AWS_PROFILE = "saml-pub"
LOCAL_ROOT = "S3_local_copy"  # local root folder to retain downloaded files
# =========================

# Counters
success_count = 0
failures = []

def copy_file_two_step(src_url, dest_bucket, dest_key, local_path):
    """Download a file from source URL (anonymous) and upload to destination bucket, keeping a local copy."""
    global success_count, failures

    # Ensure local folder exists
    os.makedirs(os.path.dirname(local_path), exist_ok=True)

    # Step 1: download anonymously using AWS CLI (or curl if URL)
    if src_url.startswith("http"):
        cmd_download = ["curl", "-sSL", src_url, "-o", local_path]
    else:
        # S3 path, use --no-sign-request
        bucket, key = src_url.split("/", 1)
        cmd_download = [
            "aws", "s3", "cp",
            f"s3://{bucket}/{key}",
            local_path,
            "--no-sign-request"
        ]

    try:
        subprocess.run(cmd_download, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error downloading {src_url}: {e}")
        failures.append(f"Download failed: {src_url}")
        return

    # Step 2: upload to destination using profile
    cmd_upload = [
        "aws", "s3", "cp",
        local_path,
        f"s3://{dest_bucket}/{dest_key}",
        "--profile", AWS_PROFILE
    ]
    try:
        subprocess.run(cmd_upload, check=True)
        print(f"Copied: {src_url} → {dest_bucket}/{dest_key}")
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

        # ---- Copy .nc files ----
        for file_info in job.get('files', []):
            filename = file_info['filename']
            src_bucket = file_info['s3']['bucket']
            src_key = file_info['s3']['key']

            product_name = os.path.splitext(filename)[0]
            dest_key = f"{DEST_PREFIX}/{year}/{eq_name}/{product_name}/{filename}"
            local_path = os.path.join(LOCAL_ROOT, year, eq_name, product_name, filename)

            src_s3_url = f"{src_bucket}/{src_key}"
            copy_file_two_step(src_s3_url, DEST_BUCKET, dest_key, local_path)

        # ---- Copy .png browse images ----
        for img_url in job.get('browse_images', []):
            filename = os.path.basename(img_url)
            product_name = os.path.splitext(filename)[0]
            dest_key = f"{DEST_PREFIX}/{year}/{eq_name}/{product_name}/{filename}"
            local_path = os.path.join(LOCAL_ROOT, year, eq_name, product_name, filename)

            copy_file_two_step(img_url, DEST_BUCKET, dest_key, local_path)

    # ---- Summary ----
    print("\n===== COPY SUMMARY =====")
    print(f"Total files successfully copied: {success_count}")
    if failures:
        print(f"Total failures: {len(failures)}")
        for f in failures:
            print(f" - {f}")
    else:
        print("No failures!")
    print("=========================")
    print("All done.")

if __name__ == "__main__":
    main()