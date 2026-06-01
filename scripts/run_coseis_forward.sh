#!/bin/bash

# Setup Environment
source ~/.bashrc

# Dynamically find the Conda/Miniforge installation in the user's home directory
if [ -f "$HOME/tools/miniforge3/etc/profile.d/mamba.sh" ]; then
    source "$HOME/tools/miniforge3/etc/profile.d/mamba.sh"
else
    source "$HOME/tools/miniforge3/etc/profile.d/conda.sh"
fi

mamba activate coseis-sar

# Dynamically set the working directory to wherever this script lives
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "$SCRIPT_DIR"

# Check the lock file BEFORE doing anything with Git
LOCK_FILE="/tmp/coseis_processing.lock"
if [ -f "$LOCK_FILE" ]; then
    echo "$(date): Previous processing run still active. Bash script exiting." >> log_tracking.txt
    exit 0
fi

# Sync with Github
git checkout main

# Pull the latest active_job_tracking.json that the GitHub Action just updated
git pull --rebase origin main 

# Execute local processing (Scouting disabled via --process_only)
python coseis.py --forward --pairing coseismic --do_processing --process_only >> log_tracking.txt 2>&1

# Sync with Github 
# Add the updated tracker JSON
git add active_job_tracking.json

# Stage any partial files that the python script deleted locally
git add -u "*_partial.json" || true
if ! git diff --cached --quiet; then
    git commit -m "Local processing: update COSEIS tracking state and remove finished partials [skip ci]"
    git pull --rebase origin main
    git push origin main
else
    echo "No processing completed this run; tracking state unchanged." >> log_tracking.txt
fi