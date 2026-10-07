#!/bin/bash
# Forward-mode processing run for the processing machine's cron job, e.g.:
#   5,15 * * * * /path/to/repo/ops/run_coseis_forward.sh
# Checks tracked events for post-event Sentinel-1 data, runs topsApp, and commits the
# updated tracker (active_jobs/, including active_jobs/partials/) to main.
#
# Optional environment overrides (see aria_coseis.config): COSEIS_DATA_DIR for topsApp
# products (default <repo>/data), COSEIS_LOCK_FILE (default /tmp/coseis_processing.lock).

# Setup Environment: ~/.bashrc normally initializes conda/mamba ("conda init")
source ~/.bashrc

# Fallback if ~/.bashrc didn't define the mamba shell function: load it from the conda install
if ! type mamba &> /dev/null; then
    CONDA_BASE="$(conda info --base 2> /dev/null)"
    if [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/etc/profile.d/mamba.sh" ]; then
        source "$CONDA_BASE/etc/profile.d/conda.sh"
        source "$CONDA_BASE/etc/profile.d/mamba.sh"
    fi
fi

mamba activate coseis-sar

# Run from the repository root (this script lives in <repo>/ops)
REPO_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." &> /dev/null && pwd )"
cd "$REPO_DIR"
mkdir -p logs
LOG_FILE="logs/forward.log"

# Check the lock file BEFORE doing anything with Git (same default as aria_coseis.config.LOCK_FILE)
LOCK_FILE="${COSEIS_LOCK_FILE:-/tmp/coseis_processing.lock}"
if [ -f "$LOCK_FILE" ]; then
    echo "$(date): Previous processing run still active. Bash script exiting." >> "$LOG_FILE"
    exit 0
fi

# Sync with Github
git checkout main

# Pull the latest tracker state that the GitHub Action just updated
git pull --rebase origin main

# Execute local processing (discovery disabled via --process_only)
python -m aria_coseis --forward --pairing coseismic --resolution 30 --do_processing --process_only >> "$LOG_FILE" 2>&1

# Sync with Github: stage tracker edits, new files and deletions (finished partials) under active_jobs/
git add -A active_jobs/ || true

if ! git diff --cached --quiet; then
    git commit -m "Local processing: update COSEIS tracking state and remove finished partials [skip ci]"
    git pull --rebase origin main
    git push origin main
else
    echo "No processing completed this run; tracking state unchanged." >> "$LOG_FILE"
fi
