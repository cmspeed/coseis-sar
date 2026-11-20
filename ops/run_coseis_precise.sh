#!/bin/bash
# Precise-orbit reruns for the processing machine's cron job, e.g.:
#   30 */6 * * * COSEIS_DATA_DIR=... COSEIS_ARCHIVE_DIR=... /path/to/repo/ops/run_coseis_precise.sh
# Reruns queued urgent-response products once precise orbits exist and places them in the
# archive's HISTORIC_EVENTS/ (aria_coseis.precise), then pushes its email notifications
# (active_jobs/outbox/) for the GitHub Action to send.
#
# Environment (see aria_coseis.config): COSEIS_ARCHIVE_DIR (the product archive),
# COSEIS_DATA_DIR (topsApp's working folders), COSEIS_PRECISE_CORES (default 17-33),
# COSEIS_PRECISE_LOCK_FILE (default /tmp/coseis_precise.lock).

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
LOG_FILE="logs/precise.log"

# The reruns (they take their own lock and exit at once if a previous run is still going). The
# code is the clone's main, which run_coseis_forward.sh keeps up to date.
python -u -m aria_coseis --precise >> "$LOG_FILE" 2>&1

# Notifications for the email runner: move them into the tracker folder and push, holding the
# git lock shared with run_coseis_forward.sh
STAGING="${COSEIS_OUTBOX_STAGING_DIR:-outputs/outbox_pending}"
OUTBOX="${COSEIS_OUTBOX_DIR:-active_jobs/outbox}"
shopt -s nullglob
NOTES=("$STAGING"/*.json)
if [ ${#NOTES[@]} -eq 0 ]; then
    exit 0
fi

export GIT_TERMINAL_PROMPT=0
exec 9> "$REPO_DIR/.git/coseis-git.lock"
if ! flock -w 900 9; then
    echo "$(date): git is busy; ${#NOTES[@]} notifications stay staged for the next run." >> "$LOG_FILE"
    exit 0
fi

mkdir -p "$OUTBOX"
mv "${NOTES[@]}" "$OUTBOX"/
git add "$OUTBOX"
git commit -q -m "Precise-orbit reruns: notifications for the email runner [skip ci]" >> "$LOG_FILE" 2>&1
source "$REPO_DIR/ops/push_main.sh"
