#!/bin/bash
# Forward-mode processing run for the processing machine's cron job, e.g.:
#   5,15 * * * * /path/to/repo/ops/run_coseis_forward.sh
# Checks tracked events for post-event Sentinel-1 data, runs topsApp, and commits the
# updated tracker (active_jobs/, including active_jobs/partials/) to main.
#
# Optional environment overrides (see aria_coseis.config): COSEIS_DATA_DIR for topsApp's working
# folders (default <repo>/data), COSEIS_ARCHIVE_DIR for the product archive (products go to its
# URGENT_RESPONSE/), COSEIS_LOCK_FILE (default /tmp/coseis_processing.lock).

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

# Check the lock file BEFORE doing anything with Git (same default as aria_coseis.config.LOCK_FILE).
# The lock holds the owning process ID; if that process is gone (e.g. a killed run), the lock is
# stale and is removed. A legacy lock without a process ID counts as stale after 12 hours.
LOCK_FILE="${COSEIS_LOCK_FILE:-/tmp/coseis_processing.lock}"
if [ -f "$LOCK_FILE" ]; then
    LOCK_PID="$(awk 'NR==1 {print $1}' "$LOCK_FILE")"
    if [[ "$LOCK_PID" =~ ^[0-9]+$ ]] && ! ps -p "$LOCK_PID" > /dev/null 2>&1; then
        STALE=yes
    elif [[ ! "$LOCK_PID" =~ ^[0-9]+$ ]] && [ -n "$(find "$LOCK_FILE" -mmin +720 2> /dev/null)" ]; then
        STALE=yes
    else
        STALE=no
    fi
    if [ "$STALE" = yes ]; then
        echo "$(date): Removing stale lock file $LOCK_FILE (the run that created it is no longer active)." >> "$LOG_FILE"
        rm -f "$LOCK_FILE"
    else
        echo "$(date): Previous processing run still active. Bash script exiting." >> "$LOG_FILE"
        exit 0
    fi
fi

# Git runs unattended here: fail instead of waiting for a password prompt, and send its
# output (including errors) to the log
export GIT_TERMINAL_PROMPT=0

# Two cron jobs use git in this clone (this one and run_coseis_precise.sh): each git step holds
# this lock, so they never run git at the same time. Processing itself runs without it (9>&-).
exec 9> "$REPO_DIR/.git/coseis-git.lock"

# Sync with Github
if flock -w 900 9; then
    git checkout -q main >> "$LOG_FILE" 2>&1

    # Pull the latest tracker state that the GitHub Action just updated (keeps any local commits
    # left from an earlier failed push on top). On a conflict, abort and leave things as they were.
    if ! git pull -q --rebase origin main >> "$LOG_FILE" 2>&1; then
        git rebase --abort 2> /dev/null || true
        echo "$(date): git pull failed; processing anyway, will retry the push below." >> "$LOG_FILE"
    fi
    flock -u 9
else
    echo "$(date): git is busy (precise runner); processing without pulling." >> "$LOG_FILE"
fi

# Execute local processing (discovery disabled via --process_only)
python -m aria_coseis --forward --pairing coseismic --resolution 30 --do_processing --process_only >> "$LOG_FILE" 2>&1 9>&-

if ! flock -w 900 9; then
    echo "$(date): git is busy (precise runner); the next run will commit and push." >> "$LOG_FILE"
    exit 0
fi

# Sync with Github: stage tracker edits, new files and deletions (finished partials) under active_jobs/
git add -A active_jobs/ || true

if ! git diff --cached --quiet; then
    git commit -q -m "Local processing: update COSEIS tracking state and remove finished partials [skip ci]" >> "$LOG_FILE" 2>&1
else
    echo "No processing completed this run; tracking state unchanged." >> "$LOG_FILE"
fi

# Push any local tracker commits, including ones left from an earlier failed push. If the GitHub
# Action pushed in between, rebase onto it and try again.
source "$REPO_DIR/ops/push_main.sh"
