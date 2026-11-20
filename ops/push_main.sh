# Sourced by the cron wrappers, with the git lock held and LOG_FILE set: push any local commits on
# main, including ones left from an earlier failed push. If someone else (the GitHub Action, or
# the other cron job) pushed in between, rebase onto it and try again.
for ATTEMPT in 1 2 3; do
    git fetch -q origin main >> "$LOG_FILE" 2>&1
    if [ -z "$(git rev-list origin/main..HEAD)" ]; then
        break
    fi
    if git pull -q --rebase origin main >> "$LOG_FILE" 2>&1 && git push -q origin main >> "$LOG_FILE" 2>&1; then
        break
    fi
    git rebase --abort 2> /dev/null || true
    echo "$(date): Push attempt $ATTEMPT failed; retrying." >> "$LOG_FILE"
    sleep 30
done
if [ -n "$(git rev-list origin/main..HEAD)" ]; then
    echo "$(date): WARNING: local commits are not on GitHub yet; the next run will retry." >> "$LOG_FILE"
fi
