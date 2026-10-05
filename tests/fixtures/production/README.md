Tracker and partial-job files copied from production commits on `main`, used as test inputs
and as references for what `main` actually produced.

- `tamarindo_ff08a73/`: written by the GitHub Action that discovered us6000tymj (2026-10-01 00:10Z).
- `ende_98ea2b5/`: Ende D61 reset to AWAITING_POST_SEISMIC ("Rerun m7.7 indonesia D61"). The partial
  job file is from before its deletion in d8729e5.
- `ende_92edeb3/`: the same track after the local cron processed it (READY_FOR_EMAIL).
