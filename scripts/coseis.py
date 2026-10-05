"""
Command-line entry point used by cron and GitHub Actions (`cd scripts && python coseis.py ...`).
The code lives in the aria_coseis package under src/.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from aria_coseis.cli import main  # noqa: E402

if __name__ == "__main__":
    main()
