"""Small shared helpers: name formatting and USGS time conversion."""

from __future__ import annotations

import re
import unicodedata
from datetime import datetime, timezone


def to_snake_case(input_string: str) -> str:
    """
    Convert the given string to snake_case with ASCII-safe characters.
    Strips accents and transliterates Unicode to closest ASCII.
    """
    # Normalize and transliterate Unicode characters to closest ASCII equivalent
    normalized = (
        unicodedata.normalize("NFKD", input_string).encode("ascii", "ignore").decode("ascii")
    )

    # Replace non-alphanumeric characters with spaces
    cleaned_string = re.sub(r"[^\w\s]", "", normalized)

    # Replace spaces with underscores and convert to lowercase
    snake_case_string = re.sub(r"\s+", "_", cleaned_string.strip()).lower()

    return snake_case_string


def convert_time(time: int | float) -> datetime:
    """
    Convert the given Unix timestamp in milliseconds to a UTC datetime object.
    :param time: Unix timestamp in milliseconds (int or float)
    :return: timezone-aware UTC datetime, truncated to whole seconds
    """
    timestamp_s = time / 1000  # Convert milliseconds to seconds
    dt = datetime.fromtimestamp(timestamp_s, tz=timezone.utc)  # Convert to datetime object in UTC
    dt = dt.replace(microsecond=0)  # Remove microseconds
    return dt
