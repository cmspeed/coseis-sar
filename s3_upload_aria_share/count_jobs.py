import json
import sys
from pathlib import Path

def count_json_items(json_path):
    path = Path(json_path)
    if not path.exists():
        raise FileNotFoundError(f"No such file: {json_path}")
    
    with path.open() as f:
        data = json.load(f)

    if not isinstance(data, list):
        raise ValueError("Top-level JSON object is not a list")
    
    return len(data)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <json_path>")
        sys.exit(1)

    json_path = sys.argv[1]
    total = count_json_items(json_path)
    print(f"Total items: {total}")