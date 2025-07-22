import json
import sys

def main(json_path):
    with open(json_path, "r") as f:
        data = json.load(f)
    print(f"Number of jobs: {len(data)}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python count_jobs.py job_list.json")
        sys.exit(1)
    main(sys.argv[1])