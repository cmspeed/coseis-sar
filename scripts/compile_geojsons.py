import json
import glob

def main():
    features = []

    for filepath in glob.glob("m_*.geojson"):
        with open(filepath) as f:
            data = json.load(f)
            if data["type"] == "FeatureCollection":
                features.extend(data["features"])
            elif data["type"] == "Feature":
                features.append(data)
            else:
                print(f"Skipping unsupported GeoJSON type in {filepath}")

    combined = {
        "type": "FeatureCollection",
        "features": features
    }

    with open("all_jobs.geojson", "w") as out_f:
        json.dump(combined, out_f, indent=2)

    print(f"Combined {len(features)} features from {len(glob.glob('m_*.geojson'))} file(s) into all_jobs.geojson")

if __name__ == "__main__":
    main()
