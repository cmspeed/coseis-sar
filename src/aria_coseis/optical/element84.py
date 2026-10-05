"""Sentinel-2 scene search and pairing via the Element84 Earth Search STAC API."""
import re
from shapely.geometry import mapping, shape
from datetime import datetime, timezone
from collections import defaultdict

from aria_coseis.config import OPTICAL_CLOUD_THRESHOLD
from aria_coseis.optical.jobs import get_utm_zone, make_optical_job_json
from aria_coseis.utils import convert_time


def search_element84_stac(aoi_polygon, start_date, end_date):
    """
    Searches Element84 Earth Search STAC API v1 for Sentinel-2 L2A COGs using pystac_client.
    """
    from pystac_client import Client

    print("Searching Element84 Earth Search for Sentinel-2 L2A...")
    
    # Connect to the API endpoint
    api_url = "https://earth-search.aws.element84.com/v1"
    client = Client.open(api_url)
    
    # Define parameters
    collection_sentinel_2_l2a = "sentinel-2-c1-l2a" 
    intersects_geom = mapping(aoi_polygon)
    
    try:
        # Perform the search using the client
        search = client.search(
            collections=[collection_sentinel_2_l2a],
            intersects=intersects_geom,
            datetime=f"{start_date}/{end_date}",
            query={"eo:cloud_cover": {"lte": OPTICAL_CLOUD_THRESHOLD}},
            limit=100
        )
        
        # Retrieve the metadata items
        items = list(search.items())
        print(f"  Found {len(items)} Sentinel-2 products.")
        
        granules = []
        
        # Each item contains information about the scene geometry, acquisition time, and properties
        for item in items:
            props = item.properties
            
            platform = props.get('platform', 'sentinel-2')
            cloud_cover = props.get('eo:cloud_cover', 0.0)
            
            # Extract Tile ID from grid:code (e.g., 'MGRS-35SNA' -> '35SNA')
            grid_code = props.get('grid:code') or props.get('s2:mgrs_tile', '')
            tile_id = grid_code.replace('MGRS-', '') if grid_code else ''
            
            # Extract UTM Zone directly from mgrs:utm_zone
            utm_zone = str(props.get('mgrs:utm_zone', 'Unknown'))
            
            # Extract Orbit ID from s2:product_uri (e.g., '..._R107_...')
            product_uri = props.get('s2:product_uri', '')
            orbit_match = re.search(r'_R(\d{3})_', product_uri)
            orbit_id = orbit_match.group(1) if orbit_match else "000"
            
            granules.append({
                "granule_id": item.id,
                "date": props.get('datetime', '').split('T')[0],
                "datetime": props.get('datetime'),
                "cloud_cover": cloud_cover,
                "platform": platform,
                "footprint": shape(item.geometry),
                "orbit_id": orbit_id,
                "tile_id": tile_id,
                "utm_zone": utm_zone
            })
            
        return granules

    except Exception as e:
        print(f"  !! Error searching Element84 STAC with pystac_client: {e}")
        return []


def process_candidate_group(orbit_key, dates_dict, rupture_dt, aoi_polygon, title, event_id, aoi_area, strategy_name, role=None):
    """
    Generic logic to select best Pre/Post pair from a grouped dictionary of candidates. Used by both 'Dominant' and 'Split' strategies.
    :param orbit_key: The key representing the group (e.g., "047" for Orbit-based, "047_Z46" for Orbit_Zone-based)
    :param dates_dict: Dictionary where keys are date strings and values are lists of scene dictionaries for that date
    :param rupture_dt: Datetime object representing the earthquake's origin time
    :param aoi_polygon: Shapely Polygon representing the Area of Interest (used for coverage calculation)
    :param title: Title of the earthquake event (used for job naming)
    :param event_id: The earthquake event ID (e.g., "us7000dflf") for job naming
    :param aoi_area: Area of the AOI polygon (used for coverage calculation)
    :param strategy_name: Name of the strategy ("Dominant" or "Split") for logging purposes
    :param role: Optional parameter to indicate if this group is 'dominant' or 'minority' in the context of mixed zones (used for logging and feature properties)
    :return: A tuple containing a list of job JSON objects and a list of GeoJSON features for visualization.
    """
    pre_candidates = []
    post_candidates = []
    
    # Evaluate candidates for every available date
    for date_str, scenes in dates_dict.items():
        if not scenes: continue
        
        date_dt = datetime.strptime(date_str, "%Y-%m-%d").replace(tzinfo=timezone.utc)
        platform = scenes[0]['platform']

        # Deduplicate tiles for this date/group
        unique_scenes = []
        seen_tiles = set()
        scenes.sort(key=lambda x: x['cloud_cover'])
        
        combined_poly = None
        for s in scenes:
            # Create a unique key for the tile (UTM + Triplet, e.g., '46QHK')
            tile_id = s.get('tile_id')
            # Check if tile_id exists and hasn't been processed yet
            if tile_id and tile_id not in seen_tiles:
                unique_scenes.append(s)
                seen_tiles.add(tile_id)
                if s.get('footprint'):
                    if combined_poly is None: 
                        combined_poly = s['footprint']
                    else: 
                        combined_poly = combined_poly.union(s['footprint'])

        if not unique_scenes:
            continue

        # Calc coverage
        coverage_pct = 0.0
        if combined_poly and aoi_polygon:
            try:
                clean_poly = combined_poly.buffer(0)
                intersection = clean_poly.intersection(aoi_polygon)
                coverage_pct = (intersection.area / aoi_area) * 100.0
            except: pass
        
        avg_cc = sum(s['cloud_cover'] for s in unique_scenes) / len(unique_scenes)
        delta = abs((date_dt - rupture_dt).days)
        
        candidate = {
            'date': date_str, 'scenes': unique_scenes, 'platform': platform,
            'coverage': coverage_pct, 'tile_count': len(unique_scenes), 'cc': avg_cc, 'delta': delta
        }
        
        if date_dt < rupture_dt: pre_candidates.append(candidate)
        elif date_dt > rupture_dt: post_candidates.append(candidate)

    jobs = []
    features = []

    # Select Best Pair
    if not pre_candidates or not post_candidates:
        return [], []

    # Sort: Max Coverage > Low Cloud > Low Delta
    pre_candidates.sort(key=lambda x: (-x['coverage'], x['cc'], x['delta']))
    best_pre = pre_candidates[0]
    
    post_candidates.sort(key=lambda x: (-x['coverage'], x['cc'], x['delta']))
    best_post = post_candidates[0]
    
    print(f"  [{strategy_name}] {orbit_key}: Pre={best_pre['date']} ({best_pre['coverage']:.1f}% Cov), Post={best_post['date']} ({best_post['coverage']:.1f}% Cov)")

    # Generate Job
    primary_ids = [s['granule_id'].replace('.SAFE', '') for s in best_pre['scenes']]
    secondary_ids = [s['granule_id'].replace('.SAFE', '') for s in best_post['scenes']]
    
    # Handle naming differences based on strategy
    zone_suffix = None
    orbit_id = orbit_key
    if "Z" in orbit_key:
        parts = orbit_key.split('_')
        orbit_id = parts[0]
        zone_suffix = parts[1]

    job = make_optical_job_json(title, event_id, orbit_id, best_pre['date'], best_post['date'], primary_ids, secondary_ids, pre_cc=best_pre['cc'], post_cc=best_post['cc'], zone_suffix=zone_suffix)
    jobs.append(job)

    # Generate GeoJSON Features for Visualization
    for stage, candidate in [("Pre-Event", best_pre), ("Post-Event", best_post)]:
        for scene in candidate['scenes']:
            if scene.get('footprint'):
                feat = {
                    "type": "Feature",
                    "geometry": mapping(scene['footprint']),
                    "properties": {
                        "job_name": job['name'],
                        "strategy": strategy_name,
                        "role": role,
                        "timing": stage,
                        "date": candidate['date'],
                        "orbit": orbit_id,
                        "zone": get_utm_zone(scene['granule_id']),
                        "granule_id": scene['granule_id']
                    }
                }
                features.append(feat)
                
    return jobs, features


def check_mixed_zones_in_group(dates_dict):
    """
    Checks if any single date/acquisition within an orbit group contains  tiles from multiple UTM zones. 
    This indicates a 'Mixed' case that requires Split/Dominant strategy comparison.
    :param dates_dict: Dictionary where keys are date strings and values are lists of scene dictionaries for that date
    :return: True if mixed zones are found within any date, False otherwise
    """
    for date, scenes in dates_dict.items():
        zones = set(get_utm_zone(s['granule_id']) for s in scenes)
        if len(zones) > 1:
            return True
    return False


def find_optical_pairs_element84(optical_scenes, rupture_time, title, event_id, aoi_polygon):
    """
    Generates optical pairs using two strategies concurrently for comparison:
    1. DOMINANT: Enforces one UTM zone per Orbit (drops minority tiles).
    2. SPLIT: Creates separate jobs for every UTM zone found.

    Uses the Element84 file nomenclature.
    
    Returns:
    - standard_jobs (Dominant jobs for all orbits, for general processing)
    - mixed_dom_features (Features for Dominant strategy, ONLY for mixed cases)
    - mixed_split_features (Features for Split strategy, ONLY for mixed cases)
    """
    rupture_dt = convert_time(rupture_time)
    aoi_area = aoi_polygon.area

    # Group by Orbit -> Date -> List of Scenes
    orbit_groups = defaultdict(lambda: defaultdict(list))
    for scene in optical_scenes:
        if scene.get('orbit_id'):
            orbit_groups[scene['orbit_id']][scene['date']].append(scene)

    # Standard production output (Dominant strategy applied to ALL)
    all_prod_jobs = []
    
    # Comparison/Debug output (ONLY for orbits that actually have mixed zones)
    mixed_dom_feats = []
    mixed_split_feats = []

    print(f"Processing {len(orbit_groups)} tracks with DOMINANT strategy (and SPLIT where applicable)...")

    for orbit_id, dates_dict in orbit_groups.items():
        
        # Check if this orbit even has a mixing issue
        is_mixed = check_mixed_zones_in_group(dates_dict)
        
        # --- STRATEGY A: DOMINANT ZONE (Applied to ALL orbits) ---
        # Note: Strategy A filters to keep only the dominant data per date.
        dom_dates_dict = defaultdict(list)
        for d, scenes in dates_dict.items():
            zone_counts = defaultdict(int)
            for s in scenes: zone_counts[get_utm_zone(s['granule_id'])] += 1
            if zone_counts:
                winner = max(zone_counts, key=zone_counts.get)
                dom_dates_dict[d] = [s for s in scenes if get_utm_zone(s['granule_id']) == winner]
        
        # Generate Dominant Jobs (Always added to production list)
        j_dom, f_dom = process_candidate_group(orbit_id, dom_dates_dict, rupture_dt, aoi_polygon, title, event_id, aoi_area, "DOMINANT", role="dominant")
        all_prod_jobs.extend(j_dom)

        # If this was a Mixed case, add to our specific debug lists
        if is_mixed:
            mixed_dom_feats.extend(f_dom)

            # --- STRATEGY B: SPLIT ZONES (Only needed for Mixed cases) ---
            # Determine the Global Dominant Zone (for the whole orbit, not just per date)
            # This allows us to label the split results as 'dominant' or 'minority'
            all_orbit_scenes = [s for scenes in dates_dict.values() for s in scenes]
            global_zone_counts = defaultdict(int)
            for s in all_orbit_scenes:
                global_zone_counts[get_utm_zone(s['granule_id'])] += 1
            
            global_winner = None
            if global_zone_counts:
                global_winner = max(global_zone_counts, key=global_zone_counts.get)

            split_groups = defaultdict(lambda: defaultdict(list))
            for d, scenes in dates_dict.items():
                for s in scenes:
                    z = get_utm_zone(s['granule_id'])
                    key = f"{orbit_id}_Z{z}"
                    split_groups[key][d].append(s)
            
            for key, z_dates_dict in split_groups.items():
                # Extract zone from key to determine role
                z_str = key.split('_Z')[-1]
                role = 'dominant' if z_str == global_winner else 'minority'

                _, f_split = process_candidate_group(key, z_dates_dict, rupture_dt, aoi_polygon, title, event_id, aoi_area, "SPLIT", role=role)
                mixed_split_feats.extend(f_split)

    return all_prod_jobs, mixed_dom_feats, mixed_split_feats
