"""Notifications: email, HTML tables, interactive frame maps and satellite overpass predictions."""
import os
import folium
import sys
import yagmail
from types import SimpleNamespace

from aria_coseis import config


def send_email(subject, body, recipients=None):
    """
    Send an email with the earthquake information to a specified list of recipients.
    :param subject: Email subject
    :param body: Email body content
    :param recipients: List of email addresses to send the email to
    """
    if recipients is None:
        recipients = config.PRIMARY_RECIPIENTS

    if not recipients:
        print("WARNING: Cannot send email. No recipients configured in environment variables.")
        return

    GMAIL_USER = os.environ.get('GMAIL_USER')
    GMAIL_PSWD = os.environ.get('GMAIL_APP_PSWD')

    if not GMAIL_USER or not GMAIL_PSWD:
        print("WARNING: Cannot send email. GMAIL_USER or GMAIL_APP_PSWD environment variables are missing.")
        return

    yag = yagmail.SMTP(GMAIL_USER, GMAIL_PSWD)

    yag.send(
             bcc=recipients,
             subject=subject,
             contents=[body]
             )
    return


def ascii_table_to_html(ascii_table):
    """Converts a raw ASCII table string into a styled HTML table."""
    if not ascii_table or "+" not in ascii_table:
        return f"<pre style='background:#2b2b2b; color:#ccc; padding:10px;'>{ascii_table}</pre>"
    
    lines = ascii_table.strip().split('\n')
    html = '<table style="width: 100%; border-collapse: collapse; font-size: 13px; text-align: left; margin-bottom: 20px;">'
    
    is_header = True
    for line in lines:
        if line.strip().startswith('+'):
            continue  # Skip the border lines
        if line.strip().startswith('|'):
            # Extract cell content, ignoring the first and last empty splits from the outer '|'
            cells = [cell.strip() for cell in line.split('|')[1:-1]]
            html += "<tr>"
            for cell in cells:
                if is_header:
                    html += f'<th style="border-bottom: 2px solid #0055a4; padding: 10px 8px; background-color: #f0f4f8; color: #003366;">{cell}</th>'
                else:
                    html += f'<td style="border-bottom: 1px solid #e0e0e0; padding: 8px; color: #444;">{cell}</td>'
            html += "</tr>"
            is_header = False
            
    html += "</table>"
    return html


def make_interactive_map(frame_dataframe, title, coords, url):
        
    # Extract latitude and longitude from coords
    lon, lat = coords[0], coords[1]
    
    # Generate an interactive map
    map_object = frame_dataframe.explore(column="pathNumber", 
                                            cmap="viridis", 
                                            tooltip=["flightDirection", "pathNumber", "frameNumber", "startTime"],
                                            )

    # Add a basemap explicitly with Folium (Esri World Imagery)
    folium.TileLayer('Esri World Imagery').add_to(map_object)

    # Create a popup with both the title and a clickable URL
    popup_content = f"""
    <b>{title}</b><br>
    <a href="{url}" target="_blank">{url}</a>
    """
    # Add AOI centroid marker with a popup that includes the title and URL
    folium.Marker(
        location=[lat, lon],
        popup=folium.Popup(popup_content, max_width=300),
        icon=folium.Icon(color='red', icon='info-sign')
    ).add_to(map_object)

    # Save to an HTML file
    map_filename = f"{title}_SLC_Map.html"
    map_object.save(map_filename)
    
    return map_filename


def get_next_pass(AOI, timestamp_dir, satellite="sentinel-1"):
    """
    Get the next satellite pass over the given AOI.
    Uses the next_pass.py script to determine the next overpass.
    :param AOI: Shapely Polygon object representing the AOI
    :param satellite: Satellite name ('sentinel-1', 'sentinel-2', or 'landsat')
    :return: Next overpass time or None if an error occurs
    """
    import os
    import sys
    import next_pass
    next_pass_dir = os.path.dirname(next_pass.__file__)
    if next_pass_dir not in sys.path:
        sys.path.append(next_pass_dir)
    try:
        from next_pass import plot_maps
    except ImportError as e:
        print(f"Could not import plot_maps: {e}")
        return None, None, None
    
    from datetime import date

    min_lon, min_lat, max_lon, max_lat = AOI.bounds
    bbox = [str(min_lat), str(max_lat), str(min_lon), str(max_lon)]

    print("=========================================")
    print(f"Querying next-pass for all satellites over AOI...")
    print("=========================================")
    
    args = SimpleNamespace(bbox=bbox, sat="all", event_date=date.today(), look_back=14, cloudiness=False)

    try:
        result = next_pass.find_next_overpass(args, timestamp_dir)
    except Exception as e:
        print(f"Next pass error: {e}")
        return None, None, None
    
    result_s1 = result.get("sentinel-1") 
    result_nisar = result.get("nisar")
    
    s1_info = result_s1.get("next_collect_info", "No S1 info available.") if result_s1 else "No S1 info."
    nisar_info = result_nisar.get("next_collect_info", "No NISAR info available.") if result_nisar else "No NISAR info."
    
    default_map_file = timestamp_dir / "satellite_overpasses_map.html"
    
    # Generate Joint S1 + NISAR Map
    map_path = timestamp_dir / "overpass_map.html"
    try:
        plot_maps.make_overpasses_map(result_s1, None, None, result_nisar, args.bbox, timestamp_dir)
        if default_map_file.exists():
            os.rename(default_map_file, map_path)
    except Exception as e:
        print(f"Could not generate overpass map: {e}")
        map_path = None

    return s1_info, nisar_info, map_path
