import sys
import xml.etree.ElementTree as ET

# Check if the required arguments are provided
if len(sys.argv) != 4 or sys.argv[2] != '--property':
    print("Usage: python extractXML.py file.xml --property property_name")
    sys.exit(1)

# Get the filename and property name from the command-line arguments
filename = sys.argv[1]
property_name = sys.argv[3].lower()  # Convert to lowercase for case-insensitive matching

try:
    # Load the XML file
    tree = ET.parse(filename)
    root = tree.getroot()

    # Initialize variable to store the found property value
    value = None

    # Search for the property element with a case-insensitive match on the name attribute
    for elem in root.findall(".//property"):
        # Convert the name attribute to lowercase and compare
        if elem.get('name') and elem.get('name').lower() == property_name:
            value_elem = elem.find('value')
            if value_elem is not None:
                value = value_elem.text
                break  # Stop after finding the first matching property

    # Print the extracted property value
    if value:
        print(f"{property_name.capitalize()}: {value}")
    else:
        print(f"Property '{property_name}' not found.")

except FileNotFoundError:
    print(f"File '{filename}' not found.")
except ET.ParseError:
    print(f"Error parsing '{filename}'. Please check if it's a valid XML file.")

