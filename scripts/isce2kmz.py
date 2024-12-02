import argparse
import numpy as np
from osgeo import gdal
import os
import zipfile

def read_isce_data(file_path):
    """Read ISCE data file (assuming it's a GeoTIFF) and return as numpy array."""
    dataset = gdal.Open(file_path)
    if dataset is None:
        raise FileNotFoundError(f"Could not open file: {file_path}")

    band = dataset.GetRasterBand(1)
    array = band.ReadAsArray()
    geo_transform = dataset.GetGeoTransform()
    
    return array, geo_transform

def create_kml_from_array(array, geo_transform):
    """Create a KML file content from a numpy array."""
    kml_content = '''<?xml version="1.0" encoding="UTF-8"?>
    <kml xmlns="http://www.opengis.net/kml/2.2">
        <Document>
            <name>ISCE Data</name>
        '''
    
    rows, cols = array.shape
    ulx = geo_transform[0]
    uly = geo_transform[3]
    x_size = geo_transform[1]
    y_size = geo_transform[5]

    for row in range(rows):
        for col in range(cols):
            value = array[row, col]
            lon = ulx + col * x_size
            lat = uly + row * y_size
            kml_content += f'''
            <Placemark>
                <name>{value}</name>
                <Point>
                    <coordinates>{lon},{lat},0</coordinates>
                </Point>
            </Placemark>
            '''
    
    kml_content += '''
        </Document>
    </kml>
    '''
    
    return kml_content

def create_kmz(kml_content, output_kmz):
    """Create a KMZ file from KML content."""
    kml_file = 'output.kml'
    with open(kml_file, 'w') as f:
        f.write(kml_content)
    
    with zipfile.ZipFile(output_kmz, 'w') as kmz:
        kmz.write(kml_file, arcname=os.path.basename(kml_file))
    
    os.remove(kml_file)  # Clean up the temporary KML file

def convert_isce_to_kmz(isce_file, output_kmz):
    array, geo_transform = read_isce_data(isce_file)
    kml_content = create_kml_from_array(array, geo_transform)
    create_kmz(kml_content, output_kmz)

def main():
    parser = argparse.ArgumentParser(description="Convert ISCE data to KMZ for Google Earth.")
    parser.add_argument('-i', '--input', required=True, help='Path to the ISCE data file (GeoTIFF)')
    parser.add_argument('-o', '--output', required=True, help='Path to the output KMZ file')
    
    args = parser.parse_args()
    
    convert_isce_to_kmz(args.input, args.output)
    print(f"KMZ file created: {args.output}")

if __name__ == "__main__":
    main()

