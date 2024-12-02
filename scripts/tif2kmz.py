# -*- coding: utf-8 -*-
"""
Produce raster ground overlay for Google Earth from GDAL raster file (color/gray interpretation)

\author paryan
\date 2018.05.17
"""
import argparse
import os
import os.path as path
import json
from osgeo import osr, ogr, gdal # min. ver 2.1
import sys
from jinja2 import Template
import zipfile

        
def main(args):
    if args.infile is not None:
        base, ext = path.splitext(args.infile)
        kmltif = base + '_kml' + ext
        kmlpng = base + '_kml' + '.png'
        kmlfile = base + '.kml'
        kmzfile = base + '.kmz'
        
        in_ds = gdal.Open(args.infile, gdal.GA_ReadOnly)
        # warp original dataset to WGS84 used by google earth
        tmp_ds = gdal.Warp(base + '_kml' + ext, in_ds, dstSRS='EPSG:4326')
        
        # get geotransform for latlonbox
        gt = tmp_ds.GetGeoTransform()
        cols = tmp_ds.RasterXSize
        rows = tmp_ds.RasterYSize
        
        # convert image format for use by Icon, nodata will be transparent
        gdal.Translate(kmlpng, tmp_ds, format="PNG")
        in_ds = None
        tmp_ds = None
        
        data = {
            'folder' : { 'name' : base, 'description': args.description },
            'overlay' : { 'name' : base, 'description': args.description },
            'imgurl' : kmlpng,
            'north': gt[3],
            'south': gt[3] + gt[5] * rows,
            'west': gt[0],
            'east': gt[0] + gt[1] * cols,
            'rotation':0,
            'color': "{}ffffff".format(hex(int(args.alpha * 255))[2:])
        }
        
        t = Template("""<?xml version="1.0" encoding="UTF-8"?>
<kml>
    <Folder>
        <name>{{ folder.name }}</name>
        <visibility>1</visibility>
        <description>{{ folder.description }}</description>
        <GroundOverlay>
            <name>{{ overlay.name }}</name>
            <visibility>1</visibility>
            <description>{{ overlay.description }}</description>
            <color>{{ color }} </color>
            <Icon>
              <href>{{imgurl}}</href>
            </Icon>
            <LatLonBox>
              <north>{{north}}</north>
              <south>{{south}}</south>
              <east>{{east}}</east>
              <west>{{west}}</west>
              <rotation>{{rotation}}</rotation>
            </LatLonBox>
        </GroundOverlay>
    </Folder>
</kml>""")

        with open(kmlfile, "wb") as fo:
            fo.write(t.render(data))
            
        if args.zip:
            zipf = zipfile.ZipFile(kmzfile, "w", zipfile.ZIP_DEFLATED)
            zipf.write(kmlfile)
            zipf.write(kmlpng)
            zipf.close()
                
        if args.removetemp:
            if path.exists(kmltif): os.unlink(kmltif)
            if path.exists(kmlfile): os.unlink(kmlfile)
            if path.exists(kmlpng): os.unlink(kmlpng)
            
        
if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-infile', help="input GDAL raster file", required=True)
    parser.add_argument('-outfile', help="output filename, when not given will be estimated")
    parser.add_argument('-zip', action="store_true", help="pack kml and image overlay as kmz file", default=True)
    parser.add_argument('-removetemp', action='store_true', help="remove intermediate files", default=True)
    parser.add_argument('-alpha', type=float, default=1.0, help="alpha transparency")
    parser.add_argument('-description', help="description of the file", default='')
    args = parser.parse_args()
    main(args)
