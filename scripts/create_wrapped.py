#!/usr/bin/env python3

'''
Description: Script to create wrapped geotiff file
./create_wrapped.py filt_topophase_watermsk.unw.geo filt_topophase_watermsk.tiff
'''

import os
import sys
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from PIL import Image

def get_cmap(mapname, N, clim):
    '''
    Get the colormap from matplotlib.
    '''

    try:
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        import matplotlib.cm as cmx
    except ImportError:
        raise Exception('Matplotlib is needed if user-defined color table is not provided.')

    cmap = plt.get_cmap(mapname)
    cNorm = colors.Normalize(vmin = clim[0], vmax = clim[1])
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
   
    vals = np.linspace(clim[0], clim[1], endpoint=True)

    outname = mapname + '.cpt'
    
    with open(outname, 'w') as fid:
        for val in vals:
            cval = scalarMap.to_rgba(val)
            fid.write('{0} {1} {2} {3} \n'.format(val,int(cval[0]*255), int(cval[1]*255), int(cval[2]*255)))
            
        fid.write('nv 0 0 0 0 \n')

    return outname


if __name__ == '__main__':
    '''
    Main driver.
    '''

    # Load and read masked unw file
    ds = gdal.Open(str(sys.argv[1]),gdal.GA_ReadOnly)
    geoTrans = ds.GetGeoTransform()
    proj = ds.GetProjection() 
    amp = ds.GetRasterBand(1).ReadAsArray()
    unw = ds.GetRasterBand(2).ReadAsArray() #* 5.5 / (4*np.pi)
    ds = None

    # Calculate wrapped 
    wrp = np.angle(np.exp(1j * unw))

    # Plot
    fig,ax = plt.subplots(1,2,figsize=(12,6), sharey=True)
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    img = ax[0].imshow(unw, cmap=cm.jet)
    ax[0].set_title('Unwrapped Phase')
    cbar = fig.colorbar(img, cax=cax, orientation='vertical')
    divider = make_axes_locatable(ax[1])
    cax = divider.append_axes("right", size="5%", pad=0.05)
    img = ax[1].imshow(wrp, cmap=cm.jet)
    #img.set_clim([-3.14,3.14])
    ax[1].set_title('Wrapped')
    cbar = fig.colorbar(img, cax=cax, orientation='vertical')
    fig.savefig('tmp_wrapped.png')

    # Save as TIFF format
    [ny,nx] = np.shape(wrp)
    tifoutname = str(sys.argv[2])
    drv = gdal.GetDriverByName('GTiff').Create(tifoutname, nx, ny, 1, gdal.GDT_Float32)
    drv.SetGeoTransform(geoTrans)
    drv.SetProjection(proj)
    drv.GetRasterBand(1).WriteArray(wrp)
    drv.GetRasterBand(1).SetNoDataValue(0)
    drv = None

    # Save to Rendered GeoTiff color
    cmap = get_cmap('jet', 64, [-3.14159, 3.14159])

    ###Build colored geotiff
    tifoutname_colored = os.path.splitext(tifoutname)[0]+'_colored.tiff'
    if os.path.exists(tifoutname_colored):
        print('TIF Colored file already exists. Cleaning it ....')
        os.remove(tifoutname_colored)
    outvrtname = gdal.DEMProcessing('',tifoutname,'color-relief', colorFilename = cmap, format = 'MEM', addAlpha =True, options = ['-b', str(1), '-of','VRT'])
    ds = gdal.Translate(tifoutname_colored,outvrtname)
    ds = None

    # kmloptions = gdal.TranslateOptions(gdal.ParseCommandLine("-of KMLSUPEROVERLAY -co format=png"))
    # kmloutname = os.path.splitext(tifoutname)[0]+'_colored.kml'
    # print('Creating KML file:', kmloutname)
    # gdal.Translate(kmloutname, tifoutname_colored, options = kmloptions)






