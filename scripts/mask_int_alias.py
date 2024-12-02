#!/usr/bin/env python3

'''
Description: Script to mask areas with high rate of fringes using connected components
./mask_int_alias.py filt_topophase_cut.unw.conncomp.geo filt_topophase_cut_msk.unw.geo filt_topophase_cut_mskAlias.unw.geo
'''

import os
import sys
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Load and read connected components
ds = gdal.Open(str(sys.argv[1]),gdal.GA_ReadOnly)
conn = ds.GetRasterBand(1).ReadAsArray()
ds = None

# Load and read masked unw file
ds = gdal.Open(str(sys.argv[2]),gdal.GA_ReadOnly)
geoTrans = ds.GetGeoTransform()
proj = ds.GetProjection() 
amp = ds.GetRasterBand(1).ReadAsArray()
unw = ds.GetRasterBand(2).ReadAsArray() #* 5.5 / (4*np.pi)
ds = None

# Mask
msk = np.array(conn)*-1
msk[msk<0] = 1

# Mask Unwrapped File
unwMsk = unw*msk

# Plot
fig,ax = plt.subplots(1,3,figsize=(12,6), sharey=True)
divider = make_axes_locatable(ax[0])
cax = divider.append_axes("right", size="5%", pad=0.05)
img = ax[0].imshow(msk, cmap=cm.jet)
ax[0].set_title('Connected Components')
cbar = fig.colorbar(img, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax[1])
cax = divider.append_axes("right", size="5%", pad=0.05)
img = ax[1].imshow(unw, cmap=cm.jet)
ax[1].set_title('Unwrapped Phase')
cbar = fig.colorbar(img, cax=cax, orientation='vertical')
divider = make_axes_locatable(ax[2])
cax = divider.append_axes("right", size="5%", pad=0.05)
img = ax[2].imshow(unwMsk, cmap=cm.jet)
#img.set_clim([0,1])
ax[2].set_title('Masked')
cbar = fig.colorbar(img, cax=cax, orientation='vertical')
plt.show()

# Save as ISCE format
[ny,nx] = np.shape(conn)
drv = gdal.GetDriverByName('ISCE').Create(str(sys.argv[3]), nx, ny, 2, gdal.GDT_Float32, options=["SCHEME=BIL"] )
drv.SetGeoTransform(geoTrans)
drv.SetProjection(proj)
drv.GetRasterBand(1).WriteArray(amp*msk)
drv.GetRasterBand(2).WriteArray(unwMsk)
drv = None
