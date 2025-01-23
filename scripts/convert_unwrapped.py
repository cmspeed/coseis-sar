from osgeo import gdal
import numpy as np

def convert_unwrapped(unwrapped_pix):

    # Open the input raster
    unwrapped = gdal.Open(unwrapped_pix)

    # Read the amplitude and phase arrays
    amp_array = unwrapped.GetRasterBand(1).ReadAsArray()
    phase_array = unwrapped.GetRasterBand(2).ReadAsArray()

    # Convert phase to meters
    phase_m = phase_array * 0.05546576/ (4 * np.pi)
    
    # Write the converted arrays to a new raster
    driver = gdal.GetDriverByName("ISCE")
    outfile = "merged/testing/filt_topophase_m.unw.geo"
    vrtfile = "merged/testing/filt_topophase_m.unw.geo.vrt"
   
    # Create the new raster
    unwrapped_m = driver.Create(outfile, unwrapped.RasterXSize, unwrapped.RasterYSize, 2, gdal.GDT_Float32)
   
    # Get the raster bands
    band1 = unwrapped_m.GetRasterBand(1)
    band2 = unwrapped_m.GetRasterBand(2)
   
    # Set the NoData value for both bands before writing data
    band1.SetNoDataValue(0)
    band2.SetNoDataValue(0)
   
    # Write the amp_array to band 1 and phase_m to band 2
    band1.WriteArray(amp_array)
    band2.WriteArray(phase_m)

    # Ensure georeferencing information is same as input
    unwrapped_m.SetGeoTransform(unwrapped.GetGeoTransform())
    unwrapped_m.SetProjection(unwrapped.GetProjection())

    # Flush data to disk
    band1.FlushCache()
    band2.FlushCache()
    unwrapped_m.FlushCache()
    unwrapped_m = None 
   
    # Create a .vrt file
    gdal.BuildVRT(vrtfile, [outfile])
    
    return

def main():
    convert_unwrapped('merged/filt_topophase.unw.geo')

if __name__ == "__main__":
    main()