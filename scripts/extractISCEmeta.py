import datetime
import glob
import logging
import os
import shelve
import time
import numpy as np

def get_processor(meta_file):
    """
    Get the name of ISCE processor (imaging mode)
    """
    meta_dir = os.path.dirname(os.path.abspath(meta_file))
    tops_meta_file = os.path.join(meta_dir, 'IW*.xml')
    alos_meta_file = os.path.join(meta_dir, '*.track.xml')
    stripmap_meta_files = [os.path.join(meta_dir, i) for i in ['data.db', 'data.dat', 'data']]

    processor = None
    if len(glob.glob(tops_meta_file)) > 0:
        # topsStack
        processor = 'tops'

    elif len(glob.glob(alos_meta_file)) > 0:
        # alosStack / alos2App
        processor = 'alosStack'

    elif any(os.path.isfile(i) for i in stripmap_meta_files):
        # stripmapStack
        processor = 'stripmap'

    elif meta_file.endswith('.xml'):
        # stripmapApp
        processor = 'stripmap'

    else:
        raise ValueError(f'Un-recognized ISCE processor for metadata file: {meta_file}')
    return processor

def extract_isce_metadata(meta_file, geom_dir=None, rsc_file=None, update_mode=True):
    """Extract metadata from ISCE stack products
    Parameters: meta_file : str, path of metadata file, reference/IW1.xml or referenceShelve/data.dat
                geom_dir  : str, path of geometry directory.
                rsc_file  : str, output file name of ROIPAC format rsc file. None for not write to disk.
    Returns:    meta      : dict
                frame     : object, isceobj.Scene.Frame.Frame / isceobj.Scene.Burst.Burst
    """
    # check existing rsc_file
    if update_mode and ut.run_or_skip(rsc_file, in_file=meta_file, readable=False) == 'skip':
        return readfile.read_roipac_rsc(rsc_file), None

    # 1. read/extract metadata from XML / shelve file
    processor = get_processor(meta_file)
    if processor == 'tops':
        print('extract metadata from ISCE/topsStack xml file:', meta_file)
        meta, frame = extract_tops_metadata(meta_file)

    elif processor == 'alosStack':
        print('extract metadata from ISCE/alosStack xml file:', meta_file)
        meta, frame = extract_alosStack_metadata(meta_file, geom_dir)

    else:
        print('extract metadata from ISCE/stripmapStack shelve file:', meta_file)
        meta, frame = extract_stripmap_metadata(meta_file)

    # 2. extract metadata from geometry file
    if geom_dir:
        if processor != 'alosStack':
            meta = extract_geometry_metadata(geom_dir, meta)

    # 3. common metadata
    meta['PROCESSOR'] = 'isce'
    if 'ANTENNA_SIDE' not in meta.keys():
        meta['ANTENNA_SIDE'] = '-1'

    # convert all value to string format
    for key, value in meta.items():
        meta[key] = str(value)

    # write to .rsc file
    meta = readfile.standardize_metadata(meta)
    if rsc_file:
        print('writing ', rsc_file)
        writefile.write_roipac_rsc(meta, rsc_file)
    return meta, frame

def extract_tops_metadata(xml_file):
    """Read metadata from xml file for Sentinel-1/TOPS
    Parameters: xml_file : str, path of the .xml file, i.e. reference/IW1.xml
    Returns:    meta     : dict, metadata
                burst    : isceobj.Sensor.TOPS.BurstSLC.BurstSLC object
    """
    import isce
    import isceobj
    from isceobj.Planet.Planet import Planet

    obj = load_product(xml_file)
    burst = obj.bursts[0]
    burstEnd = obj.bursts[-1]

    meta = {}
    meta['prf']             = burst.prf
    meta['startUTC']        = burst.burstStartUTC
    meta['stopUTC']         = burstEnd.burstStopUTC
    meta['radarWavelength'] = burst.radarWavelength
    meta['startingRange']   = burst.startingRange
    meta['passDirection']   = burst.passDirection
    meta['polarization']    = burst.polarization
    meta['trackNumber']     = burst.trackNumber
    meta['orbitNumber']     = burst.orbitNumber

    try:
        meta['PLATFORM'] = sensor.standardize_sensor_name(obj.spacecraftName)
    except:
        if os.path.basename(xml_file).startswith('IW'):
            meta['PLATFORM'] = 'sen'

    time_seconds = (burst.sensingMid.hour * 3600.0 +
                    burst.sensingMid.minute * 60.0 +
                    burst.sensingMid.second)
    meta['CENTER_LINE_UTC'] = time_seconds

    orbit = burst.orbit
    peg = orbit.interpolateOrbit(burst.sensingMid, method='hermite')

    # Sentinel-1 TOPS pixel spacing
    Vs = np.linalg.norm(peg.getVelocity())   #satellite speed
    meta['satelliteSpeed'] = Vs
    meta['azimuthPixelSize'] = Vs * burst.azimuthTimeInterval
    meta['rangePixelSize'] = burst.rangePixelSize

    # Sentinel-1 TOPS spatial resolution
    iw_str = 'IW2'
    if os.path.basename(xml_file).startswith('IW'):
        iw_str = os.path.splitext(os.path.basename(xml_file))[0]
    meta['azimuthResolution'] = sensor.SENSOR_DICT['sen'][iw_str]['azimuth_resolution']
    meta['rangeResolution']   = sensor.SENSOR_DICT['sen'][iw_str]['range_resolution']

    elp = Planet(pname='Earth').ellipsoid
    llh = elp.xyz_to_llh(peg.getPosition())
    elp.setSCH(llh[0], llh[1], orbit.getENUHeading(burst.sensingMid))
    meta['HEADING'] = orbit.getENUHeading(burst.sensingMid)
    meta['earthRadius'] = elp.pegRadCur
    meta['altitude'] = llh[2]

    # for Sentinel-1
    meta['beam_mode'] = 'IW'
    meta['swathNumber'] = burst.swathNumber
    # 1. multipel subswaths
    xml_files = glob.glob(os.path.join(os.path.dirname(xml_file), 'IW*.xml'))
    if len(xml_files) > 1:
        swath_num = [load_product(fname).bursts[0].swathNumber for fname in xml_files]
        meta['swathNumber'] = ''.join(str(i) for i in sorted(swath_num))

    # 2. calculate ASF frame number for Sentinel-1
    meta['firstFrameNumber'] = int(0.2 * (burst.burstStartUTC   - obj.ascendingNodeTime).total_seconds())
    meta['lastFrameNumber']  = int(0.2 * (burstEnd.burstStopUTC - obj.ascendingNodeTime).total_seconds())
    return meta, burst

def load_product(xml_name):
    """Load the product using Product Manager."""
    import isce
    from iscesys.Component.ProductManager import ProductManager
    pm = ProductManager()
    pm.configure()
    return pm.loadProduct(xml_name)
