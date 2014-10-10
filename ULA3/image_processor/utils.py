"""
Provides some utilities that are used in the image processor. Most of the functions contained herein are not likely
to be of use outside of the image processor at this stage as they still rely heavily on
:py:class:`ULA3.image_processor.ProcessorConfig` and :py:class:`ULA3.DataManager`.

:todo:
    Get Rid of dependencies on :py:class:`ULA3.image_processor.ProcessorConfig` and :py:class:`ULA3.DataManager`.
"""
import logging, errno, os, re, numpy, ephem
from copy import copy
from datetime import datetime

from ULA3 import DataGrid
from ULA3.utils import create_output_image as create_image
from ULA3.utils import log_multiline, execute
from ULA3.meta import print_call
from ULA3.dataset import SceneDataset
from ULA3.metadata import Metadata, XMLMetadata, ReportMetadata
from ULA3.geodesic import Satellite, eval_centre_data, eval_sat_grids
from ULA3.solar import eval_sol_grids

# Imports for newer angle calculations method
from ULA3.geodesic import calculate_angles as ca
from ULA3.tests import unittesting_tools as ut

logger = logging.getLogger('root.' + __name__)

@print_call(logger.info)
def create_output_image(DATA, CONFIG, rasterXSize=None, work_dir_prefix=None, image_postfix=None):
    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve input scene dataset'
    logger.debug('SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    nbar_dataset_id = DATA.get_item('nbar_dataset_id.dat', str)
    assert nbar_dataset_id, 'Unable to retrieve nbar_dataset_id string'
    logger.debug('string for nbar_dataset_id retrieved')

    nbar_temp_output = DATA.get_item('nbar_temp_output.dat', str)
    assert nbar_temp_output, 'Unable to retrieve nbar_temp_output string'
    logger.debug('string for nbar_temp_output retrieved')

    nbar_output_path = DATA.get_item('nbar_output_path.dat', str)
    assert nbar_output_path, 'Unable to retrieve nbar_output_path string'
    logger.debug('string for nbar_output_path retrieved')

    pfiles = [os.path.join(nbar_temp_output, 'scene01',
        nbar_dataset_id + '_B' + str(band_file_number) + '.tif') for
        band_file_number in l1t_input_dataset.satellite.rgb_bands]

    if not rasterXSize:
        rasterXSize = l1t_input_dataset.RasterXSize
        image_postfix = image_postfix or '_FR'
        work_dir_prefix = work_dir_prefix or 'full_res'
    else:
        assert work_dir_prefix, "argument work_dir_prefix must not be None"
        image_postfix = image_postfix or ''

    browse_image_path = os.path.join(nbar_temp_output, nbar_dataset_id + image_postfix + '.jpg')
    thumbnail_work_dir = os.path.join(CONFIG.work_path, work_dir_prefix)

    try:
        os.makedirs(thumbnail_work_dir)
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise

    res = {
        'image_dim':create_image(pfiles, browse_image_path, thumbnail_work_dir, rasterXSize),
        'l1t_input_dataset':l1t_input_dataset,
        'browse_image_path':browse_image_path}

    return res




import shutil
@print_call(logger.info)
def move_outputs(from_dir, to_dir):
    if from_dir == to_dir:
        logger.info('Source and destinations directories the same (%s). No files moved', from_dir)
        return
    shutil.rmtree(to_dir, ignore_errors=True)
    shutil.move(from_dir, to_dir)






@print_call(logger.info)
def calc_lat_long_grids(DATA, CONFIG):
    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug('SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    lon_grid = DataGrid(shape=l1t_input_dataset.shape, eval_func=l1t_input_dataset.geo_lon, depth=CONFIG.blrb_depth)
    lat_grid = DataGrid(shape=l1t_input_dataset.shape, eval_func=l1t_input_dataset.geo_lat, depth=CONFIG.blrb_depth)

    # Register grids required for downstream processes (i.e. grid computation)
    # These do NOT need to be saved to disk
    DATA.set_item('LON_RAD.bin', lon_grid)
    DATA.set_item('LAT_RAD.bin', lat_grid)

    # These grids are required for modtran input
    lon_grid.save_binary(filename=os.path.join(CONFIG.work_path, 'LON.bin'), convert_angles=True, dtype=CONFIG.FORTRAN_GRID_DATATYPE)
    lat_grid.save_binary(filename=os.path.join(CONFIG.work_path, 'LAT.bin'), convert_angles=True, dtype=CONFIG.FORTRAN_GRID_DATATYPE)

    if CONFIG.debug:
        lon_grid.save_image(filename=os.path.join(CONFIG.work_path, 'LON.tif'), convert_angles=True)
        lat_grid.save_image(filename=os.path.join(CONFIG.work_path, 'LAT.tif'), convert_angles=True)






@print_call(logger.info)
def calc_satellite_grids(DATA, CONFIG):
    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    satellite = l1t_input_dataset.satellite
    assert satellite, 'Unable to find Satellite object'

    # Load TLE file for the scene date and time.
    # N.B: EarthSatellite objects cannot be pickled, so we need to reload this from
    # a Satellite object if it is needed in another module
    sat_ephem = satellite.load_tle(l1t_input_dataset.scene_centre_datetime, CONFIG.EPHEM_DATA_ROOT, CONFIG.TLE_SEARCH_RANGE)
    assert sat_ephem, "Unable to load %s TLE for %s" % (satellite.NAME, l1t_input_dataset.scene_centre_datetime)

    if CONFIG.debug:
        # Test sat_ephem object
        sat_ephem.compute(l1t_input_dataset.scene_centre_datetime)
        observer = ephem.Observer()
        observer.long = sat_ephem.sublong
        observer.lat = sat_ephem.sublat
        observer.elevation = 0.0
        sat_ephem.compute(observer)
        logger.debug('sat_ephem.alt = %s', sat_ephem.alt)
        logger.debug('sat_ephem.az = %s', sat_ephem.az)
        logger.debug('sat_ephem.elevation = %s', sat_ephem.elevation)

    lat_grid = DATA.get_item('LAT_RAD.bin', DataGrid)
    assert lat_grid, "Unable to load DataGrid for LAT_RAD.bin"
    lon_grid = DATA.get_item('LON_RAD.bin', DataGrid)
    assert lon_grid, "Unable to load DataGrid for LON_RAD.bin"

    centre_data = eval_centre_data(l1t_input_dataset, sat_ephem)
    logger.debug('Centre data dict created')
    log_multiline(logger.debug, centre_data, 'CENTRE INFO', '\t')
    DATA.set_item('centre_data.dat', centre_data) # Register for use downstream

    (time, sat_v, sat_az) = eval_sat_grids(lon_grid.array, lat_grid.array, centre_data, l1t_input_dataset, sat_ephem, CONFIG.view_angle_max, CONFIG.work_path, CONFIG.debug)

    # Register in-memory grids for use in downstream processes and save them to disk
    # N.B: Angles stored as radians
    sat_az_grid = DataGrid(array=sat_az)
    DATA.set_item('SAT_AZ_RAD.bin', sat_az_grid)
    DATA.set_item('SAT_AZ_DEG.bin', DataGrid(array=numpy.degrees(sat_az)))

    # Need to save sat_az to disk
    sat_az_grid.save_binary(os.path.join(CONFIG.work_path, 'SAT_AZ.bin'), convert_angles=True, dtype=CONFIG.FORTRAN_GRID_DATATYPE)

    # Needed for solar grid calculation downstream - no angular conversion required
    time_grid = DataGrid(array=time)
    DATA.set_item('TIME.bin', time_grid)

    # Need to save sat_v to disk but we don't need an in-memory reference kept
    sat_v_grid = DataGrid(array=sat_v)
    sat_v_grid.save_binary(os.path.join(CONFIG.work_path, 'SAT_V.bin'), convert_angles=True, dtype=CONFIG.FORTRAN_GRID_DATATYPE)
    DATA.set_item('SAT_V_DEG.bin', DataGrid(array=numpy.degrees(sat_v)))

    if CONFIG.debug:
        # Save debug files to disk with angles in degrees
        sat_v_grid.save_image(filename=os.path.join(CONFIG.work_path, 'SAT_V.tif'), convert_angles=True)
        sat_az_grid.save_image(filename=os.path.join(CONFIG.work_path, 'SAT_AZ.tif'), convert_angles=True)
        time_grid.save_image(filename=os.path.join(CONFIG.work_path, 'TIME.tif'), convert_angles=False)





@print_call(logger.info)
def calc_solar_grids(DATA, CONFIG):
    """
    Generate all solar grids.
    """
    #l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    #assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    #logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    lon_grid = DATA.get_item('LON_RAD.bin', DataGrid)
    assert(lon_grid), 'Unable to retrieve lon_grid'

    lat_grid = DATA.get_item('LAT_RAD.bin', DataGrid)
    assert(lat_grid), 'Unable to retrieve lat_grid'

    centre_data = DATA.get_item('centre_data.dat', dict)
    assert(centre_data), 'Unable to retrieve centre_data'

    time_grid = DATA.get_item('TIME.bin', DataGrid)
    assert(time_grid), 'Unable to retrieve time_grid'

    sat_az_grid = DATA.get_item('SAT_AZ_RAD.bin', DataGrid)
    assert(sat_az_grid), 'Unable to retrieve sat_az_grid'

    (sol_z, sol_az) = eval_sol_grids(lon_grid.array, lat_grid.array, time_grid.array, centre_data['day_of_year'])

    sol_z_grid = DataGrid(array=sol_z)
    rel_az_grid = DataGrid(array=(sat_az_grid.array - sol_az))

    DATA.set_item('SOL_AZ_DEG.bin', DataGrid(array=numpy.degrees(sol_az)))
    DATA.set_item('SOL_Z_DEG.bin', DataGrid(array=numpy.degrees(sol_z)))
    DATA.set_item('REL_AZ_DEG.bin', DataGrid(array=numpy.degrees(rel_az_grid.array)))

    # Save solar & relative azimuth grids as binary files containing angles in degrees
    # Required by brdf_sim_bin
    sol_z_grid.save_binary(os.path.join(CONFIG.work_path, 'SOL_Z.bin'), convert_angles=True, dtype=CONFIG.FORTRAN_GRID_DATATYPE)
    rel_az_grid.save_binary(os.path.join(CONFIG.work_path, 'REL_AZ.bin'), convert_angles=True, dtype=CONFIG.FORTRAN_GRID_DATATYPE)

    if CONFIG.debug:
        # Save debug files to disk with angles in degrees
        sol_z_grid.save_image(os.path.join(CONFIG.work_path, 'SOL_Z.tif'), convert_angles=True)
        rel_az_grid.save_image(os.path.join(CONFIG.work_path, 'REL_AZ.tif'), convert_angles=True)
        DataGrid(array=sol_az).save_image(os.path.join(CONFIG.work_path, 'SOL_AZ.tif'), convert_angles=True)

    # Delete grids which are not required for downstream computation.
    DATA.free_item('LON_RAD.bin', DataGrid)
    DATA.free_item('LAT_RAD.bin', DataGrid)
    DATA.free_item('TIME.bin', DataGrid)
    #DATA.free_item('SAT_AZ_RAD.bin', DataGrid) // this is now used in terrain correction.
    DATA.free_item('centre_data.dat', dict)
    DATA.free_item('satellite.dat', Satellite)


@print_call(logger.info)
def create_centreline_file(y, x, n, cols, view_max, outdir, outfname='CENTRELINE'):
    """
    :param y:
        A 1D NumPy array of type int with the same shape as x & n.
        Details the row number starting at 1.

    :param x:
        A 1D NumPy array of type int with the same shape as y & n.
        Details the column number starting at 0.

    :param n:
        A 1D NumPy array of type int with the same shaoe as y & x.
        Details whether or not the track point coordinate is
        averaged.

    :param cols:
        An integer indicating the number of columns in the original
        array.

    :param view_max:
        An float indicating the maximum view angle of the satellite
        used for determining the centreline.

    :param outdir:
        A string containing the output directory.

    :param outfname:
        A string containing the output filename.
        Default is CENTRELINE.
    """

    rows = y.shape[0]

    # I'm not sure of the formatting used by FORTRAN, nothing was specified
    # but the outputs had spaces.
    # It might be more ideal to created it as a csv???
    fname = os.path.join(outdir, outfname)
    outf  = open(fname, 'w')
    
    outf.write('%f\n' %view_max)
    outf.write('%i      %i\n' %(rows, cols))

    for r in range(rows):
        outf.write('%i      %i       %f     %f     %f\n'%(y[r], x[r], n[r], n[r], n[r]))    

    outf.close()


@print_call(logger.info)
def calc_sat_sol_angle_grids(DATA, CONFIG):
    """
    Generate the satellite and solar grids.
    """

    # Track our working directory so we can save files
    work_dir = CONFIG.work_path

    # Find and open the longitude and lattitude files
    # Avoiding DataManger here. find_file will be used sparingly until a proper
    # workflow is written.
    # NOTE: find_file() will call sys.exit(2) if the file isn't found
    lon_fname = ut.find_file(work_dir, 'LON.tif')
    lat_fname = ut.find_file(work_dir, 'LAT.tif')

    lon_arr = ut.read_img(lon_fname)
    lat_arr = ut.read_img(lat_fname)

    # Get the array dimensions
    dims = lon_arr.shape
    cols = dims[1]
    rows = dims[0]

    # Get the L1T data. (Code borrowed from above)
    L1T_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert L1T_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', L1T_dataset.pathname)

    # Get the angles, time, & satellite track coordinates
    (satellite_zenith, satellite_azimuth, solar_zenith,
     solar_azimuth, relative_azimuth, time,
     Y_cent, X_cent, N_cent) = ca.calculate_angles(L1T_dataset, lon_arr,
                                                   lat_arr, npoints=12)

    # Image projection, geotransform
    prj  = L1T_dataset.GetProjection()
    geoT = L1T_dataset.GetGeoTransform()

    # Define the output file names
    sat_view_zenith_fname  = os.path.join(work_dir, 'SAT_V.bin')
    sat_azimuth_fname      = os.path.join(work_dir, 'SAT_AZ.bin')
    solar_zenith_fname     = os.path.join(work_dir, 'SOL_Z.bin')
    solar_azimuth_fname    = os.path.join(work_dir, 'SOL_AZ.bin')
    relative_azimuth_fname = os.path.join(work_dir, 'REL_AZ.bin')
    time_fname             = os.path.join(work_dir, 'TIME.bin')

    # Write the image files to disk
    logger.debug("Writing satelite view zenith angle: %s", sat_view_zenith_fname)
    ut.write_img(satellite_zenith, sat_view_zenith_fname, projection=prj,
                 geotransform=geoT)

    logger.debug("Writing satellite azimuth angle: %s", sat_azimuth_fname)
    ut.write_img(satellite_azimuth, sat_azimuth_fname, projection=prj,
                 geotransform=geoT)

    logger.debug("Writing solar zenith angle: %s", solar_zenith_fname)
    ut.write_img(solar_zenith, solar_zenith_fname, projection=prj,
                 geotransform=geoT)

    logger.debug("Writing solar azimith angle: %s", solar_azimuth_fname)
    ut.write_img(solar_azimuth, solar_azimuth_fname, projection=prj,
                 geotransform=geoT)

    logger.debug("Writing relative azimuth angle: %s", relative_azimuth_fname)
    ut.write_img(relative_azimuth, relative_azimuth_fname, projection=prj,
                 geotransform=geoT)

    logger.debug("Writing time array: %s", time_fname)
    ut.write_img(time, time_fname, projection=prj, geotransform=geoT)


    # Close files
    del satellite_zenith, satellite_azimuth, solar_zenith, solar_azimuth
    del relative_azimuth, time

    # Write out the CENTRELINE file
    create_centreline_file(Y_cent, X_cent, N_cent, cols, view_max=9.0,
                           outdir=work_dir)


