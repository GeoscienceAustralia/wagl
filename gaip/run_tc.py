"""
Runs the terrain correction. This code runs the terrain correction algorithm.
"""

import logging, os, numpy, gc
from osgeo import gdal
from ULA3 import DataManager, DataGrid
from ULA3.tc import clip_dsm, filter_dsm, run_slope, run_castshadow, run_brdfterrain # , run_brdfterrain_LS8
from ULA3.utils import Buffers, dump_array, load_bin_file, as_array
from ULA3.dataset import SceneDataset
from ULA3.image_processor import ProcessorConfig

from ULA3.geodesic import calculate_angles as ca
from ULA3.tests import unittesting_tools as ut

from rasterio.warp import RESAMPLING
from gaip import write_img
from gaip import find_file
from gaip import read_img

logger = logging.getLogger('root.' + __name__)





def write_tif_file(l1t_input_dataset, band_data, filename, file_type):
    dump_array(
        array=band_data,
        output_path=filename,
        template_dataset=l1t_input_dataset,
        geoc=l1t_input_dataset.geotransform,
        proj=l1t_input_dataset.spatial_ref.ExportToWkt(),
        no_data_value=-999,
        convert_to_byte=False,
        file_format=file_type)





def write_new_brdf_file(file_name, *args):
    output = open(file_name, 'w')
    output.write("%f\n%f %f %f\n%f %f %f %f\n%f\n" % args)
    output.close()





# args to the command line
# - image-buffer width (singular at moment, but may need to be multiple later). will be accessible as CONFIG.get_item('tc_dsm_buffer_width', int)
dsm_buffer_width = 250
# - width of submatrix used in run_castshadow. will be accessible as CONFIG.get_item('shadow_dsm_submatrix_width', int)
shadow_sub_matrix_height = 500
# - height of submatrix used in run_castshadow. will be accessible as CONFIG.get_item('shadow_dsm_submatrix_height', int)
shadow_sub_matrix_width = 500
# one of the parameters required for terrain correction - need to figure out how this is specified.
rori = 0.52
# these were copied from the files files in /g/data/v10/ULA3-TESTDATA/brdf_modis_band%i.txt,
# they were contained in the last line of those files.
ave_reflectance_values_LS5_7 = {10:0.0365, 20:0.0667, 30:0.0880, 40:0.2231, 50:0.2512, 70:0.1648}
ave_reflectance_values_LS8   = {1:0.0365, 2:0.0365, 3:0.0667, 4:0.0880, 5:0.2231, 6:0.2512, 7:0.1648}





def process(subprocess_list=[], resume=False):
    #logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    #CONFIG = ProcessorConfig()
    #DATA = DataManager()

    dump_path = False
    output_format = 'ENVI'
    output_extension = '.img'
    #l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    #is_utm =  not l1t_input_dataset.IsGeographic()
    output_path = DATA.get_item('nbar_temp_output.dat', str)
    work_path = CONFIG.work_path
    nbar_dataset_id = DATA.get_item('nbar_dataset_id.dat', str)

    #dsm data
    national_dsm_path = os.path.join(CONFIG.DIR_DEM_TC, 'dsm1sv1_0_Clean.img')
    pixel_buf = Buffers(dsm_buffer_width)
    output_dsm_path = os.path.join(work_path, 'region_dsm_image' + output_extension)

    # TODO Rework the clip and filter to not use the SceneDataset    
    dsm_data = filter_dsm(clip_dsm(l1t_input_dataset, national_dsm_path, output_dsm_path, pixel_buf, output_format))

    pref = ''


    ###################################################################

    # Use the 1st acquisition to setup the geobox
    geobox = gridded_geo_box(aqcuistions[0])

    # Retrive the spheroid parameters
    # (used in calculating pixel size in metres per lat/lon)
    spheroid = ca.setup_spheroid(geobox.crs.ExportToWkt())

    # Are we in projected or geographic space
    is_utm = not geobox.crs.IsGeographic()

    # National DSM data path
    national_dsm_path = os.path.join(CONFIG.DIR_DEM_TC, 'dsm1sv1_0_Clean.img')
    pixel_buf = Buffers(dsm_buffer_width)

    # Define Top, Bottom, Left, Right pixel buffers
    pixel_buf = Buffers(dsm_buffer_width)

    # Get the dimensions and geobox of the new image
    dem_cols = geobox.getShapeXY()[0] + pixel_buf.left + pixel_buf.right
    dem_rows = geobox.getShapeXY()[1] + pixel_buf.top + pixel_buf.bottom
    dem_shape = (dem_rows, dem_cols)
    dem_origin = geobox.convert_coordinates((0 - pixel_buf.left,
        0 - pixel_buf.top))
    dem_geobox = GriddedGeoBox(dem_shape, origin=dem_origin,
        pixelsize=geobox.pixelsize, crs=geobox.crs.ExportToWkt())

    # Retrive the DSM data
    dsm_data = reprojectFile2Array(national_dsm_path, dst_geobox=dem_geobox,
        resampling=RESAMPLING.bilinear)

    # Output the reprojected result
    write_img(dsm_data, output_dsm_path, geobox=dem_geobox)

    # Smooth the DSM
    dsm_data = filter_dsm(dsm_data)

    # Output the smoothed DSM
    output_smDSM_path = os.path.join(work_path, 'region_dsm_image_smoothed' + output_extension)
    write_img(dsm_data, output_smDSM_path, geobox=dem_geobox)

    #TODO Need a new debug flag for the luigi workflow
    if CONFIG.debug:
        # the location to dump the data.
        # TODO change important_intermediates to TC_intermediates??
        dump_path = os.path.join(work_path, 'important_intermediates')
        try:
            os.mkdir(dump_path)
        except OSError:
            if not os.path.exists(dump_path):
                dump_path = False
                logger.error("error creating directory for important intermediates.")

    if dump_path:
        # write the equivalent input file for Fuqin.
        l1_shape = l1t_input_dataset.bounds_getter(l1t_input_dataset)
        l1_shape.write_header_slope_file(os.path.join(dump_path, 'SLOPE_ANGLE_INPUTS'), pixel_buf)

    # solar angle data
    fname = find_file(work_path, 'SOL_Z.bin')
    solar_angle = read_img(fname)
    fname = find_file(work_path, 'SOL_AZ.bin')
    sazi_angle = read_img(fname)
    

    # satellite angle data
    fname = find_file(work_path, 'SAT_V.bin')
    view_angle = read_img(fname)
    fname = find_file(work_path, 'SAT_AZ.bin')
    azi_angle = read_img(fname)
    fname = find_file(work_path, 'REL_AZ.bin')
    rela_angle = read_img(fname)


    # TODO re-work this routine
    # calculate the slope and angle
    slope_results = run_slope(l1t_input_dataset, dsm_data, solar_angle,
                              view_angle, sazi_angle, azi_angle, pixel_buf,
                              is_utm, spheroid)

    # TODO find out what shadow_s & shadow_v are
    shadow_s = run_castshadow(acquisition[0], dsm_data, solar_angle,
        sazi_angle, pixel_buf, shadow_sub_matrix_height,
        shadow_sub_matrix_width, spheroid)

    shadow_v = run_castshadow(acquisition[0], dsm_data, view_angle,
        azi_angle, pixel_buf, shadow_sub_matrix_height,
        shadow_sub_matrix_width, spheroid)

    write_img(shadow_s, 'shadow_s.img', geobox=geobox)
    write_img(shadow_v, 'shadow_v.img', geobox=geobox)

    if dump_path:
        write_tif_file(l1t_input_dataset, shadow_s, os.path.join(dump_path, 'shadow_s.img'), file_type = 'ENVI')
        write_tif_file(l1t_input_dataset, shadow_v, os.path.join(dump_path, 'shadow_v.img'), file_type = 'ENVI')
        # TODO re-work this slope_results routine
        slope_results.dump_arrays(dump_path, l1t_input_dataset, "ENVI", ".img")

    # load the line starts and ends.
    region_nrow, region_ncol = l1t_input_dataset.shape

    boo = DATA.get_item("bilinear_ortho_outputs", item_type=dict)

    # this process runs close to the wind on memory... get rid of everything for now.
    DATA.clean()
    gc.collect()

    if (l1t_input_dataset.satellite.NAME == 'Landsat-8'):
        ave_reflectance_values = ave_reflectance_values_LS8
    else:
        ave_reflectance_values = ave_reflectance_values_LS5_7

    # We need a reverse lookup for the band file names, ie for LS_5/7: 10, 20, 30 etc
    # to correspond with the 1 based index order retrieved from scene_dataset_class.bands('REFLECTIVE')
    # Was unsure if there was already a lookup defined so this will suffice.
    # This will be used for output filenames as well as average reflectance lookups.
    band_fname_lookup = {}
    for key in l1t_input_dataset._band_number_map.keys():
        val = l1t_input_dataset._band_number_map[key]
        band_fname_lookup[val] = key # Basically turn keys into values and vice versa. Only works for 1 to 1 mapping dicts.

    #for band_number in (2, 3, 4): #(1, 2, 3, 4, 5, 7):
    for band_number in l1t_input_dataset.bands('REFLECTIVE'):
        # Get the band file name, eg 10, 20, 30 etc, used in L1T. LS8 uses 1, 2, 3 etc
        out_bn_name = band_fname_lookup[band_number]

        # not sure where these get created.
        param_file = open(os.path.join(CONFIG.work_path, 'brdf_modis_band%i.txt' % band_number), 'r')
        brdf0, brdf1, brdf2, bias, slope_ca, esun, dd = map(float, ' '.join(param_file.readlines()).split())
        param_file.close()

        if dump_path:
            write_new_brdf_file(
                os.path.join(dump_path, 'new_brdf_modis_band%i.txt' % band_number),
                rori, brdf0, brdf1, brdf2, bias, slope_ca, esun, dd, ave_reflectance_values[out_bn_name])

        # need to check that these are OK.
        band_data = l1t_input_dataset.band_read_as_array(band_number)

	ref_lm, ref_brdf, ref_terrain = run_brdfterrain(
	    rori,
	    brdf0, brdf1, brdf2,
	    bias, slope_ca, esun, dd,
	    ave_reflectance_values[out_bn_name],
	    band_data,
	    slope_results.mask_self,
	    shadow_s,
	    shadow_v,
	    solar_angle,
	    sazi_angle,
	    view_angle,
	    rela_angle,
	    slope_results.slope,
	    slope_results.aspect,
	    slope_results.incident,
	    slope_results.exiting,
	    slope_results.rela_slope,
	    load_bin_file(boo[(band_number, 'a')], region_nrow, region_ncol, dtype=numpy.float32),
	    load_bin_file(boo[(band_number, 'b')], region_nrow, region_ncol, dtype=numpy.float32),
	    load_bin_file(boo[(band_number, 's')], region_nrow, region_ncol, dtype=numpy.float32),
	    load_bin_file(boo[(band_number, 'fs')], region_nrow, region_ncol, dtype=numpy.float32),
	    load_bin_file(boo[(band_number, 'fv')], region_nrow, region_ncol, dtype=numpy.float32),
	    load_bin_file(boo[(band_number, 'ts')], region_nrow, region_ncol, dtype=numpy.float32),
	    load_bin_file(boo[(band_number, 'dir')], region_nrow, region_ncol, dtype=numpy.float32),
	    load_bin_file(boo[(band_number, 'dif')], region_nrow, region_ncol, dtype=numpy.float32))


        write_tif_file(l1t_input_dataset, ref_lm, os.path.join(work_path, pref + 'ref_lm_' + str(out_bn_name) + output_extension), file_type = output_format)
        write_tif_file(l1t_input_dataset, ref_brdf, os.path.join(work_path, pref + 'ref_brdf_' + str(out_bn_name) + output_extension), file_type = output_format)
        write_tif_file(l1t_input_dataset, ref_terrain, os.path.join(work_path, pref + 'ref_terrain_' + str(out_bn_name) + output_extension), file_type = output_format)

        outfname = os.path.join(output_path, 'scene01', '%s_B%d%s' % (nbar_dataset_id, out_bn_name, '.tif'))
        write_tif_file(l1t_input_dataset, ref_terrain, outfname, file_type="GTiff")

        ref_lm = ref_brdf = ref_terrain = None
        gc.collect()

