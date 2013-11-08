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
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    dump_path = False
    output_format = 'ENVI'
    output_extension = '.img'
    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    is_utm =  not l1t_input_dataset.IsGeographic()
    output_path = DATA.get_item('nbar_temp_output.dat', str)
    work_path = CONFIG.work_path
    nbar_dataset_id = DATA.get_item('nbar_dataset_id.dat', str)

    #dsm data
    national_dsm_path = os.path.join(CONFIG.DIR_DEM_TC, 'dsm1sv1_0_Clean.img')
    pixel_buf = Buffers(dsm_buffer_width)
    output_dsm_path = os.path.join(work_path, 'region_dsm_image' + output_extension)
    dsm_data = filter_dsm(clip_dsm(l1t_input_dataset, national_dsm_path, output_dsm_path, pixel_buf, output_format))
    pref = ''

    if CONFIG.debug:
        # the location to dump the data.
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

    # solar angle data (these are set in ULA3.image_processor.utils.calc_solar_grids).
    solar_angle = as_array(DATA.get_item('SOL_Z_DEG.bin', DataGrid).array, dtype=numpy.float32)
    sazi_angle = as_array(DATA.get_item('SOL_AZ_DEG.bin', DataGrid).array, dtype=numpy.float32)

    # satellite angle data (these are set in ULA3.image_processor.utils.calc_satellite_grids)
    view_angle = as_array(DATA.get_item('SAT_V_DEG.bin', DataGrid).array, dtype=numpy.float32)
    azi_angle = as_array(DATA.get_item('SAT_AZ_DEG.bin', DataGrid).array, dtype=numpy.float32)
    rela_angle = as_array(DATA.get_item('REL_AZ_DEG.bin', DataGrid).array, dtype=numpy.float32)

    if dump_path:
        write_tif_file(l1t_input_dataset, solar_angle, os.path.join(dump_path, 'solar_angle.img'), file_type = 'ENVI')
        write_tif_file(l1t_input_dataset, sazi_angle, os.path.join(dump_path, 'sazi_angle.img'), file_type = 'ENVI')
        write_tif_file(l1t_input_dataset, rela_angle, os.path.join(dump_path, 'rela_angle.img'), file_type = 'ENVI')
        write_tif_file(l1t_input_dataset, view_angle, os.path.join(dump_path, 'view_angle.img'), file_type = 'ENVI')
        write_tif_file(l1t_input_dataset, azi_angle, os.path.join(dump_path, 'azi_angle.img'), file_type = 'ENVI')

    def do_run_castshadow(zen_angle, azi_angle):
        return run_castshadow(
            l1t_input_dataset,
            dsm_data,
            zen_angle,
            azi_angle,
            pixel_buf,
            shadow_sub_matrix_height,
            shadow_sub_matrix_width,
            is_utm)

    # calculate the slope and angle
    slope_results = run_slope(l1t_input_dataset, dsm_data, solar_angle, view_angle, sazi_angle, azi_angle, pixel_buf, is_utm)
    shadow_s = do_run_castshadow(solar_angle, sazi_angle)
    shadow_v = do_run_castshadow(view_angle, azi_angle)

    if dump_path:
        write_tif_file(l1t_input_dataset, shadow_s, os.path.join(dump_path, 'shadow_s.img'), file_type = 'ENVI')
        write_tif_file(l1t_input_dataset, shadow_v, os.path.join(dump_path, 'shadow_v.img'), file_type = 'ENVI')
        slope_results.dump_arrays(dump_path, l1t_input_dataset, "ENVI", ".img")

    # load the line starts and ends.
    region_nrow, region_ncol = l1t_input_dataset.shape
    # the istart and iend params are no longer passed to brdf_terrain as the function doesn't use them
    istart = numpy.zeros(region_nrow, dtype=numpy.int32)
    iend = numpy.zeros(region_nrow, dtype=numpy.int32)
    start_end_file = os.path.join(CONFIG.work_path, 'STARTEND')
    start_end = open(start_end_file)
    line_no = 0
    for line in start_end:
        bits = line.split()
        istart[line_no] = int(bits[1])
        iend[line_no] = int(bits[3])
        line_no = line_no + 1
    start_end.close()
    assert region_nrow == line_no, "read wrong number of items from start_end file."

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

        """
        # At this point in time we have different instances of brdfterrain
        # One for LS5/7 and one for LS8. They're different datatypes and we're
        # preserving the science code which is FORTRAN. THIS NEEDS TO CHANGE!!!
        if (l1t_input_dataset.satellite.NAME == 'Landsat-8'):
            ref_lm, ref_brdf, ref_terrain = run_brdfterrain_LS8(
                rori,
                brdf0, brdf1, brdf2,
                bias, slope_ca, esun, dd,
                ave_reflectance_values[band_number],
                istart, iend,
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
        else:
        """
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


        print
        print 'band_number: ', band_number
        print 'out_bn_name: ', out_bn_name
        print

        #write_tif_file(l1t_input_dataset, ref_lm, os.path.join(work_path, pref + 'ref_lm_' + str(band_number) + output_extension), file_type = output_format)
        write_tif_file(l1t_input_dataset, ref_lm, os.path.join(work_path, pref + 'ref_lm_' + str(out_bn_name) + output_extension), file_type = output_format)
        #write_tif_file(l1t_input_dataset, ref_brdf, os.path.join(work_path, pref + 'ref_brdf_' + str(band_number) + output_extension), file_type = output_format)
        write_tif_file(l1t_input_dataset, ref_brdf, os.path.join(work_path, pref + 'ref_brdf_' + str(out_bn_name) + output_extension), file_type = output_format)
        #write_tif_file(l1t_input_dataset, ref_terrain, os.path.join(work_path, pref + 'ref_terrain_' + str(band_number) + output_extension), file_type = output_format)
        write_tif_file(l1t_input_dataset, ref_terrain, os.path.join(work_path, pref + 'ref_terrain_' + str(out_bn_name) + output_extension), file_type = output_format)

        #outfname = os.path.join(output_path, 'scene01', '%s_B%d%s' % (nbar_dataset_id, band_number, '.tif'))
        outfname = os.path.join(output_path, 'scene01', '%s_B%d%s' % (nbar_dataset_id, out_bn_name, '.tif'))
        write_tif_file(l1t_input_dataset, ref_terrain, outfname, file_type="GTiff")

        ref_lm = ref_brdf = ref_terrain = None
        gc.collect()

