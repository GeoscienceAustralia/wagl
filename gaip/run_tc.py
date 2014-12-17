"""
Runs the terrain correction. This code runs the terrain correction algorithm.
"""

import gc
import os
from os.path import join as pjoin
import pickle
import numpy

from rasterio.warp import RESAMPLING

from gaip import write_img
from gaip import find_file
from gaip import read_img
from gaip import load_2D_bin_file
from gaip import calculate_angles as ca
from gaip import run_slope
from gaip import run_castshadow
from gaip import run_brdfterrain
from gaip import GriddedGeoBox
from gaip import reprojectFile2Array
from gaip import constants
from gaip import Buffers
from gaip import filter_dsm
from gaip import write_header_slope_file
from gaip import write_new_brdf_file


def run_tc(acquisitions, dsm_buffer_width, shadow_sub_matrix_height,
    shadow_sub_matrix_width, rori, national_dsm, work_path):
    """
    The terrain correction workflow.

    :param acquisitions:
        A list of acquisition class objects that will be run through
        the terrain correction workflow.

    :param dsm_buffer_width:
        The buffer in pixels around the acquisition dimensions. Used
        for the extracting a subset from the national digital surface
        model.

    :param shadow_sub_matrix_height:
        The height (rows) of the window/submatrix used in the cast
        shadow algorithm.

    :param shadow_sub_matrix_width:
        The width (rows) of the window/submatrix used in the cast
        shadow algorithm.

    :param rori:
        Threshold for terrain correction.

    :param national_dsm:
        A string containing the full file system path to the national
        digital surface model image file.

    :param work_path:
        A full file system path to the working directory.
        Intermediate files will be saved to
        `work_path/tc_intermediates`.

    :return:
        None.
        The terrain correction algorithm will output 3 files for every
        band.
        ref_lm_{band_number}.img -> Lambertian reflectance.
        ref_brdf_{band_number}.img -> BRDF corrected reflectance.
        ref_terrain_{band_number}.img -> Terrain corrected
            reflectance.
    """
    # Terrain correction working path
    tc_work_path = pjoin(work_path, 'tc_intermediates')
    try:
        os.mkdir(tc_work_path)
    except OSError:
        if not os.path.exists(tc_work_path):
            msg = ("Error creating directory for terrain correction "
                   "intermediates.")
            print msg
            raise

    # Use the 1st acquisition to setup the geobox
    geobox = acquisitions[0].gridded_geo_box()

    # Retrive the spheroid parameters
    # (used in calculating pixel size in metres per lat/lon)
    spheroid = ca.setup_spheroid(geobox.crs.ExportToWkt())

    # Are we in projected or geographic space
    is_utm = not geobox.crs.IsGeographic()

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
    dsm_data = reprojectFile2Array(national_dsm, dst_geobox=dem_geobox,
        resampling=RESAMPLING.bilinear)

    # Output the reprojected result
    fname_DSM_subset = pjoin(tc_work_path, 'region_dsm_image.img')
    write_img(dsm_data, fname_DSM_subset, geobox=dem_geobox)

    # Smooth the DSM
    dsm_data = filter_dsm(dsm_data)

    # Output the smoothed DSM
    fname_smDSM = pjoin(tc_work_path, 'region_dsm_image_smoothed.img')
    write_img(dsm_data, fname_smDSM, geobox=dem_geobox)


    # write the equivalent input file for Fuqin.
    write_header_slope_file(pjoin(tc_work_path, 'SLOPE_ANGLE_INPUTS'),
        pixel_buf, geobox)

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


    # calculate the slope and angle
    slope_results = run_slope(acquisitions[0], dsm_data, solar_angle,
        view_angle, sazi_angle, azi_angle, pixel_buf, is_utm, spheroid)

    # Output slope results
    slope_results.write_arrays(tc_work_path, geobox, "ENVI", ".img")

    # TODO find out what shadow_s & shadow_v are
    shadow_s = run_castshadow(acquisitions[0], dsm_data, solar_angle,
        sazi_angle, pixel_buf, shadow_sub_matrix_height,
        shadow_sub_matrix_width, spheroid)

    shadow_v = run_castshadow(acquisitions[0], dsm_data, view_angle,
        azi_angle, pixel_buf, shadow_sub_matrix_height,
        shadow_sub_matrix_width, spheroid)

    # Output the two shadow masks to disk
    fname_shadow_s = pjoin(tc_work_path, 'shadow_s.img')
    fname_shadow_v = pjoin(tc_work_path, 'shadow_v.img')
    write_img(shadow_s, fname_shadow_s, geobox=geobox, nodata=-999)
    write_img(shadow_v, fname_shadow_v, geobox=geobox, nodata=-999)

    # load the line starts and ends.
    rows, cols = geobox.shape

    # load the pickled dict containing the bilinear outputs
    boo_fname = pjoin(work_path, '"bilinear_ortho_outputs')
    with open(boo_fname) as boo_file:
        boo = pickle.load(boo_file)

    # this process runs close to the wind on memory... get rid of everything for now.
    gc.collect()

    # Retrieve the satellite and sensor for the acquisition
    satellite = acquisitions[0].spacecraft_id
    sensor = acquisitions[0].sensor_id

    # Get the required nbar bands list for processing
    nbar_constants = constants.NBARConstants(satellite, sensor)
    bands_list = nbar_constants.getNBARlut()
    ave_reflectance_values = nbar_constants.getAvgReflut()

    for acq in acquisitions:
        band_number = acq.band_num
        if band_number not in bands_list:
            # skip
            continue

        # Read the BRDF modis file for a given band
        brdf_modis_file = 'brdf_modis_band{0}.txt'.format(band_number)
        brdf_modis_file = pjoin(work_path, brdf_modis_file)

        # Read the BRDF modis file for a given band
        brdf_modis_file = 'brdf_modis_band{0}.txt'.format(band_number)
        brdf_modis_file = os.path.join(work_path, brdf_modis_file)
        with open(brdf_modis_file, 'r') as param_file:
            brdf0, brdf1, brdf2, bias, slope_ca, esun, dd = map(float,
                ' '.join(param_file.readlines()).split())

        write_new_brdf_file(pjoin(tc_work_path,
            'new_brdf_modis_band{band_num}.txt'.format(band_num=band_number)),
            rori, brdf0, brdf1, brdf2, bias, slope_ca, esun, dd,
            ave_reflectance_values[band_number])

        # Read the data
        band_data = acq.data()

        # Run terrain correction
	ref_lm, ref_brdf, ref_terrain = run_brdfterrain(
	    rori,
	    brdf0, brdf1, brdf2,
	    bias, slope_ca, esun, dd,
	    ave_reflectance_values[band_number],
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
	    load_2D_bin_file(boo[(band_number, 'a')], rows, cols,
                dtype=numpy.float32),
	    load_2D_bin_file(boo[(band_number, 'b')], rows, cols,
                dtype=numpy.float32),
	    load_2D_bin_file(boo[(band_number, 's')], rows, cols,
                dtype=numpy.float32),
	    load_2D_bin_file(boo[(band_number, 'fs')], rows, cols,
                dtype=numpy.float32),
	    load_2D_bin_file(boo[(band_number, 'fv')], rows, cols,
                dtype=numpy.float32),
	    load_2D_bin_file(boo[(band_number, 'ts')], rows, cols,
                dtype=numpy.float32),
	    load_2D_bin_file(boo[(band_number, 'dir')], rows, cols,
                dtype=numpy.float32),
	    load_2D_bin_file(boo[(band_number, 'dif')], rows, cols,
                dtype=numpy.float32))


        # Output filenames for lambertian, brdf and terrain corrected reflectance
        lmbrt_fname = pjoin(work_path, 'ref_lm_{}.img'.format(band_number))
        brdf_fname = pjoin(work_path, 'ref_brdf_{}.img'.format(band_number))
        tc_fname = pjoin(work_path, 'ref_terrain_{}.img'.format(band_number))

        # Output the files.
        write_img(ref_lm, lmbrt_fname, geobox=geobox, nodata=-999)
        write_img(ref_brdf, brdf_fname, geobox=geobox, nodata=-999)
        write_img(ref_terrain, tc_fname, geobox=geobox, nodata=-999)
        
        ref_lm = ref_brdf = ref_terrain = None
        gc.collect()
