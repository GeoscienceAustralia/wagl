#!/usr/bin/env python

import gc

from gaip import constants
from gaip import load_2D_bin_file
from gaip import read_img
from gaip import terrain_correction
from gaip import write_img
from gaip import write_new_brdf_file


def run_tc(acquisitions, bilinear_ortho_filenames, rori, self_shadow_fname,
        cast_shadow_sun_fname, cast_shadow_satellite_fname,
        solar_zenith_fname, solar_azimuth_fname, satellite_view_fname,
        relative_angle_fname, slope_fname, aspect_fname, incident_angle_fname,
        exiting_angle_fname, relative_slope_fname, reflectance_filenames):
    """
    The terrain correction workflow.

    :param acquisitions:
        A list of acquisition class objects that will be run through
        the terrain correction workflow.

    :param bilinear_ortho_filenames:
        A dictionary with keys specified via a tuple of
        (band_number, factor) and the value corresponding to a full
        file pathname to the bilinearly interpolated flaot32 array.
        Valid factor strings are:
            * fv: MODTRAN output (fv).
            * fs: MODTRAN output (fs).
            * b: MODTRAN output (b).
            * s: MODTRAN output (s).
            * a: MODTRAN output (a).
            * dir: MODTRAN output (direct irradiance).
            * dif: MODTRAN output (diffuse irradiance).
            * ts: MODTRAN output (ts).

    :param rori:
        Threshold for terrain correction.

    :param self_shadow_fname:
        A string containing the full file path name to the self
        shadow mask image.

    :param cast_shadow_sun_fname:
        A string containing the full file path name to the sun cast
        shadow mask image.

    :param cast_shadow_satellite_fname:
        A string containing the full file path name to the satellite
        cast shadow mask image.

    :param solar_zenith_fname:
        A string containing the full file path name to the solar
        zenith angle image.

    :param solar_azimuth_fname:
        A string containing the full file path name to the solar
        azimuth angle image.

    :param satellite_view_fname:
        A string containing the full file path name to the satellite
        view angle image.

    :param relative_angle_fname:
        A string containing the full file path name to the relative
        angle image.

    :param slope_fname:
        A string containing the full file path name to the slope
        image.

    :param aspect_fname:
        A string containing the full file path name to the aspect
        image.

    :param incident_angle_fname:
        A string containing the full file path name to the incident
        angle image.

    :param exiting_angle_fname:
        A string containing the full file path name to the exiting
        angle image.

    :param relative_slope_fname:
        A string containing the full file path name to the relative
        slope image.

    :param reflectance_filenames:
        A dictionary with keys specified via a tuple of
        (band, reflectance_level) and the value corresponding to a
        full file pathname.
        Valid reflectance level string are:
            ref_lm -> Lambertian reflectance
            ref_brdf -> BRDF corrected reflectance
            ref_terrain -> Terrain corrected reflectance

    :return:
        None.
        The terrain correction algorithm will output 3 files for every
        band in the following format.
        reflectance_lambertian_{band_number}.bin -> Lambertian
            reflectance.
        reflectance_brdf_{band_number}.bin -> BRDF corrected
            reflectance.
        reflectance_terrain_{band_number}.bin -> Terrain corrected
            reflectance.

    :notes:
        Arrays will be converted to the required datatype and
        transposed. The trnasposing should prevent array copies
        being made by F2Py. The results are transposed back before
        being written to disk.
        All arrays should have the same dimensions.
        Required datatypes are as follows:
            * acquisitions: `numpy.int16`
            * self_shadow: `numpy.int16`
            * cast_shadow_sun: `numpy.int16`
            * cast_shadow_satellite: `numpy.int16`
            * solar_zenith: `numpy.float32`
            * solar_azimuth: `numpy.float32`
            * satellite_view: `numpy.float32`
            * relative_angle: `numpy.float32`
            * slope: `numpy.float32`
            * aspect: `numpy.float32`
            * incident_angle: `numpy.float32`
            * exiting_angle: `numpy.float32`
            * relative_slope: `numpy.float32`
            * MODTRAN outputs: `numpy.float32`

        The acquisitions will be converted internally to int32 on a
        row by row basis.
    """
    # load the line starts and ends.
    rows, cols = geobox.shape

    # Specify the biliner binary files datatype
    boo_fnames = bilinear_ortho_filenames
    bilinear_dtype = 'float32'

    # Retrieve the satellite and sensor for the acquisition
    satellite = acquisitions[0].spacecraft_id
    sensor = acquisitions[0].sensor_id

    # Get the average reflectance values per band
    nbar_constants = constants.NBARConstants(satellite, sensor)
    avg_reflectance_values = nbar_constants.getAvgReflut()

    # Read arrays into memory
    # Convert to the appropriate datatype and transpose the array to convert
    # to Fortran contiguous memory. This should prevent any array copying
    self_shadow = as_array(read_image(self_shadow_fname), dtype=numpy.int16,
        transpose=True)
    cast_shadow_sun = as_array(read_image(cast_shadow_sun_fname),
        dtype=numpy.int16, transpose=True)
    cast_shadow_satellite = as_array(read_image(cast_shadow_satellite_fname),
        dtype=numpy.int16, transpose=True)
    solar_zenith = as_array(read_image(solar_zenith_fname),
        dtype=numpy.float32, transpose=True)
    solar_azimuth = as_array(read_image(solar_azimuth_fname),
        dtype=numpy.float32, transpose=True)
    satellite_view = as_array(read_image(satellite_view_fname),
        dtype=numpy.float32, transpose=True)
    relative_angle = as_array(read_image(relative_angle_fname),
        dtype=numpy.float32, transpose=True)
    slope = as_array(read_image(slope_fname), dtype=numpy.float32,
        transpose=True)
    aspect = as_array(read_image(aspect_fname), dtype=numpy.float32,
        transpose=True)
    incident_angle = as_array(read_image(incident_angle_fname),
        dtype=numpy.float32, transpose=True)
    exiting_angle = as_array(read_image(exiting_angle_fname),
        dtype=numpy.float32, transpose=True)
    relative_slope = as_array(read_image(relative_slope_fname),
        dtype=numpy.float32, transpose=True)

    # Loop over each acquisition and compute various reflectance arrays
    for acq in acquisitions:
        band_number = acq.band_num

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
            avg_reflectance_values[band_number])

        # Read the data; convert to required dtype and transpose
        band_data = make_array(acq.data(), dtype=numpy.int16, transpose=True)

        # Run terrain correction
	ref_lm, ref_brdf, ref_terrain = terrain_correction(
	    rori,
	    brdf0, brdf1, brdf2,
	    bias, slope_ca, esun, dd,
	    avg_reflectance_values[band_number],
	    band_data,
	    self_shadow,
	    cast_shadow_sun,
	    cast_shadow_satellite,
	    solar_zenith,
	    solar_azimuth,
	    satellite_view,
	    relative_angle,
	    slope,
	    aspect,
	    incident_angle,
	    exiting_angle,
	    relative_slope,
	    load_2D_bin_file(boo_fnames[(band_number, 'a')], rows, cols,
                dtype=bilinear_dtype, transpose=True),
	    load_2D_bin_file(boo_fnames[(band_number, 'b')], rows, cols,
                dtype=bilinear_dtype, transpose=True),
	    load_2D_bin_file(boo_fnames[(band_number, 's')], rows, cols,
                dtype=bilinear_dtype, transpose=True),
	    load_2D_bin_file(boo_fnames[(band_number, 'fs')], rows, cols,
                dtype=bilinear_dtype, transpose=True),
	    load_2D_bin_file(boo_fnames[(band_number, 'fv')], rows, cols,
                dtype=bilinear_dtype, transpose=True),
	    load_2D_bin_file(boo_fnames[(band_number, 'ts')], rows, cols,
                dtype=bilinear_dtype, transpose=True),
	    load_2D_bin_file(boo_fnames[(band_number, 'dir')], rows, cols,
                dtype=bilinear_dtype, transpose=True),
	    load_2D_bin_file(boo_fnames[(band_number, 'dif')], rows, cols,
                dtype=bilinear_dtype, transpose=True))


        # Filenames for lambertian, brdf & terrain corrected reflectance
        lmbrt_fname = reflectance_filenames[(band_number,
            'reflectance_lambertian')]
        brdf_fname = reflectance_filenames[(band_number, 'reflectance_brdf')]
        tc_fname = reflectance_filenames[(band_number, 'reflectance_terrain')]

        # Output the files, remember to transpose back
        write_img(ref_lm.transpose(), lmbrt_fname, geobox=geobox, nodata=-999)
        write_img(ref_brdf.transpose(), brdf_fname, geobox=geobox, nodata=-999)
        write_img(ref_terrain.transpose(), tc_fname, geobox=geobox, nodata=-999)

        # Free the memory
        ref_lm = ref_brdf = ref_terrain = None
        gc.collect()
