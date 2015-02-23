"""
Reflectance Calculations
------------------------
"""
import gc
import numpy

from gaip import as_array
from gaip import constants
from gaip import load_2d_bin_file
from gaip import read_img
from gaip import reflectance
from gaip import write_img
from gaip import write_new_brdf_file


def calculate_reflectance(acquisitions, bilinear_ortho_filenames, rori,
                          self_shadow_fname, cast_shadow_sun_fname,
                          cast_shadow_satellite_fname, solar_zenith_fname,
                          solar_azimuth_fname, satellite_view_fname,
                          relative_angle_fname, slope_fname, aspect_fname,
                          incident_angle_fname, exiting_angle_fname,
                          relative_slope_fname, reflectance_filenames,
                          brdf_fname_format, new_brdf_fname_format):
    """
    The workflow used to calculate lambertian, BRDF corrected and
    terrain corrected surface reflectance.

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
        Valid reflectance level strings are:

            * 1. ref_lm -> Lambertian reflectance
            * 2. ref_brdf -> BRDF corrected reflectance
            * 3. ref_terrain -> Terrain corrected reflectance

    :param brdf_fname_format:
        A string containing the brdf filename format eg:
        brdf_modis_band_{band_num}.txt, where {band_num} will be
        substituted for the current band number.

    :param new_brdf_fname_format:
        A string containing the new brdf filename format eg:
        new_brdf_modis_band_{band_num}.txt, where {band_num} will be
        substituted for the current band number.

    :return:
        None.
        The terrain correction algorithm will output 3 files for every
        band in the following format:

            * 1. reflectance_lambertian_{band_number}.bin -> Lambertian
                 reflectance.
            * 2. reflectance_brdf_{band_number}.bin -> BRDF corrected
                 reflectance.
            * 3. reflectance_terrain_{band_number}.bin -> Terrain corrected
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
    # Specify the biliner binary files datatype
    boo_fnames = bilinear_ortho_filenames
    bilinear_dtype = 'float32'

    # Retrieve the satellite and sensor for the acquisition
    satellite = acquisitions[0].spacecraft_id
    sensor = acquisitions[0].sensor_id

    # Get the average reflectance values per band
    nbar_constants = constants.NBARConstants(satellite, sensor)
    avg_reflectance_values = nbar_constants.get_avg_ref_lut()

    # Read all required arrays into memory
    # Convert to the appropriate datatype and transpose the array to convert
    # to Fortran contiguous memory. This should prevent any array copying
    self_shadow = as_array(read_img(self_shadow_fname), dtype=numpy.int16,
                           transpose=True)
    cast_shadow_sun = as_array(read_img(cast_shadow_sun_fname),
                               dtype=numpy.int16, transpose=True)
    cast_shadow_satellite = as_array(read_img(cast_shadow_satellite_fname),
                                     dtype=numpy.int16, transpose=True)
    solar_zenith = as_array(read_img(solar_zenith_fname),
                            dtype=numpy.float32, transpose=True)
    solar_azimuth = as_array(read_img(solar_azimuth_fname),
                             dtype=numpy.float32, transpose=True)
    satellite_view = as_array(read_img(satellite_view_fname),
                              dtype=numpy.float32, transpose=True)
    relative_angle = as_array(read_img(relative_angle_fname),
                              dtype=numpy.float32, transpose=True)
    slope = as_array(read_img(slope_fname), dtype=numpy.float32,
                     transpose=True)
    aspect = as_array(read_img(aspect_fname), dtype=numpy.float32,
                      transpose=True)
    incident_angle = as_array(read_img(incident_angle_fname),
                              dtype=numpy.float32, transpose=True)
    exiting_angle = as_array(read_img(exiting_angle_fname),
                             dtype=numpy.float32, transpose=True)
    relative_slope = as_array(read_img(relative_slope_fname),
                              dtype=numpy.float32, transpose=True)

    # Loop over each acquisition and compute various reflectance arrays
    for acq in acquisitions:
        rows = acq.lines
        cols = acq.samples
        band_number = acq.band_num
        geobox = acq.gridded_geo_box()

        # Read the BRDF modis file for a given band
        brdf_modis_file = brdf_fname_format.format(band_num=acq.band_num)
        with open(brdf_modis_file, 'r') as param_file:
            (brdf0, brdf1, brdf2, bias, slope_ca,
             esun, dd) = map(float, ' '.join(param_file.readlines()).split())

        # Output the new format of the brdf file (QA/QC)
        write_new_brdf_file(
            new_brdf_fname_format.format(band_num=band_number),
            rori, brdf0, brdf1, brdf2, bias, slope_ca, esun, dd,
            avg_reflectance_values[band_number])

        # Read the data; convert to required dtype and transpose
        band_data = as_array(acq.data(), dtype=numpy.int16, transpose=True)

        # Allocate the output arrays
        # They'll be transposed upon input to the function
        # The output and work arrays could be allocated once outside
        # of the loop. But for future proof, we may get a product where
        # each acquisition could have different resolutions.
        ref_lm = numpy.zeros((rows, cols), dtype='int16')
        ref_brdf = numpy.zeros((rows, cols), dtype='int16')
        ref_terrain = numpy.zeros((rows, cols), dtype='int16')

        # Allocate the work arrays (single row of data)
        band_work = numpy.zeros(cols, dtype='int32')
        ref_lm_work = numpy.zeros(cols, dtype='float32')
        ref_brdf_work = numpy.zeros(cols, dtype='float32')
        ref_terrain_work = numpy.zeros(cols, dtype='float32')

        # Read the bilinear ortho files for the current band
        a_mod = as_array(read_img(boo_fnames[(band_number, 'a')]),
                         dtype=numpy.float32, transpose=True)
        b_mod = as_array(read_img(boo_fnames[(band_number, 'b')]),
                         dtype=numpy.float32, transpose=True)
        s_mod = as_array(read_img(boo_fnames[(band_number, 's')]),
                         dtype=numpy.float32, transpose=True)
        fs = as_array(read_img(boo_fnames[(band_number, 'fs')]),
                      dtype=numpy.float32, transpose=True)
        fv = as_array(read_img(boo_fnames[(band_number, 'fv')]),
                      dtype=numpy.float32, transpose=True)
        ts = as_array(read_img(boo_fnames[(band_number, 'ts')]),
                      dtype=numpy.float32, transpose=True)
        edir_h = as_array(read_img(boo_fnames[(band_number, 'dir')]),
                          dtype=numpy.float32, transpose=True)
        edif_h = as_array(read_img(boo_fnames[(band_number, 'dif')]),
                          dtype=numpy.float32, transpose=True)

        # Run terrain correction
        # We use transposed arrays; rows become cols and cols become rows
        reflectance(cols, rows, rori, brdf0, brdf1, brdf2, bias, slope_ca,
                    esun, dd, avg_reflectance_values[band_number], band_data,
                    self_shadow, cast_shadow_sun, cast_shadow_satellite,
                    solar_zenith, solar_azimuth, satellite_view,
                    relative_angle, slope, aspect, incident_angle,
                    exiting_angle, relative_slope, a_mod, b_mod, s_mod, fs, fv,
                    ts, edir_h, edif_h, band_work, ref_lm_work, ref_brdf_work,
                    ref_terrain_work, ref_lm.transpose(), ref_brdf.transpose(),
                    ref_terrain.transpose())

        # Filenames for lambertian, brdf & terrain corrected reflectance
        lmbrt_fname = reflectance_filenames[(band_number,
                                             'reflectance_lambertian')]
        brdf_fname = reflectance_filenames[(band_number, 'reflectance_brdf')]
        tc_fname = reflectance_filenames[(band_number, 'reflectance_terrain')]

        # Output the files
        write_img(ref_lm, lmbrt_fname, geobox=geobox, nodata=-999)
        write_img(ref_brdf, brdf_fname, geobox=geobox, nodata=-999)
        write_img(ref_terrain, tc_fname, geobox=geobox, nodata=-999)

        # Free the memory
        ref_lm = ref_brdf = ref_terrain = None
        band_work = ref_lm_work = ref_brdf_work = ref_terrain_work = None
        gc.collect()
