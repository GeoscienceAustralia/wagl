"""
Reflectance Calculations
------------------------
"""
import numpy
import rasterio

from gaip import as_array
from gaip import constants
from gaip import reflectance
from eotools import tiling


def calculate_reflectance(acquisitions, bilinear_ortho_filenames, rori,
                          brdf_data,
                          self_shadow_fname, cast_shadow_sun_fname,
                          cast_shadow_satellite_fname, solar_zenith_fname,
                          solar_azimuth_fname, satellite_view_fname,
                          relative_angle_fname, slope_fname, aspect_fname,
                          incident_angle_fname, exiting_angle_fname,
                          relative_slope_fname, reflectance_filenames,
                          x_tile=None, y_tile=None):
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

    :param brdf_data:
        A dictionary with keys specified via a tuple of
        (band_number, factor), with valid factors being (iso, vol, geo).

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

    :param x_tile:
        Defines the tile size along the x-axis. Default is None which
        equates to all elements along the x-axis.

    :param y_tile:
        Defines the tile size along the y-axis. Default is None which
        equates to all elements along the y-axis.

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

    # Retrieve the satellite and sensor for the acquisition
    satellite = acquisitions[0].spacecraft_id
    sensor = acquisitions[0].sensor_id

    # Get the average reflectance values per band
    nbar_constants = constants.NBARConstants(satellite, sensor)
    avg_reflectance_values = nbar_constants.get_avg_ref_lut()

    # Open all the images for reading
    self_shadow_ds = rasterio.open(self_shadow_fname)
    cast_shadow_sun_ds = rasterio.open(cast_shadow_sun_fname)
    cast_shadow_satellite_ds = rasterio.open(cast_shadow_satellite_fname) 
    solar_zenith_ds = rasterio.open(solar_zenith_fname)
    solar_azimuth_ds = rasterio.open(solar_azimuth_fname)
    satellite_view_ds = rasterio.open(satellite_view_fname)
    relative_angle_ds = rasterio.open(relative_angle_fname)
    slope_ds = rasterio.open(slope_fname)
    aspect_ds = rasterio.open(aspect_fname)
    incident_angle_ds = rasterio.open(incident_angle_fname)
    exiting_angle_ds = rasterio.open(exiting_angle_fname)
    relative_slope_ds = rasterio.open(relative_slope_fname)

    # Loop over each acquisition and compute various reflectance arrays
    for acq in acquisitions:
        rows = acq.lines
        cols = acq.samples
        band_number = acq.band_num
        geobox = acq.gridded_geo_box()
        crs = geobox.crs.ExportToWkt()

        # Filenames for lambertian, brdf & terrain corrected reflectance
        lmbrt_fname = reflectance_filenames[(band_number,
                                             'reflectance_lambertian')]
        brdf_fname = reflectance_filenames[(band_number,
                                            'reflectance_brdf')]
        tc_fname = reflectance_filenames[(band_number,
                                          'reflectance_terrain')]

        # Initialise the output files
        kwargs = {'driver': 'GTiff',
                  'width': cols,
                  'height': rows,
                  'count': 1,
                  'crs': crs,
                  'transform': geobox.affine,
                  'dtype': 'int16',
                  'nodata': -999,
                  'compress': 'deflate',
                  'zlevel': 1,
                  'predictor': 2}
        outds_lmbrt = rasterio.open(lmbrt_fname, 'w', **kwargs)
        outds_brdf = rasterio.open(brdf_fname, 'w', **kwargs)
        outds_tc = rasterio.open(tc_fname, 'w', **kwargs)

        # Initialise the tiling scheme for processing
        if x_tile is None:
            x_tile = cols
        if y_tile is None:
            y_tile = rows
        tiles = tiling.generate_tiles(cols, rows, x_tile, y_tile,
                                      generator=False)

        # get the brdf values for each factor
        brdf0 = brdf_data[(band_number, 'iso')]['value']
        brdf1 = brdf_data[(band_number, 'vol')]['value']
        brdf2 = brdf_data[(band_number, 'geo')]['value']

        # Open all the bilinear interpolated files for the current band
        with rasterio.open(boo_fnames[(band_number, 'a')]) as a_mod_ds,\
            rasterio.open(boo_fnames[(band_number, 'b')]) as b_mod_ds,\
            rasterio.open(boo_fnames[(band_number, 's')]) as s_mod_ds,\
            rasterio.open(boo_fnames[(band_number, 'fs')]) as fs_ds,\
            rasterio.open(boo_fnames[(band_number, 'fv')]) as fv_ds,\
            rasterio.open(boo_fnames[(band_number, 'ts')]) as ts_ds,\
            rasterio.open(boo_fnames[(band_number, 'dir')]) as edir_h_ds,\
            rasterio.open(boo_fnames[(band_number, 'dif')]) as edif_h_ds:

            # Loop over each tile
            for tile in tiles:
                # Row and column start and end locations
                ystart = tile[0][0]
                xstart = tile[1][0]
                yend = tile[0][1]
                xend = tile[1][1]

                # Tile size
                ysize = yend - ystart
                xsize = xend - xstart

                # define some static arguments
                acq_args = {'window': tile,
                            'masked': False,
                            'apply_gain_offset': acq.scaled_radiance}
                args = {'window': tile,
                        'masked': False}
                i16_args = {'dtype': numpy.int16,
                            'transpose': True}
                f32_args = {'dtype': numpy.float32,
                            'transpose': True}

                # Read the data corresponding to the current tile for all
                # files
                # Convert the datatype if required and transpose
                band_data = as_array(acq.data(**acq_args), **f32_args)
                
                self_shadow = as_array(self_shadow_ds.read(1, **args),
                                       **i16_args)
                cast_shadow_sun = as_array(cast_shadow_sun_ds.read(1, **args),
                                           **i16_args)
                cast_shadow_satellite = as_array(
                    cast_shadow_satellite_ds.read(1, **args), **i16_args)
                solar_zenith = as_array(solar_zenith_ds.read(1, **args),
                                        **f32_args)
                solar_azimuth = as_array(solar_azimuth_ds.read(1, **args),
                                         **f32_args)
                satellite_view = as_array(satellite_view_ds.read(1, **args),
                                          **f32_args)
                relative_angle = as_array(relative_angle_ds.read(1, **args),
                                          **f32_args)
                slope = as_array(slope_ds.read(1, **args), **f32_args)
                aspect = as_array(aspect_ds.read(1, **args), **f32_args)
                incident_angle = as_array(incident_angle_ds.read(1, **args),
                                          **f32_args)
                exiting_angle = as_array(exiting_angle_ds.read(1, **args),
                                         **f32_args)
                relative_slope = as_array(relative_slope_ds.read(1, **args),
                                          **f32_args)
                a_mod = as_array(a_mod_ds.read(1, **args), **f32_args)
                b_mod = as_array(b_mod_ds.read(1, **args), **f32_args)
                s_mod = as_array(s_mod_ds.read(1, **args), **f32_args)
                fs = as_array(fs_ds.read(1, **args), **f32_args)
                fv = as_array(fv_ds.read(1, **args), **f32_args)
                ts = as_array(ts_ds.read(1, **args), **f32_args)
                edir_h = as_array(edir_h_ds.read(1, **args), **f32_args)
                edif_h = as_array(edif_h_ds.read(1, **args), **f32_args)

                # Allocate the output arrays
                # The output and work arrays could be allocated once
                # outside of the loop. But for future proof, we may get a
                # product where each acquisition could have different
                # resolutions.
                ref_lm = numpy.zeros((ysize, xsize), dtype='int16')
                ref_brdf = numpy.zeros((ysize, xsize), dtype='int16')
                ref_terrain = numpy.zeros((ysize, xsize), dtype='int16')

                # Allocate the work arrays (single row of data)
                ref_lm_work = numpy.zeros(xsize, dtype='float32')
                ref_brdf_work = numpy.zeros(xsize, dtype='float32')
                ref_terrain_work = numpy.zeros(xsize, dtype='float32')

                # Run terrain correction
                reflectance(xsize, ysize, rori, brdf0, brdf1, brdf2,
                            avg_reflectance_values[band_number],
                            band_data, self_shadow, cast_shadow_sun,
                            cast_shadow_satellite,
                            solar_zenith, solar_azimuth, satellite_view,
                            relative_angle, slope, aspect, incident_angle,
                            exiting_angle, relative_slope, a_mod, b_mod,
                            s_mod, fs, fv, ts, edir_h, edif_h,
                            ref_lm_work, ref_brdf_work, ref_terrain_work,
                            ref_lm.transpose(), ref_brdf.transpose(),
                            ref_terrain.transpose())


                # Write the current tile to disk
                outds_lmbrt.write(ref_lm, 1, window=tile)
                outds_brdf.write(ref_brdf, 1, window=tile)
                outds_tc.write(ref_terrain, 1, window=tile)

        # Close the files to complete the writing
        outds_lmbrt.close()
        outds_brdf.close()
        outds_tc.close()

    # close all the opened image files
    self_shadow_ds.close()
    cast_shadow_sun_ds.close()
    cast_shadow_satellite_ds.close()
    solar_zenith_ds.close()
    solar_azimuth_ds.close()
    satellite_view_ds.close()
    relative_angle_ds.close()
    slope_ds.close()
    aspect_ds.close()
    incident_angle_ds.close()
    exiting_angle_ds.close()
    relative_slope_ds.close()
