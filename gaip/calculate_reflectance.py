"""
Calculates the Lambertian, BRDF corrected and BRDF + Terrain corrected
----------------------------------------------------------------------

reflectance
-----------
"""

from __future__ import absolute_import
import numpy
import h5py

from gaip import constants
from gaip.data import as_array
from gaip.hdf5 import dataset_compression_kwargs
from gaip.hdf5 import attach_image_attributes
from gaip.hdf5 import create_external_link
from gaip.tiling import generate_tiles
from gaip.__surface_reflectance import reflectance


DATASET_NAME_FMT = '{product}-reflectance-band-{band}'


def _calculate_reflectance(acquisition, bilinear_fname,
                           satellite_solar_angles_fname, slope_aspect_fname,
                           relative_slope_fname, incident_angles_fname,
                           exiting_angles_fname, shadow_masks_fname,
                           ancillary_fname, rori, out_fname,
                           compression='lzf', x_tile=None, y_tile=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    band_num = acquisition.band_num
    with h5py.File(bilinear_fname, 'r') as fid_bil,\
        h5py.File(satellite_solar_angles_fname, 'r') as fid_sat_sol,\
        h5py.File(slope_aspect_fname, 'r') as fid_slp_asp,\
        h5py.File(relative_slope_fname, 'r') as fid_rel_slp,\
        h5py.File(incident_angles_fname, 'r') as fid_inc,\
        h5py.File(exiting_angles_fname, 'r') as fid_exi,\
        h5py.File(shadow_masks_fname, 'r') as fid_shadow,\
        h5py.File(ancillary_fname, 'r') as fid_anc:

        fv_dset = fid_bil['fv-band-{band}'.format(band=band_num)]
        fs_dset = fid_bil['fs-band-{band}'.format(band=band_num)]
        b_dset = fid_bil['b-band-{band}'.format(band=band_num)]
        s_dset = fid_bil['s-band-{band}'.format(band=band_num)]
        a_dset = fid_bil['a-band-{band}'.format(band=band_num)]
        dir_dset = fid_bil['dir-band-{band}'.format(band=band_num)]
        dif_dset = fid_bil['dif-band-{band}'.format(band=band_num)]
        ts_dset = fid_bil['ts-band-{band}'.format(band=band_num)]
        sol_zen_dset = fid_sat_sol['solar-zenith']
        sol_azi_dset = fid_sat_sol['solar-azimuth']
        sat_view_dset = fid_sat_sol['satellite-view']
        rel_ang_dset = fid_sat_sol['relative-azimuth']
        slope_dset = fid_slp_asp['slope']
        aspect_dset = fid_slp_asp['aspect']
        rel_slp_dset = fid_rel_slp['relative-slope']
        inc_dset = fid_inc['incident-angle']
        exi_dset = fid_exi['exiting-angle']
        shad_dset = fid_shadow['combined-shadow']

        dname = "BRDF-Band-{band}-{factor}"
        brdf_iso = fid_anc[dname.format(band=band_num, factor='iso')][()]
        brdf_vol = fid_anc[dname.format(band=band_num, factor='vol')][()]
        brdf_geo = fid_anc[dname.format(band=band_num, factor='geo')][()]

    fid = calculate_reflectance(acquisition, fv_dset, fs_dset, b_dset, s_dset,
                                a_dset, dir_dset, dif_dset, ts_dset,
                                sol_zen_dset, sol_azi_dset, sat_view_dset,
                                rel_ang_dset, slope_dset, aspect_dset,
                                rel_slp_dset, inc_dset, exi_dset, shad_dset,
                                rori, brdf_iso, brdf_vol, brdf_geo, out_fname,
                                compression, x_tile, y_tile)

    fid.close()
    return


def calculate_reflectance(acquisition, fv_dataset, fs_dataset, b_dataset,
                          s_dataset, a_dataset, dir_dataset, dif_dataset,
                          ts_dataset, solar_zenith_dataset,
                          solar_azimuth_dataset, satellite_view_dataset,
                          relative_angle_dataset, slope_dataset,
                          aspect_dataset, relative_slope_dataset,
                          incident_angle_dataset, exiting_angle_dataset,
                          shadow_dataset, rori, brdf_iso,
                          brdf_vol, brdf_geo, out_fname=None,
                          compression='lzf', x_tile=None, y_tile=None):
    """
    Calculates Lambertian, BRDF corrected and BRDF + terrain
    corrected surface reflectance.

    :param acquisition:
        An instance of an acquisition object.

    :param fv_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the MODTRAN
        factor `fv` data values when indexed/sliced.

    :param fs_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the MODTRAN
        factor `fs` data values when indexed/sliced.

    :param b_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the MODTRAN
        factor `b` data values when indexed/sliced.

    :param s_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the MODTRAN
        factor `s` data values when indexed/sliced.

    :param a_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the MODTRAN
        factor `a` data values when indexed/sliced.

    :param dir_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the MODTRAN
        factor `dir` (direct irradiance) data values when
        indexed/sliced.

    :param dif_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the MODTRAN
        factor `dif` (diffuse irradiance) data values when
        indexed/sliced.

    :param ts_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the MODTRAN
        factor `ts` data values when indexed/sliced.

    :param solar_zenith_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the solar
        zenith angles when indexed/sliced.

    :param solar_azimuth_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the solar
        azimuth angles when indexed/sliced.

    :param satellite_view_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the satellite
        view angles when indexed/sliced.

    :param relative_angle_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the relative
        azimuth angles when indexed/sliced.

    :param slope_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the slope values
        when index/sliced.

    :param aspect_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the aspect angles
        when index/sliced.

    :param relative_slope_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the relative
        slope values when index/sliced.

    :param incident_angle_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the incident angles
        when index/sliced.

    :param exiting_angle_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the exiting angles
        when index/sliced.

    :param shadow_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the shadow
        mask when index/sliced.

    :param rori:
        Threshold for terrain correction. Fuqin to document.

    :param brdf_iso:
        A point value representing the Bidirectional Reflectance
        Distribution Function for the isometric scattering fraction.

    :param brdf_vol:
        A point value representing the Bidirectional Reflectance
        Distribution Function for the volumetric scattering fraction.

    :param brdf_geo:
        A point value representing the Bidirectional Reflectance
        Distribution Function for the geometric scattering fraction.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset names will be as follows:

        * lambertian-reflectance-band-{number}
        * brdf-reflectance-band-{number}
        * terrain-reflectance-band-{number}

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param x_tile:
        Defines the tile size along the x-axis. Default is None which
        equates to all elements along the x-axis.

    :param y_tile:
        Defines the tile size along the y-axis. Default is None which
        equates to all elements along the y-axis.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    # Retrieve info about the acquisition
    rows = acquisition.lines
    cols = acquisition.samples
    geobox = acquisition.gridded_geo_box()

    # Get the average reflectance values per band
    nbar_constants = constants.NBARConstants(acquisition.spacecraft_id,
                                             acquisition.sensor_id)
    avg_reflectance_values = nbar_constants.get_avg_ref_lut()

    # Initialise the output file
    if out_fname is None:
        fid = h5py.File('satellite-solar-angles.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, cols))
    kwargs['shape'] = (rows, cols)
    kwargs['fillvalue'] = -999
    kwargs['dtype'] = 'int16'

    # create the datasets
    dataset_name = DATASET_NAME_FMT.format(product='lambertian',
                                           band=acquisition.band_num)
    lmbrt_dset = fid.create_dataset(dataset_name, **kwargs)

    dataset_name = DATASET_NAME_FMT.format(product='brdf',
                                           band=acquisition.band_num)
    brdf_dset = fid.create_dataset(dataset_name, **kwargs)

    dataset_name = DATASET_NAME_FMT.format(product='terrain',
                                           band=acquisition.band_num)
    tc_dset = fid.create_dataset(dataset_name, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.affine.to_gdal(),
             'no_data_value': kwargs['fillvalue'],
             'rori threshold setting': rori,
             'sattelite': acquisition.spacecraft_id,
             'sensor': acquisition.sensor_id,
             'band number': acquisition.band_num}

    desc = "Contains the lambertian reflectance data scaled by 10000."
    attrs['Description'] = desc
    attach_image_attributes(lmbrt_dset, attrs)

    desc = "Contains the brdf corrected reflectance data scaled by 10000."
    attrs['Description'] = desc
    attach_image_attributes(brdf_dset, attrs)

    desc = ("Contains the brdf and terrain corrected reflectance data scaled "
            "by 10000.")
    attrs['Description'] = desc
    attach_image_attributes(tc_dset, attrs)

    # Initialise the tiling scheme for processing
    if x_tile is None:
        x_tile = cols
    if y_tile is None:
        y_tile = rows
    tiles = generate_tiles(cols, rows, x_tile, y_tile)

    # Loop over each tile
    for tile in tiles:
        # tile indices
        idx = (slice(tile[0][0], tile[0][1]), slice(tile[1][0], tile[1][1]))

        # define some static arguments
        acq_args = {'window': tile,
                    'masked': False,
                    'apply_gain_offset': acquisition.scaled_radiance,
                    'out_no_data': kwargs['fillvalue']}
        i16_args = {'dtype': numpy.int16, 'transpose': True}
        f32_args = {'dtype': numpy.float32, 'transpose': True}

        # Read the data corresponding to the current tile for all dataset
        # Convert the datatype if required and transpose
        band_data = as_array(acquisition.data(**acq_args), **f32_args)
        
        shadow = as_array(shadow_dataset[idx], numpy.int8, transpose=True)
        solar_zenith = as_array(solar_zenith_dataset[idx], **f32_args)
        solar_azimuth = as_array(solar_azimuth_dataset[idx], **f32_args)
        satellite_view = as_array(satellite_view_dataset[idx], **f32_args)
        relative_angle = as_array(relative_angle_dataset[idx], **f32_args)
        slope = as_array(slope_dataset[idx], **f32_args)
        aspect = as_array(aspect_dataset[idx], **f32_args)
        incident_angle = as_array(incident_angle_dataset[idx], **f32_args)
        exiting_angle = as_array(exiting_angle_dataset[idx], **f32_args)
        relative_slope = as_array(relative_slope_dataset[idx], **f32_args)
        a_mod = as_array(a_dataset[idx], **f32_args)
        b_mod = as_array(b_dataset[idx], **f32_args)
        s_mod = as_array(s_dataset[idx], **f32_args)
        fs = as_array(fs_dataset[idx], **f32_args)
        fv = as_array(fv_dataset[idx], **f32_args)
        ts = as_array(ts_dataset[idx], **f32_args)
        direct = as_array(dir_dataset[idx], **f32_args)
        diffuse = as_array(dif_dataset[idx], **f32_args)

        # Allocate the output arrays
        ysize, xsize = band_data.shape
        ref_lm = numpy.zeros((ysize, xsize), dtype='int16')
        ref_brdf = numpy.zeros((ysize, xsize), dtype='int16')
        ref_terrain = numpy.zeros((ysize, xsize), dtype='int16')

        # Allocate the work arrays (single row of data)
        ref_lm_work = numpy.zeros(xsize, dtype='float32')
        ref_brdf_work = numpy.zeros(xsize, dtype='float32')
        ref_terrain_work = numpy.zeros(xsize, dtype='float32')

        # Run terrain correction
        reflectance(xsize, ysize, rori, brdf_iso, brdf_vol, brdf_geo,
                    avg_reflectance_values[acquisition.band_num],
                    kwargs['fillvalue'], band_data, shadow,
                    solar_zenith, solar_azimuth, satellite_view,
                    relative_angle, slope, aspect, incident_angle,
                    exiting_angle, relative_slope, a_mod, b_mod,
                    s_mod, fs, fv, ts, direct, diffuse,
                    ref_lm_work, ref_brdf_work, ref_terrain_work,
                    ref_lm.transpose(), ref_brdf.transpose(),
                    ref_terrain.transpose())


        # Write the current tile to disk
        lmbrt_dset[idx] = ref_lm
        brdf_dset[idx] = ref_brdf
        tc_dset[idx] = ref_terrain

    return fid


def link_reflectance_data(input_fnames, out_fname):
    # TODO: incorporate linking for multi-granule and multi-group
    #       datasets
    """
    Links the individual reflectance results into a
    single file for easier access.
    """
    for fname in input_fnames:
        with h5py.File(fname, 'r') as fid:
            dataset_names = fid.keys()
        for dname in dataset_names:
            create_external_link(fname, dname, out_fname, dname)
