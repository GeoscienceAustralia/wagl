#!/usr/bin/env python

"""
Calculates the Lambertian, BRDF corrected and BRDF + Terrain corrected
----------------------------------------------------------------------

reflectance
-----------
"""

from __future__ import absolute_import, print_function
import numpy
import h5py

from gaip import constants
from gaip.constants import DatasetName
from gaip.data import as_array
from gaip.hdf5 import dataset_compression_kwargs
from gaip.hdf5 import attach_image_attributes
from gaip.hdf5 import create_external_link
from gaip.metadata import create_ard_yaml
from gaip.tiling import generate_tiles
from gaip.__surface_reflectance import reflectance


def _calculate_reflectance(acquisition, bilinear_fname,
                           satellite_solar_angles_fname, slope_aspect_fname,
                           relative_slope_fname, incident_angles_fname,
                           exiting_angles_fname, shadow_masks_fname,
                           ancillary_fname, rori, out_fname,
                           compression='lzf', y_tile=None):
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

        dname_fmt = DatasetName.interpolation_fmt.value
        fv_dset = fid_bil[dname_fmt.format(factor='fv', band=band_num)]
        fs_dset = fid_bil[dname_fmt.format(factor='fs', band=band_num)]
        b_dset = fid_bil[dname_fmt.format(factor='b', band=band_num)]
        s_dset = fid_bil[dname_fmt.format(factor='s', band=band_num)]
        a_dset = fid_bil[dname_fmt.format(factor='a', band=band_num)]
        dir_dset = fid_bil[dname_fmt.format(factor='dir', band=band_num)]
        dif_dset = fid_bil[dname_fmt.format(factor='dif', band=band_num)]
        ts_dset = fid_bil[dname_fmt.format(factor='ts', band=band_num)]
        sol_zen_dset = fid_sat_sol[DatasetName.solar_zenith.value]
        sol_azi_dset = fid_sat_sol[DatasetName.solar_azimuth.value]
        sat_view_dset = fid_sat_sol[DatasetName.satellite_view.value]
        rel_ang_dset = fid_sat_sol[DatasetName.relative_azimuth.value]
        slope_dset = fid_slp_asp[DatasetName.slope.value]
        aspect_dset = fid_slp_asp[DatasetName.aspect.value]
        rel_slp_dset = fid_rel_slp[DatasetName.relative_slope.value]
        inc_dset = fid_inc[DatasetName.incident.value]
        exi_dset = fid_exi[DatasetName.exiting.value]
        shad_dset = fid_shadow[DatasetName.combined_shadow.value]

        dname = DatasetName.brdf_fmt.value
        brdf_iso = fid_anc[dname.format(band=band_num, factor='iso')][()]
        brdf_vol = fid_anc[dname.format(band=band_num, factor='vol')][()]
        brdf_geo = fid_anc[dname.format(band=band_num, factor='geo')][()]

        fid = calculate_reflectance(acquisition, fv_dset, fs_dset, b_dset,
                                    s_dset, a_dset, dir_dset, dif_dset,
                                    ts_dset, sol_zen_dset, sol_azi_dset,
                                    sat_view_dset, rel_ang_dset, slope_dset,
                                    aspect_dset, rel_slp_dset, inc_dset,
                                    exi_dset, shad_dset, rori, brdf_iso,
                                    brdf_vol, brdf_geo, out_fname, compression,
                                    y_tile)

    create_ard_yaml(acquisition, ancillary_fname, fid)

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
                          compression='lzf', y_tile=100):
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

    :param y_tile:
        Defines the tile size along the y-axis. Default is 100.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    acq = acquisition
    geobox = acq.gridded_geo_box()

    # Get the average reflectance values per band
    nbar_constants = constants.NBARConstants(acq.spacecraft_id, acq.sensor_id)
    avg_reflectance_values = nbar_constants.get_avg_ref_lut()

    # Initialise the output file
    if out_fname is None:
        fid = h5py.File('surface-reflectance.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, acq.samples))
    kwargs['shape'] = (acq.lines, acq.samples)
    kwargs['fillvalue'] = -999
    kwargs['dtype'] = 'int16'

    # create the datasets
    dname_fmt = DatasetName.reflectance_fmt.value
    dname = dname_fmt.format(product='lambertian', band=acq.band_num)
    lmbrt_dset = fid.create_dataset(dname, **kwargs)

    dname = dname_fmt.format(product='brdf', band=acq.band_num)
    brdf_dset = fid.create_dataset(dname, **kwargs)

    dname = dname_fmt.format(product='terrain', band=acq.band_num)
    tc_dset = fid.create_dataset(dname, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal(),
             'no_data_value': kwargs['fillvalue'],
             'rori_threshold_setting': rori,
             'sattelite': acq.spacecraft_id,
             'sensor': acq.sensor_id,
             'band_number': acq.band_num}

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
    tiles = generate_tiles(acq.samples, acq.lines, acq.samples, y_tile)

    # Loop over each tile
    for tile in tiles:
        # tile indices
        idx = (slice(tile[0][0], tile[0][1]), slice(tile[1][0], tile[1][1]))

        # define some static arguments
        acq_args = {'window': tile,
                    'masked': False,
                    'apply_gain_offset': acq.scaled_radiance,
                    'out_no_data': kwargs['fillvalue']}
        f32_args = {'dtype': numpy.float32, 'transpose': True}

        # Read the data corresponding to the current tile for all dataset
        # Convert the datatype if required and transpose
        band_data = as_array(acq.data(**acq_args), **f32_args)
        
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
        xsize, ysize = band_data.shape # band_data has been transposed
        ref_lm = numpy.zeros((ysize, xsize), dtype='int16')
        ref_brdf = numpy.zeros((ysize, xsize), dtype='int16')
        ref_terrain = numpy.zeros((ysize, xsize), dtype='int16')

        # Allocate the work arrays (single row of data)
        ref_lm_work = numpy.zeros(xsize, dtype='float32')
        ref_brdf_work = numpy.zeros(xsize, dtype='float32')
        ref_terrain_work = numpy.zeros(xsize, dtype='float32')

        # Run terrain correction
        reflectance(xsize, ysize, rori, brdf_iso, brdf_vol, brdf_geo,
                    avg_reflectance_values[acq.band_num], kwargs['fillvalue'],
                    band_data, shadow, solar_zenith, solar_azimuth,
                    satellite_view, relative_angle, slope, aspect,
                    incident_angle, exiting_angle, relative_slope, a_mod,
                    b_mod, s_mod, fs, fv, ts, direct, diffuse, ref_lm_work,
                    ref_brdf_work, ref_terrain_work, ref_lm.transpose(),
                    ref_brdf.transpose(), ref_terrain.transpose())


        # Write the current tile to disk
        lmbrt_dset[idx] = ref_lm
        brdf_dset[idx] = ref_brdf
        tc_dset[idx] = ref_terrain

    return fid


def link_standard_data(input_fnames, out_fname):
    # TODO: incorporate linking for multi-granule and multi-group
    #       datasets
    """
    Links the individual reflectance and surface temperature
    results into a single file for easier access.
    """
    def exclude(obj):
        """
        A simple function to test an object against a
        h5py.Group object.
        """
        return isinstance(obj, h5py.Group)

    for fname in input_fnames:
        with h5py.File(fname, 'r') as fid:
            dataset_names = [k for k, v in fid.items() if not exclude(v)]

        for dname in dataset_names:
            if isinstance(dname, h5py.Group):
                continue
            create_external_link(fname, dname, out_fname, dname)

        # metadata
        with h5py.File(fname, 'r') as fid:
            with h5py.File(out_fname) as out_fid:
                yaml_dname = DatasetName.nbar_yaml.value
                if yaml_dname in fid and yaml_dname not in out_fid:
                    fid.copy(yaml_dname, out_fid, name=yaml_dname)

                yaml_dname = DatasetName.sbt_yaml.value
                if yaml_dname in fid and yaml_dname not in out_fid:
                    fid.copy(yaml_dname, out_fid, name=yaml_dname)
