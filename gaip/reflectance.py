#!/usr/bin/env python

"""
Calculates the Lambertian, BRDF corrected and BRDF + Terrain corrected
----------------------------------------------------------------------

reflectance
-----------
"""

from __future__ import absolute_import, print_function
from posixpath import join as ppjoin
import numpy
import h5py

from gaip import constants
from gaip.constants import DatasetName, Model, GroupName
from gaip.data import as_array
from gaip.hdf5 import dataset_compression_kwargs
from gaip.hdf5 import attach_image_attributes
from gaip.hdf5 import create_external_link
from gaip.metadata import create_ard_yaml
from gaip.tiling import generate_tiles
from gaip.__surface_reflectance import reflectance


def _calculate_reflectance(acquisition, interpolation_fname,
                           satellite_solar_angles_fname, slope_aspect_fname,
                           relative_slope_fname, incident_angles_fname,
                           exiting_angles_fname, shadow_masks_fname,
                           ancillary_fname, rori, out_fname,
                           compression='lzf', y_tile=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(interpolation_fname, 'r') as fid_interp,\
        h5py.File(satellite_solar_angles_fname, 'r') as fid_sat_sol,\
        h5py.File(slope_aspect_fname, 'r') as fid_slp_asp,\
        h5py.File(relative_slope_fname, 'r') as fid_rel_slp,\
        h5py.File(incident_angles_fname, 'r') as fid_inc,\
        h5py.File(exiting_angles_fname, 'r') as fid_exi,\
        h5py.File(shadow_masks_fname, 'r') as fid_shadow,\
        h5py.File(ancillary_fname, 'r') as fid_anc,\
        h5py.File(out_fname, 'w') as fid:

        grp1 = fid_interp[GroupName.interp_group.value]
        grp2 = fid_sat_sol[GroupName.sat_sol_group.value]
        grp3 = fid_slp_asp[GroupName.slp_asp_group.value]
        grp4 = fid_rel_slp[GroupName.rel_slp_group.value]
        grp5 = fid_inc[GroupName.incident_group.value]
        grp6 = fid_exi[GroupName.exiting_group.value]
        grp7 = fid_shadow[GroupName.shadow_group.value]
        grp8 = fid_anc[GroupName.ancillary_group.value]
        calculate_reflectance(acquisition, grp1, grp2, grp3, grp4, grp5, grp6,
                              grp7, grp8, rori, fid, compression, y_tile)

        create_ard_yaml(acquisition, grp8, fid)


def calculate_reflectance(acquisition, interpolation_group,
                          satellite_solar_group, slope_aspect_group,
                          relative_slope_group, incident_angles_group,
                          exiting_angles_group, shadow_masks_group,
                          ancillary_group, rori, out_group=None,
                          compression='lzf', y_tile=100):
    """
    Calculates Lambertian, BRDF corrected and BRDF + terrain
    illumination corrected surface reflectance.

    :param acquisition:
        An instance of an acquisition object.

    :param interpolation_group:
        The root HDF5 `Group` that contains the interpolated
        atmospheric coefficients.
        The dataset pathnames are given by:

        * DatasetName.interpolation_fmt

    :param satellite_solar_group:
        The root HDF5 `Group` that contains the solar zenith and
        solar azimuth datasets specified by the pathnames given by:

        * DatasetName.solar_zenith
        * DatasetName.solar_azimuth
        * DatasetName.satellite_view
        * DatasetName.satellite_azimuth
        * DatasetName.relative_azimuth
        
    :param slope_aspect_group:
        The root HDF5 `Group` that contains the slope and aspect
        datasets specified by the pathnames given by:

        * DatasetName.slope
        * DatasetName.aspect

    :param relative_slope_group:
        The root HDF5 `Group` that contains the relative slope dataset
        specified by the pathname given by:

        * DatasetName.relative_slope

    :param incident_angles_group:
        The root HDF5 `Group` that contains the incident
        angle dataset specified by the pathname given by:

        * DatasetName.incident

    :param exiting_angles_group:
        The root HDF5 `Group` that contains the exiting
        angle dataset specified by the pathname given by:

        * DatasetName.exiting

    :param shadow_masks_group:
        The root HDF5 `Group` that contains the combined shadow
        masks; self shadow, cast shadow (solar),
        cast shadow (satellite), dataset specified by the pathname
        given by:

        * DatasetName.combined_shadow

    :param ancillary_group:
        The root HDF5 `Group` that contains the Isotropic (iso),
        RossThick (vol), and LiSparseR (geo) BRDF scalar coefficients.
        The dataset pathnames are given by:

        * DatasetName.brdf_fmt

    :param rori:
        Threshold for terrain correction. Fuqin to document.

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be given by the format string detailed
        by:

        * DatasetName.reflectance_fmt

        The reflectance products are:

        * lambertian
        * brdf
        * terrain (brdf + terrain illumination correction)

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
    bn = acq.band_num

    dname_fmt = DatasetName.interpolation_fmt.value
    fv_dataset = interpolation_group[dname_fmt.format(factor='fv', band=bn)]
    fs_dataset = interpolation_group[dname_fmt.format(factor='fs', band=bn)]
    b_dataset = interpolation_group[dname_fmt.format(factor='b', band=bn)]
    s_dataset = interpolation_group[dname_fmt.format(factor='s', band=bn)]
    a_dataset = interpolation_group[dname_fmt.format(factor='a', band=bn)]
    dir_dataset = interpolation_group[dname_fmt.format(factor='dir', band=bn)]
    dif_dataset = interpolation_group[dname_fmt.format(factor='dif', band=bn)]
    ts_dataset = interpolation_group[dname_fmt.format(factor='ts', band=bn)]
    solar_zenith_dset = satellite_solar_group[DatasetName.solar_zenith.value]
    solar_azimuth_dset = satellite_solar_group[DatasetName.solar_azimuth.value]
    satellite_v_dset = satellite_solar_group[DatasetName.satellite_view.value]
    relative_a_dset = satellite_solar_group[DatasetName.relative_azimuth.value]
    slope_dataset = slope_aspect_group[DatasetName.slope.value]
    aspect_dataset = slope_aspect_group[DatasetName.aspect.value]
    relative_s_dset = relative_slope_group[DatasetName.relative_slope.value]
    incident_angle_dataset = incident_angles_group[DatasetName.incident.value]
    exiting_angle_dataset = exiting_angles_group[DatasetName.exiting.value]
    shadow_dataset = shadow_masks_group[DatasetName.combined_shadow.value]

    dname_fmt = DatasetName.brdf_fmt.value
    brdf_iso = ancillary_group[dname_fmt.format(band=bn, factor='iso')][()]
    brdf_vol = ancillary_group[dname_fmt.format(band=bn, factor='vol')][()]
    brdf_geo = ancillary_group[dname_fmt.format(band=bn, factor='geo')][()]

    # Get the average reflectance values per band
    nbar_constants = constants.NBARConstants(acq.spacecraft_id, acq.sensor_id)
    avg_reflectance_values = nbar_constants.get_avg_ref_lut()

    # Initialise the output file
    if out_group is None:
        fid = h5py.File('surface-reflectance.h5', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    grp = fid.create_group(GroupName.standard_group.value)
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, acq.samples))
    kwargs['shape'] = (acq.lines, acq.samples)
    kwargs['fillvalue'] = -999
    kwargs['dtype'] = 'int16'

    # create the datasets
    dname_fmt = DatasetName.reflectance_fmt.value
    dname = dname_fmt.format(product='lambertian', band=bn)
    lmbrt_dset = grp.create_dataset(dname, **kwargs)

    dname = dname_fmt.format(product='brdf', band=bn)
    brdf_dset = grp.create_dataset(dname, **kwargs)

    dname = dname_fmt.format(product='terrain', band=bn)
    tc_dset = grp.create_dataset(dname, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal(),
             'no_data_value': kwargs['fillvalue'],
             'rori_threshold_setting': rori,
             'sattelite': acq.spacecraft_id,
             'sensor': acq.sensor_id,
             'band_number': bn}

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
        solar_zenith = as_array(solar_zenith_dset[idx], **f32_args)
        solar_azimuth = as_array(solar_azimuth_dset[idx], **f32_args)
        satellite_view = as_array(satellite_v_dset[idx], **f32_args)
        relative_angle = as_array(relative_a_dset[idx], **f32_args)
        slope = as_array(slope_dataset[idx], **f32_args)
        aspect = as_array(aspect_dataset[idx], **f32_args)
        incident_angle = as_array(incident_angle_dataset[idx], **f32_args)
        exiting_angle = as_array(exiting_angle_dataset[idx], **f32_args)
        relative_slope = as_array(relative_s_dset[idx], **f32_args)
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

    if out_group is None:
        return fid


def link_standard_data(input_fnames, out_fname, model):
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

    group_path = GroupName.standard_group.value
    for fname in input_fnames:
        with h5py.File(fname, 'r') as fid:
            base_group = fid[group_path]
            dataset_names = []
            for group in model.ard_products:
                if group not in base_group:
                    continue
                grp = base_group[group]
                dnames = [k for k, v in grp.items() if not exclude(v)]
                base_path = ppjoin(group_path, group)
                dataset_names.extend([ppjoin(base_path, d) for d in dnames])

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
