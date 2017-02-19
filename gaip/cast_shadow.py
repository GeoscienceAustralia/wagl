#!/usr/bin/env python

"""
Cast shadow calculation
-----------------------
"""

import h5py
from gaip import ImageMargins
from gaip import setup_spheroid
from gaip import run_castshadow
from gaip import dataset_compression_kwargs
from gaip import attach_image_attributes


def _calculate_cast_shadow(acquisition, dsm_fname, margins, block_height,
                           block_width, satellite_solar_angles_fname,
                           out_fname, compression='lzf', solar_source=True):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(dsm_fname, 'r') as dsm_src,\
        h5py.File(satellite_solar_angles_fname, 'r') as sat_sol:

        dsm_dset = dsm_src['dsm-smoothed']

        if solar_source:
            view_name = 'solar-zenith'
            azimuth_name = 'solar-azimuth'
        else:
            view_name = 'satellite-view'
            azimuth_name = 'satellite-azimuth'

        view_dset = sat_sol[view_name]
        azi_dset = sat_sol[azimuth_name]

    fid = calculate_cast_shadow(acquisition, dsm_dset, margins, block_height,
                                block_width, view_dset, azi_dset, out_fname,
                                compression, solar_source)

    fid.close()
    return


def calculate_cast_shadow(acquisition, dsm_dataset, margins, block_height,
                          block_width, view_angle_dataset,
                          azimuth_angle_dataset, out_fname=None,
                          compression='lzf', solar_source=True):
    """
    :param acquisition:
        An instance of an acquisition object.

    :param dsm_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the Digital Surface
        Model data when index/sliced.

    :param margins:
        An object with members top, bottom, left and right giving the
        size of the margins (in pixels) which have been added to the
        corresponding sides of dsm.

    :param block_height:
        The height (rows) of the window/submatrix used in the cast
        shadow algorithm.

    :param block_width:
        The width (rows) of the window/submatrix used in the cast
        shadow algorithm.

    :param view_angle_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the solar zenith
        or the satellite view angles when index/sliced.

    :param azimuth_angle_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the solar azimuth
        or the satellite azimuth angles when index/sliced.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset names will be as follows:

        * cast-shadow-{source} where source is either the sun or satellite

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param solar_source:
        A `bool` indicating whether or not the source for the line
        of sight comes from the sun (True; Default), or False
        indicating the satellite.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    # Setup the geobox
    geobox = acquisition.gridded_geo_box()

    # Retrive the spheroid parameters
    # (used in calculating pixel size in metres per lat/lon)
    spheroid = setup_spheroid(geobox.crs.ExportToWkt())

    # Read the dsm and angle arrays into memory
    dsm = dsm_dataset[:]
    view_angle = view_angle_dataset[:]
    azimuth_angle = azimuth_angle_dataset[:]

    # Define Top, Bottom, Left, Right pixel margins
    pixel_buf = ImageMargins(margins)

    # Compute the cast shadow mask
    mask = run_castshadow(acquisition, dsm, view_angle, azimuth_angle,
                          pixel_buf, block_height, block_width, spheroid)

    source_dir = 'sun' if solar_source else 'satellite'

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('cast-shadow-{}.h5'.format(source_dir), driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))
    no_data = -999
    kwargs['fillvalue'] = no_data
    kwargs['no_data_value'] = no_data

    out_dset = fid.create_dataset('cast-shadow-{}'.format(source_dir),
                                  data=mask, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.affine.to_gdal()}
    desc = ("The cast shadow mask determined using the {} "
            "as the source direction.").format(source_dir)
    attrs['Description'] = desc
    attach_image_attributes(out_dset, attrs)

    fid.flush()
    return fid
