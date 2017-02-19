#!/usr/bin/env python

"""
Functions for calculating cast shadow from both the sun and satellite
---------------------------------------------------------------------

as source directions, as well as self shadow masks.
---------------------------------------------------
"""

import numpy
import h5py

from eotools import tiling
from gaip import GriddedGeoBox
from gaip import ImageMargins
from gaip import setup_spheroid
from gaip import run_castshadow
from gaip import dataset_compression_kwargs
from gaip import attach_image_attributes


def _self_shadow(incident_angles_fname, exiting_angles_fname, out_fname,
                 compression='lzf', x_tile=None, y_tile=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(incident_angles_fname, 'r') as inci_angles,\
        h5py.File(exiting_angles_fname, 'r') as exit_angles:

        inci_dset = inci_angles['incident']
        exit_dset = exit_angles['exiting']

        shape = inci_dset.shape
        transform = inci_dset.attrs['geotransform']
        origin = (transform[0], transform[3])
        crs = inci_dset.attrs['crs_wkt']
        pixelsize = (abs(transform[1]), abs(transform[5]))
        geobox = GriddedGeoBox(shape, origin, pixelsize, crs)

        fid = self_shadow(inci_dset, exit_dset, geobox, out_fname, x_tile,
                          y_tile)

    fid.close()
    return


def self_shadow(incident_dataset, exiting_dataset, geobox, out_fname,
                compression='lzf', x_tile=None, y_tile=None):
    """
    Computes the self shadow mask.

    :param incident_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the incident
        angles when index/sliced.

    :param exiting_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the exiting
        angles when index/sliced.

    :param geobox:
        An instance of a `GriddedGeoBox` object.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset names will be as follows:

        * self-shadow

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
    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('self-shadow.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))
    cols, rows = geobox.get_shape_xy()
    kwargs['shape'] = (rows, cols)
    kwargs['dtype'] = 'bool'

    # output dataset
    out_dset = fid.create_dataset('self-shadow', **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.affine.to_gdal()}
    desc = "Self shadow mask derived using the incident and exiting angles."
    attrs['Description'] = desc
    attach_image_attributes(out_dset, attrs)

    # Initialise the tiling scheme for processing
    if x_tile is None:
        x_tile = cols
    if y_tile is None:
        y_tile = rows
    tiles = tiling.generate_tiles(cols, rows, x_tile, y_tile, generator=False)

    # Loop over each tile
    for tile in tiles:
        # Row and column start locations
        ystart, yend = tile[0]
        xstart, xend = tile[1]
        idx = (slice(ystart, yend, slice(xstart, xend)))

        # Read the data for the current tile
        inc = numpy.radians(incident_dataset[idx])
        exi = numpy.radians(exiting_dataset[idx])

        # Process the tile
        mask = numpy.ones(inc.shape, xsize), dtype='uint8')
        mask[numpy.cos(inc) <= 0.0] = 0
        mask[numpy.cos(exi) <= 0.0] = 0

        # Write the current tile to disk
        out_dset[idx] = mask

    fid.flush()
    return fid


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
    kwargs['dtype'] = 'bool'

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
