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
from gaip import cast_shadow_main
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


class FortranError(Exception):

    """
    Base class for errors thrown from the Fortran code used in this module.
    """

    def __init__(self, function_name, code, msg):
        self.function_name = function_name
        self.code = code
        self.msg = msg or "Unknown error"

    def __str__(self):
        """
        Return a string representation of this Error.
        """
        err = "Error in Fotran code {0} (code {1}): {2}"
        err = err.format(self.function_name, self.code, self.msg)
        return err


class CastShadowError(FortranError):

    """
    Class that deals with errors from :py:func:`run_castshadow`.
    """

    def __init__(self, code):
        super(CastShadowError,
              self).__init__("cast_shadow_main",
                             code,
                             CastShadowError.get_error_message(code))

    @staticmethod
    def get_error_message(code):
        """
        Generate an error message for a specific code. It is OK for this have
        non-returning control paths, as this will results in ``None``, which
        is handled in the super class.
        """
        def tmpt(d, n):
            """Generate message."""
            err = "attempt to access invalid {0} of {1}".format(d, n)
            return err

        if code == 20:
            return tmpt('x', 'dem')
        if code == 21:
            return tmpt('x', 'dem_data')
        if code == 22:
            return tmpt('x', 'solar and sazi')
        if code == 23:
            return tmpt('x', 'solar_data')
        if code == 24:
            return tmpt('x', 'a')
        if code == 25:
            return tmpt('y', 'dem_data')
        if code == 26:
            return tmpt('y', 'a')
        if code == 27:
            return tmpt('x', 'mask_all')
        if code == 28:
            return tmpt('y', 'mask_all')
        if code == 29:
            return tmpt('x', 'mask')
        if code == 30:
            return tmpt('y', 'mask')
        if code == 31:
            return tmpt('X', 'dem and a')
        if code == 32:
            return tmpt('y', 'a')
        if code == 33:
            return tmpt('y', 'dem')
        if code == 34:
            return tmpt('x', 'mask_all')
        if code == 35:
            return tmpt('x', 'mask')
        if code == 36:
            return tmpt('y', 'mask_all')
        if code == 37:
            return tmpt('y', 'mask')
        if code == 38:
            return tmpt('x', 'dem')
        if code == 39:
            return tmpt('x', 'dem_data')
        if code == 40:
            return tmpt('x', 'solar')
        if code == 41:
            return tmpt('x', 'solar_data')
        if code == 42:
            return tmpt('x', 'a and dem')
        if code == 43:
            return tmpt('y', 'a')
        if code == 44:
            return tmpt('y', 'dem')
        if code == 45:
            return tmpt('x', 'mask_all')
        if code == 46:
            return tmpt('x', 'mask')
        if code == 47:
            return tmpt('y', 'mask_alll')
        if code == 48:
            return tmpt('y', 'mask')
        if code == 49:
            return tmpt('x', 'a and dem')
        if code == 50:
            return tmpt('y', 'a')
        if code == 51:
            return tmpt('y', 'dem')
        if code == 52:
            return tmpt('x', 'mask_all')
        if code == 53:
            return tmpt('x', 'mask')
        if code == 54:
            return tmpt('y', 'mask_all')
        if code == 55:
            return tmpt('y', 'mask')
        if code == 61:
            return "azimuth case not possible - phi_sun must be in 0 to 360 deg"
        if code == 62:
            return "k_max gt k_setting"
        if code == 63:
            return "add outside add_max ranges"
        if code == 71:
            return "Parameters defining A are invalid"
        if code == 72:
            return "Matrix A not embedded in image"
        if code == 73:
            return "matrix A does not have sufficient y margin"
        if code == 74:
            return "matrix A does not have sufficient x margin"


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
            zenith_name = 'solar-zenith'
            azimuth_name = 'solar-azimuth'
        else:
            zenith_name = 'satellite-view'
            azimuth_name = 'satellite-azimuth'

        view_dset = sat_sol[view_name]
        azi_dset = sat_sol[azimuth_name]

    fid = calculate_cast_shadow(acquisition, dsm_dset, margins, block_height,
                                block_width, zenith_dset, azi_dset, out_fname,
                                compression, solar_source)

    fid.close()
    return


def calculate_cast_shadow(acquisition, dsm_dataset, margins, block_height,
                          block_width, zenith_angle_dataset,
                          azimuth_angle_dataset, out_fname=None,
                          compression='lzf', solar_source=True):
    """
    This code is an interface to the fortran code
    cast_shadow_main.f90 written by Fuqin (and modified to
    work with F2py).

    The following was taken from the top of the Fotran program:
    "cast_shadow_main.f90":

    Creates a shadow mask for a standard Landsat scene
    the program was originally written by DLB Jupp in Oct. 2010
    for a small sub_matrix and was modified by Fuqin Li in Oct.
    2010 so that the program can be used for large landsat scene.

    Basically, a sub-matrix A is embedded in a larger DEM image
    and the borders must be large enough to find the shaded pixels.
    If we assume the solar azimuth and zenith angles change very
    little within the sub-matrix A, then the Landsat scene can be
    divided into several sub_matrix.
    For Australian region, with 0 .00025 degree resolution, the
    sub-marix A is set to 500x500

    we also need to set extra DEM lines/columns to run the Landsat
    scene (see parameter ``pix_buf``. This will change with elevation
    difference within the scene and solar zenith angle. For
    Australian region and Landsat scene with 0.00025 degree
    resolution, the maximum extra lines are set to 250 pixels/lines
    for each direction. This figure shold be sufficient for everywhere
    and anytime in Australia. Thus the DEM image will be larger than
    landsat image for 500 lines x 500 columns

    :param acquisition:
        An instance of an acquisition object.

    :param dsm_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the Digital Surface
        Model data when index/sliced.
        This must have the same dimensions as `acquisition`
        plus a margin of widths specified by margin

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

    :param zenith_angle_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the solar zenith
        or the satellite view angles when index/sliced.
        Must be of the same dimensions as `acquisition`.

    :param azimuth_angle_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the solar azimuth
        or the satellite azimuth angles when index/sliced.
        Must be of the same dimensions as `acquisition`.

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

    :warning:
        The Fortran code cannot be compiled with ``-O3`` as it
        produces incorrect results if it is.
    """
    # Setup the geobox
    geobox = acquisition.gridded_geo_box()
    x_res, y_res = geobox.pixelsize
    x_origin, y_origin = geobox.origin

    # Are we in UTM or geographics?
    is_utm = not geobox.crs.IsGeographic()

    # Retrive the spheroid parameters
    # (used in calculating pixel size in metres per lat/lon)
    spheroid = setup_spheroid(geobox.crs.ExportToWkt())

    # Read the dsm and angle arrays into memory
    dsm = dsm_dataset[:]
    zenith_angle = zenith_angle_dataset[:]
    azimuth_angle = azimuth_angle_dataset[:]

    # Define Top, Bottom, Left, Right pixel buffer margins
    pixel_buf = ImageMargins(margins)

    # Compute the cast shadow mask
    ierr, mask = cast_shadow_main(dsm, zenith_angle, azimuth_angle, x_res,
                                  y_res, spheroid, y_origin, x_origin,
                                  pixel_buf.left, pixel_buf.right,
                                  pixel_buf.top, pixel_buf.bottom,
                                  block_height, block_width, is_utm)

    if ierr:
        raise CastShadowError(ierr)

    source_dir = 'sun' if solar_source else 'satellite'

    # Initialise the output file
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
