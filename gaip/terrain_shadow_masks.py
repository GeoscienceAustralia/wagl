#!/usr/bin/env python

"""
Functions for calculating cast shadow from both the sun and satellite
---------------------------------------------------------------------

as source directions, as well as self shadow masks.
---------------------------------------------------
"""

from __future__ import absolute_import, print_function
from posixpath import join as ppjoin
import numpy
import h5py

from gaip.constants import DatasetName, GroupName
from gaip.geobox import GriddedGeoBox
from gaip.margins import ImageMargins
from gaip.satellite_solar_angles import setup_spheroid
from gaip.hdf5 import dataset_compression_kwargs
from gaip.hdf5 import attach_image_attributes
from gaip.hdf5 import create_external_link
from gaip.tiling import generate_tiles
from gaip.__cast_shadow_mask import cast_shadow_main


def _self_shadow(incident_angles_fname, exiting_angles_fname, out_fname,
                 compression='lzf', y_tile=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(incident_angles_fname, 'r') as fid_incident,\
        h5py.File(exiting_angles_fname, 'r') as fid_exiting,\
        h5py.File(out_fname, 'w') as fid:

        grp1 = fid_incident[GroupName.incident_group.value]
        grp2 = fid_exiting[GroupName.exiting_group.value]
        self_shadow(grp1, grp2, fid, compression, y_tile)


def self_shadow(incident_angles_group, exiting_angles_group, out_group=None,
                compression='lzf', y_tile=None):
    """
    Computes the self shadow mask.

    :param incident_angles_group:
        The root HDF5 `Group` that contains the incident
        angle dataset specified by the pathname given by:

        * DatasetName.incident

    :param exiting_angles_group:
        The root HDF5 `Group` that contains the exiting
        angle dataset specified by the pathname given by:

        * DatasetName.exiting

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.
        The dataset name will be given by:

        * DatasetName.self_shadow

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param y_tile:
        Defines the tile size along the y-axis. Default is None which
        equates to all elements along the y-axis.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    incident_angle = incident_angles_group[DatasetName.incident.value]
    exiting_angle = exiting_angles_group[DatasetName.exiting.value]
    geobox = GriddedGeoBox.from_dataset(incident_angle)

    # Initialise the output file
    if out_group is None:
        fid = h5py.File('self-shadow.h5', driver='core', backing_store=False)
    else:
        fid = out_group

    grp = fid.create_group(GroupName.shadow_group.value)

    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))
    cols, rows = geobox.get_shape_xy()
    kwargs['shape'] = (rows, cols)
    kwargs['dtype'] = 'bool'

    # output dataset
    dataset_name = DatasetName.self_shadow.value
    out_dset = grp.create_dataset(dataset_name, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal()}
    desc = "Self shadow mask derived using the incident and exiting angles."
    attrs['Description'] = desc
    attach_image_attributes(out_dset, attrs)

    # Initialise the tiling scheme for processing
    tiles = generate_tiles(cols, rows, cols, y_tile)

    # Loop over each tile
    for tile in tiles:
        # Row and column start locations
        ystart, yend = tile[0]
        xstart, xend = tile[1]
        idx = (slice(ystart, yend), slice(xstart, xend))

        # Read the data for the current tile
        inc = numpy.radians(incident_angle[idx])
        exi = numpy.radians(exiting_angle[idx])

        # Process the tile
        mask = numpy.ones(inc.shape, dtype='uint8')
        mask[numpy.cos(inc) <= 0.0] = 0
        mask[numpy.cos(exi) <= 0.0] = 0

        # Write the current tile to disk
        out_dset[idx] = mask

    if out_group is None:
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
    Class that deals with errors from :py:func:`calculate_cast_shadow`.
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
                           out_fname, compression='lzf', y_tile=100,
                           solar_source=True):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(dsm_fname, 'r') as dsm_fid,\
        h5py.File(satellite_solar_angles_fname, 'r') as fid_sat_sol,\
        h5py.File(out_fname, 'w') as fid:

        grp1 = dsm_fid[GroupName.elevation_group.value]
        grp2 = fid_sat_sol[GroupName.sat_sol_group.value]
        calculate_cast_shadow(acquisition, grp1, grp2, margins, block_height,
                              block_width, fid, compression, y_tile,
                              solar_source)


def calculate_cast_shadow(acquisition, dsm_group, satellite_solar_group,
                          margins, block_height, block_width, out_group=None,
                          compression='lzf', y_tile=100, solar_source=True):
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

    :param dsm_group:
        The root HDF5 `Group` that contains the Digital Surface Model
        data.
        The dataset pathnames are given by:

        * DatasetName.dsm_smoothed

        The dataset must have the same dimensions as `acquisition`
        plus a margin of widths specified by margin.

    :param satellite_solar_group:
        The root HDF5 `Group` that contains the satellite and solar
        datasets specified by the pathnames given by:

        * DatasetName.solar_zenith
        * DatasetName.solar_azimuth
        * DatasetName.satellite_view
        * DatasetName.satellite_azimuth
        
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

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be given by the format string detailed
        by:

        * DatasetName.cast_shadow_fmt

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param y_tile:
        Defines the tile size along the y-axis. Default is 100.

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
    spheroid, _ = setup_spheroid(geobox.crs.ExportToWkt())

    # Define Top, Bottom, Left, Right pixel buffer margins
    pixel_buf = ImageMargins(margins)

    if solar_source:
        zenith_name = DatasetName.solar_zenith.value
        azimuth_name = DatasetName.solar_azimuth.value
    else:
        zenith_name = DatasetName.satellite_view.value
        azimuth_name = DatasetName.satellite_azimuth.value

    zenith_angle = satellite_solar_group[zenith_name][:]
    azimuth_angle = satellite_solar_group[azimuth_name][:]
    elevation = dsm_group[DatasetName.dsm_smoothed.value][:]

    # Compute the cast shadow mask
    ierr, mask = cast_shadow_main(elevation, zenith_angle, azimuth_angle,
                                  x_res, y_res, spheroid, y_origin, x_origin,
                                  pixel_buf.left, pixel_buf.right,
                                  pixel_buf.top, pixel_buf.bottom,
                                  block_height, block_width, is_utm)

    if ierr:
        raise CastShadowError(ierr)

    source_dir = 'sun' if solar_source else 'satellite'

    # Initialise the output file
    if out_group is None:
        fid = h5py.File('cast-shadow-{}.h5'.format(source_dir), driver='core',
                        backing_store=False)
    else:
        fid = out_group

    grp = fid.create_group(GroupName.shadow_group.value)
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))
    kwargs['dtype'] = 'bool'

    dname_fmt = DatasetName.cast_shadow_fmt.value
    out_dset = grp.create_dataset(dname_fmt.format(source=source_dir),
                                  data=mask, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal()}
    desc = ("The cast shadow mask determined using the {} "
            "as the source direction.").format(source_dir)
    attrs['Description'] = desc
    attach_image_attributes(out_dset, attrs)

    if out_group is None:
        return fid


def _combine_shadow(self_shadow_fname, cast_shadow_sun_fname,
                    cast_shadow_satellite_fname, out_fname,
                    compression='lzf', y_tile=100):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(self_shadow_fname, 'r') as fid_self,\
        h5py.File(cast_shadow_sun_fname, 'r') as fid_sun,\
        h5py.File(cast_shadow_satellite_fname, 'r') as fid_sat,\
        h5py.File(out_fname, 'w') as fid:

        grp1 = fid_self[GroupName.shadow_group.value]
        grp2 = fid_sun[GroupName.shadow_group.value]
        grp3 = fid_sat[GroupName.shadow_group.value]
        combine_shadow_masks(grp1, grp2, grp3, fid, compression, y_tile)

    link_shadow_datasets(self_shadow_fname, cast_shadow_sun_fname,
                         cast_shadow_satellite_fname, out_fname)


def combine_shadow_masks(self_shadow_group, cast_shadow_sun_group,
                         cast_shadow_satellite_group, out_group=None,
                         compression='lzf', y_tile=100):
    """
    A convienice function for combining the shadow masks into a single
    boolean array.

    :param self_shadow_group:
        The root HDF5 `Group` that contains the self shadow
        dataset specified by the pathname given by:

        * DatasetName.self_shadow

    :param cast_shadow_sun_group:
        The root HDF5 `Group` that contains the cast shadow
        (solar direction) dataset specified by the pathname
        given by:

        * DatasetName.cast_shadow_fmt

    :param cast_shadow_sun_group:
        The root HDF5 `Group` that contains the cast shadow
        (satellite direction) dataset specified by the pathname
        given by:

        * DatasetName.cast_shdadow_fmt

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be given by the format string detailed
        by:

        * DatasetName.combined_shadow

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
    # access the datasets
    dname_fmt = DatasetName.cast_shadow_fmt.value
    self_shad = self_shadow_group[DatasetName.self_shadow.value]
    cast_sun = cast_shadow_sun_group[dname_fmt.format(source='sun')]
    dname = dname_fmt.format(source='satellite')
    cast_sat = cast_shadow_satellite_group[dname]
    geobox = GriddedGeoBox.from_dataset(self_shad)

    # Initialise the output files
    if out_group is None:
        fid = h5py.File('combined-shadow.h5', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    grp = fid.create_group(GroupName.shadow_group.value)
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))
    cols, rows = geobox.get_shape_xy()
    kwargs['shape'] = (rows, cols)
    kwargs['dtype'] = 'bool'

    # output dataset
    out_dset = grp.create_dataset(DatasetName.combined_shadow.value, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal()}
    desc = ("Combined shadow masks: 1. self shadow, "
            "2. cast shadow (solar direction), "
            "3. cast shadow (satellite direction).")
    attrs['Description'] = desc
    attrs['mask_values'] = "False = Shadow; True = Non Shadow"
    attach_image_attributes(out_dset, attrs)

    # Initialise the tiling scheme for processing
    tiles = generate_tiles(cols, rows, cols, y_tile)

    # Loop over each tile
    for tile in tiles:
        # Row and column start locations
        ystart, yend = tile[0]
        xstart, xend = tile[1]
        idx = (slice(ystart, yend), slice(xstart, xend))

        out_dset[idx] = (self_shad[idx] & cast_sun[idx] & cast_sat[idx])

    if out_group is None:
        return fid


def link_shadow_datasets(self_shadow_fname, cast_shadow_sun_fname,
                         cast_shadow_satellite_fname, out_fname):
    """
    Link the self shadow mask, and the two cast shadow masks into a
    single file for easier access.
    """
    group_path = GroupName.shadow_group.value
    dname_fmt = DatasetName.cast_shadow_fmt.value
    dname = ppjoin(group_path, DatasetName.self_shadow.value)
    create_external_link(self_shadow_fname, dname, out_fname, dname)

    dname = ppjoin(group_path, dname_fmt.format(source='sun'))
    create_external_link(cast_shadow_sun_fname, dname, out_fname, dname)

    dname = ppjoin(group_path, dname_fmt.format(source='satellite'))
    create_external_link(cast_shadow_satellite_fname, dname, out_fname, dname)
