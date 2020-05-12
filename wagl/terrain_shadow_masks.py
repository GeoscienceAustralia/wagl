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

from wagl.constants import DatasetName, GroupName
from wagl.geobox import GriddedGeoBox
from wagl.margins import pixel_buffer
from wagl.satellite_solar_angles import setup_spheroid
from wagl.hdf5 import H5CompressionFilter, attach_image_attributes
from wagl.hdf5 import create_external_link
from wagl.tiling import generate_tiles
from wagl.__cast_shadow_mask import cast_shadow_main


def _self_shadow(incident_angles_fname, exiting_angles_fname, out_fname,
                 compression=H5CompressionFilter.LZF, filter_opts=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(incident_angles_fname, 'r') as fid_incident,\
        h5py.File(exiting_angles_fname, 'r') as fid_exiting,\
        h5py.File(out_fname, 'w') as fid:

        grp1 = fid_incident[GroupName.INCIDENT_GROUP.value]
        grp2 = fid_exiting[GroupName.EXITING_GROUP.value]
        self_shadow(grp1, grp2, fid, compression, filter_opts)


def self_shadow(incident_angles_group, exiting_angles_group, out_group=None,
                compression=H5CompressionFilter.LZF, filter_opts=None):
    """
    Computes the self shadow mask.

    :param incident_angles_group:
        The root HDF5 `Group` that contains the incident
        angle dataset specified by the pathname given by:

        * DatasetName.INCIDENT

    :param exiting_angles_group:
        The root HDF5 `Group` that contains the exiting
        angle dataset specified by the pathname given by:

        * DatasetName.EXITING

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.
        The dataset name will be given by:

        * DatasetName.SELF_SHADOW

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF 

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    incident_angle = incident_angles_group[DatasetName.INCIDENT.value]
    exiting_angle = exiting_angles_group[DatasetName.EXITING.value]
    geobox = GriddedGeoBox.from_dataset(incident_angle)

    # Initialise the output file
    if out_group is None:
        fid = h5py.File('self-shadow.h5', driver='core', backing_store=False)
    else:
        fid = out_group

    if filter_opts is None:
        filter_opts = {}
    else:
        filter_opts = filter_opts.copy()

    if GroupName.SHADOW_GROUP.value not in fid:
        fid.create_group(GroupName.SHADOW_GROUP.value)

    grp = fid[GroupName.SHADOW_GROUP.value]

    tile_size = exiting_angle.chunks
    filter_opts['chunks'] = tile_size
    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    cols, rows = geobox.get_shape_xy()
    kwargs['shape'] = (rows, cols)
    kwargs['dtype'] = 'bool'

    # output dataset
    dataset_name = DatasetName.SELF_SHADOW.value
    out_dset = grp.create_dataset(dataset_name, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal()}
    desc = "Self shadow mask derived using the incident and exiting angles."
    attrs['description'] = desc
    attrs['alias'] = 'self-shadow'
    attach_image_attributes(out_dset, attrs)

    # process by tile
    for tile in generate_tiles(cols, rows, tile_size[1], tile_size[0]):
        # Read the data for the current tile
        inc = numpy.radians(incident_angle[tile])
        exi = numpy.radians(exiting_angle[tile])

        # Process the tile
        mask = numpy.ones(inc.shape, dtype='uint8')
        mask[numpy.cos(inc) <= 0.0] = 0
        mask[numpy.cos(exi) <= 0.0] = 0

        # Write the current tile to disk
        out_dset[tile] = mask

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


def _calculate_cast_shadow(acquisition, dsm_fname, buffer_distance,
                           satellite_solar_angles_fname, out_fname,
                           compression=H5CompressionFilter.LZF,
                           filter_opts=None, solar_source=True):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(dsm_fname, 'r') as dsm_fid,\
        h5py.File(satellite_solar_angles_fname, 'r') as fid_sat_sol,\
        h5py.File(out_fname, 'w') as fid:

        grp1 = dsm_fid[GroupName.ELEVATION_GROUP.value]
        grp2 = fid_sat_sol[GroupName.SAT_SOL_GROUP.value]
        calculate_cast_shadow(acquisition, grp1, grp2, buffer_distance, fid,
                              compression, filter_opts, solar_source)


def calculate_cast_shadow(acquisition, dsm_group, satellite_solar_group,
                          buffer_distance, out_group=None,
                          compression=H5CompressionFilter.LZF,
                          filter_opts=None, solar_source=True):
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
    scene. This will change with elevation
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

        * DatasetName.DSM_SMOOTHED

        The dataset must have the same dimensions as `acquisition`
        plus a margin of widths specified by margin.

    :param satellite_solar_group:
        The root HDF5 `Group` that contains the satellite and solar
        datasets specified by the pathnames given by:

        * DatasetName.SOLAR_ZENITH
        * DatasetName.SOLAR_AZIMUTH
        * DatasetName.SATELLITE_VIEW
        * DatasetName.SATELLITE_AZIMUTH

    :param buffer_distance:
        A number representing the desired distance (in the same
        units as the acquisition) in which to calculate the extra
        number of pixels required to buffer an image.
        Default is 8000.

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be given by the format string detailed
        by:

        * DatasetName.CAST_SHADOW_FMT

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF 

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

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
    margins = pixel_buffer(acquisition, buffer_distance)

    if solar_source:
        zenith_name = DatasetName.SOLAR_ZENITH.value
        azimuth_name = DatasetName.SOLAR_AZIMUTH.value
    else:
        zenith_name = DatasetName.SATELLITE_VIEW.value
        azimuth_name = DatasetName.SATELLITE_AZIMUTH.value

    zenith_angle = satellite_solar_group[zenith_name][:]
    azimuth_angle = satellite_solar_group[azimuth_name][:]
    elevation = dsm_group[DatasetName.DSM_SMOOTHED.value][:]

    # block height and width of the window/submatrix used in the cast
    # shadow algorithm
    block_width = margins.left + margins.right
    block_height = margins.top + margins.bottom

    # Compute the cast shadow mask
    ierr, mask = cast_shadow_main(elevation, zenith_angle, azimuth_angle,
                                  x_res, y_res, spheroid, y_origin, x_origin,
                                  margins.left, margins.right, margins.top,
                                  margins.bottom, block_height, block_width,
                                  is_utm)

    if ierr:
        raise CastShadowError(ierr)

    source_dir = 'SUN' if solar_source else 'SATELLITE'

    # Initialise the output file
    if out_group is None:
        fid = h5py.File('cast-shadow-{}.h5'.format(source_dir), driver='core',
                        backing_store=False)
    else:
        fid = out_group

    if GroupName.SHADOW_GROUP.value not in fid:
        fid.create_group(GroupName.SHADOW_GROUP.value)

    if filter_opts is None:
        filter_opts = {}
    else:
        filter_opts = filter_opts.copy()

    grp = fid[GroupName.SHADOW_GROUP.value]
    tile_size = satellite_solar_group[zenith_name].chunks
    filter_opts['chunks'] = tile_size
    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    kwargs['dtype'] = 'bool'

    dname_fmt = DatasetName.CAST_SHADOW_FMT.value
    out_dset = grp.create_dataset(dname_fmt.format(source=source_dir),
                                  data=mask, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal()}
    desc = ("The cast shadow mask determined using the {} "
            "as the source direction.").format(source_dir)
    attrs['description'] = desc
    attrs['alias'] = 'cast-shadow-{}'.format(source_dir).lower()
    attach_image_attributes(out_dset, attrs)

    if out_group is None:
        return fid


def _combine_shadow(self_shadow_fname, cast_shadow_sun_fname,
                    cast_shadow_satellite_fname, out_fname,
                    compression=H5CompressionFilter.LZF, filter_opts=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(self_shadow_fname, 'r') as fid_self,\
        h5py.File(cast_shadow_sun_fname, 'r') as fid_sun,\
        h5py.File(cast_shadow_satellite_fname, 'r') as fid_sat,\
        h5py.File(out_fname, 'w') as fid:

        grp1 = fid_self[GroupName.SHADOW_GROUP.value]
        grp2 = fid_sun[GroupName.SHADOW_GROUP.value]
        grp3 = fid_sat[GroupName.SHADOW_GROUP.value]
        combine_shadow_masks(grp1, grp2, grp3, fid, compression, filter_opts)

    link_shadow_datasets(self_shadow_fname, cast_shadow_sun_fname,
                         cast_shadow_satellite_fname, out_fname)


def combine_shadow_masks(self_shadow_group, cast_shadow_sun_group,
                         cast_shadow_satellite_group, out_group=None,
                         compression=H5CompressionFilter.LZF,
                         filter_opts=None):
    """
    A convienice function for combining the shadow masks into a single
    boolean array.

    :param self_shadow_group:
        The root HDF5 `Group` that contains the self shadow
        dataset specified by the pathname given by:

        * DatasetName.SELF_SHADOW

    :param cast_shadow_sun_group:
        The root HDF5 `Group` that contains the cast shadow
        (solar direction) dataset specified by the pathname
        given by:

        * DatasetName.CAST_SHADOW_FMT

    :param cast_shadow_sun_group:
        The root HDF5 `Group` that contains the cast shadow
        (satellite direction) dataset specified by the pathname
        given by:

        * DatasetName.CAST_SHDADOW_FMT

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be given by the format string detailed
        by:

        * DatasetName.COMBINED_SHADOW

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF 

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    # access the datasets
    dname_fmt = DatasetName.CAST_SHADOW_FMT.value
    self_shad = self_shadow_group[DatasetName.SELF_SHADOW.value]
    cast_sun = cast_shadow_sun_group[dname_fmt.format(source='SUN')]
    dname = dname_fmt.format(source='SATELLITE')
    cast_sat = cast_shadow_satellite_group[dname]
    geobox = GriddedGeoBox.from_dataset(self_shad)

    # Initialise the output files
    if out_group is None:
        fid = h5py.File('combined-shadow.h5', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    if GroupName.SHADOW_GROUP.value not in fid:
        fid.create_group(GroupName.SHADOW_GROUP.value)

    if filter_opts is None:
        filter_opts = {}
    else:
        filter_opts = filter_opts.copy()

    grp = fid[GroupName.SHADOW_GROUP.value]
    tile_size = cast_sun.chunks
    filter_opts['chunks'] = tile_size
    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    cols, rows = geobox.get_shape_xy()
    kwargs['shape'] = (rows, cols)
    kwargs['dtype'] = 'bool'

    # output dataset
    out_dset = grp.create_dataset(DatasetName.COMBINED_SHADOW.value, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal()}
    desc = ("Combined shadow masks: 1. self shadow, "
            "2. cast shadow (solar direction), "
            "3. cast shadow (satellite direction).")
    attrs['description'] = desc
    attrs['mask_values'] = "False = Shadow; True = Non Shadow"
    attrs['alias'] = 'terrain-shadow'
    attach_image_attributes(out_dset, attrs)

    # process by tile
    for tile in generate_tiles(cols, rows, tile_size[1], tile_size[0]):
        out_dset[tile] = (self_shad[tile] & cast_sun[tile] & cast_sat[tile])

    if out_group is None:
        return fid


def link_shadow_datasets(self_shadow_fname, cast_shadow_sun_fname,
                         cast_shadow_satellite_fname, out_fname):
    """
    Link the self shadow mask, and the two cast shadow masks into a
    single file for easier access.
    """
    group_path = GroupName.SHADOW_GROUP.value
    dname_fmt = DatasetName.CAST_SHADOW_FMT.value
    dname = ppjoin(group_path, DatasetName.SELF_SHADOW.value)
    create_external_link(self_shadow_fname, dname, out_fname, dname)

    dname = ppjoin(group_path, dname_fmt.format(source='SUN'))
    create_external_link(cast_shadow_sun_fname, dname, out_fname, dname)

    dname = ppjoin(group_path, dname_fmt.format(source='SATELLITE'))
    create_external_link(cast_shadow_satellite_fname, dname, out_fname, dname)
