#!/usr/bin/env python

"""
Longitude and Latitude 2D grid creation.
"""

from __future__ import absolute_import, print_function
from functools import partial
import numpy
import osr
import h5py

from gaip.constants import DatasetName, GroupName
from gaip.interpolation import interpolate_grid
from gaip.hdf5 import dataset_compression_kwargs
from gaip.hdf5 import attach_image_attributes

CRS = "EPSG:4326"
LON_DESC = "Contains the longitude values for each pixel."
LAT_DESC = "Contains the latitude values for each pixel."


def get_lon_coordinate(y, x, geobox, geo_crs=None, centre=False):
    """
    Given an image/array y & x co-ordinate return the corresponding
    longitude co-ordinate. The y, x style mimics Python indices.

    :param y:
        An integer representing an image/array row coordinate.

    :param x:
        An integer representing an image/array column coordinate.

    :param geobox:
        An instance of a GriddedGeoBox object.

    :param geo_crs:
        An instance of a defined geographic osr.SpatialReference
        object. If set to None (Default), then geo_crs will be set
        to WGS84.

    :param centre:
        A boolean indicating whether or not the returned co-ordinate
        should reference the centre of a pixel, in which case a 0.5
        offset is applied in the x & y directions. Default is False.

    :return:
        A floating point value representing the longitude
        co-ordinate.
    """
    if geo_crs is None:
        geo_crs = osr.SpatialReference()
        geo_crs.SetFromUserInput(CRS)

    xy = (x, y)
    mapx, mapy = geobox.convert_coordinates(xy, to_map=True, centre=centre)
    x = geobox.transform_coordinates((mapx, mapy), geo_crs)[0]

    return x


def get_lat_coordinate(y, x, geobox, geo_crs=None, centre=False):
    """
    Given an image/array x & y co-ordinate return the corresponding
    latitude co-ordinate. The y, x style mimics Python indices.

    :param y:
        An integer representing an image/array row coordinate.

    :param x:
        An integer representing an image/array column coordinate.

    :param geobox:
        An instance of a GriddedGeoBox object.

    :param geo_crs:
        An instance of a defined geographic osr.SpatialReference
        object. If set to None (Default), then geo_crs will be set
        to WGS84.

    :param centre:
        A boolean indicating whether or not the returned co-ordinate
        should reference the centre of a pixel, in which case a 0.5
        offset is applied in the x & y directions. Default is False.

    :return:
        A floating point value representing the latitude
        co-ordinate.
    """
    if geo_crs is None:
        geo_crs = osr.SpatialReference()
        geo_crs.SetFromUserInput(CRS)

    xy = (x, y)
    mapx, mapy = geobox.convert_coordinates(xy, to_map=True, centre=centre)
    y = geobox.transform_coordinates((mapx, mapy), geo_crs)[1]

    return y


def _create_lon_lat_grids(acquisition, out_fname=None, compression='lzf',
                          depth=7):
    """
    A private wrapper for dealing with the internal custom workings of the
    multifile workflow.
    """
    with h5py.File(out_fname, 'w') as fid:
        create_lon_lat_grids(acquisition, fid, compression, depth)


def create_lon_lat_grids(acquisition, out_group=None, compression='lzf',
                         depth=7):
    """
    Creates 2 by 2D NumPy arrays containing longitude and latitude
    co-ordinates for each array element.

    :param acquisition:
        An instance of an `Acquisition` object.

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be given by:

        * contants.DatasetName.lon.value
        * contants.DatasetName.lat.value

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    geobox = acquisition.gridded_geo_box()
    # Define the lon and lat transform funtions
    lon_func = partial(get_lon_coordinate, geobox=geobox, centre=True)
    lat_func = partial(get_lat_coordinate, geobox=geobox, centre=True)

    # Get some basic info about the image
    shape = geobox.get_shape_yx()

    # Initialise the array to contain the result
    result = numpy.zeros(shape, dtype='float64')
    interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
                     eval_func=lon_func, grid=result)

    # Initialise the output files
    if out_group is None:
        fid = h5py.File('longitude-latitude.h5', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    if GroupName.lon_lat_group.value not in fid:
        fid.create_group(GroupName.lon_lat_group.value)

    grp = fid[GroupName.lon_lat_group.value]

    # define some base attributes for the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal()}

    attrs['description'] = LON_DESC
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=acquisition.tile_size)
    lon_dset = grp.create_dataset(DatasetName.lon.value, data=result, **kwargs)
    attach_image_attributes(lon_dset, attrs)

    result = numpy.zeros(shape, dtype='float64')
    interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
                     eval_func=lat_func, grid=result)

    attrs['description'] = LAT_DESC
    lat_dset = grp.create_dataset(DatasetName.lat.value, data=result, **kwargs)
    attach_image_attributes(lat_dset, attrs)

    return fid


def create_grid(geobox, coord_fn, depth=7):
    """
    Interpolates a `NumPy` array based on the input coordinate function
    `coord_fn`.

    :param geobox:
        An instance of an `GriddedGeoBox` object.

    :param coord_fn:
        A function that maps coordinates.

    :return:
        A `NumPy` array.
    """
    # Define the transform funtions
    func = partial(coord_fn, geobox=geobox, centre=True)

    # Get some basic info about the image
    shape = geobox.get_shape_yx()

    # Initialise the array to contain the result
    arr = numpy.zeros(shape, dtype='float64')
    interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
                     eval_func=func, grid=arr)

    return arr


def create_lon_grid(acquisition, out_fname=None, compression='lzf', depth=7):
    """Create longitude grid.

    :param acquisition:
        An instance of an `Acquisition` object.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset path names will be as follows:

        * contants.DatasetName.lon.value

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('longitude.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    geobox = acquisition.gridded_geo_box()

    # define some base attributes for the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal()}
    attrs['description'] = LON_DESC
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=acquisition.tile_size)

    lon_grid = create_grid(geobox, get_lon_coordinate, depth)

    grp = fid.create_group(GroupName.lon_lat_group.value)
    dset = grp.create_dataset(DatasetName.lon.value, data=lon_grid, **kwargs)
    attach_image_attributes(dset, attrs)

    return fid


def create_lat_grid(acquisition, out_fname=None, compression='lzf', depth=7):
    """Create latitude grid.

    :param acquisition:
        An instance of an `Acquisition` object.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset path names will be as follows:

        * contants.DatasetName.lat.value

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('latitude.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    geobox = acquisition.gridded_geo_box()

    # define some base attributes for the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal()}
    attrs['description'] = LAT_DESC
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=acquisition.tile_size)

    lat_grid = create_grid(geobox, get_lat_coordinate, depth)

    grp = fid.create_group(GroupName.lon_lat_group.value)
    dset = grp.create_dataset(DatasetName.lat.value, data=lat_grid, **kwargs)
    attach_image_attributes(dset, attrs)

    return fid
