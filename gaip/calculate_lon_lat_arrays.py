#!/usr/bin/env python

"""
Longitude and Latitude 2D grid creation.
"""

from __future__ import absolute_import
from functools import partial
import numpy
import osr
import h5py

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


def create_lon_lat_grids(acquisition, depth=7, dtype='float64',
                         out_fname=None, compression='lzf'):
    """
    Creates 2 by 2D NumPy arrays containing longitude and latitude
    co-ordinates for each array element.

    :param acquisition:
        An instance of an acquisitions object.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset path names will be as follows:

        * longitude
        * latitude

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
    # retrieve the geobox
    geobox = acquisition.gridded_geo_box()

    # Define the lon and lat transform funtions
    lon_func = partial(get_lon_coordinate, geobox=geobox, centre=True)
    lat_func = partial(get_lat_coordinate, geobox=geobox, centre=True)

    # Get some basic info about the image
    shape = geobox.get_shape_yx()

    # Initialise the array to contain the result
    result = numpy.zeros(shape, dtype=dtype)
    interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
                     eval_func=lon_func, grid=result)

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('longitude-latitude.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    # define some base attributes for the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.affine.to_gdal()}

    attrs['Description'] = LON_DESC
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))
    lon_dset = fid.create_dataset('longitude', data=result, **kwargs)
    attach_image_attributes(lon_dset, attrs)

    result = numpy.zeros(shape, dtype=dtype)
    interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
                     eval_func=lat_func, grid=result)

    attrs['Description'] = LAT_DESC
    lat_dset = fid.create_dataset('latitude', data=result, **kwargs)
    attach_image_attributes(lat_dset, attrs)

    fid.flush()
    return fid


def create_grid(geobox, coord_fn, depth=7, dtype='float64'):
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
    arr = numpy.zeros(shape, dtype=dtype)
    interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
                     eval_func=func, grid=arr)

    return arr


def create_lon_grid(acquisition, out_fname=None, compression='lzf', depth=7,
                    dtype='float64'):
    """Create longitude grid.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset path names will be as follows:

        * longitude

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
    # retrieve the geobox
    geobox = acquisition.gridded_geo_box()

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('longitude.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    # define some base attributes for the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.affine.to_gdal()}
    attrs['Description'] = LON_DESC
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))

    lon_grid = create_grid(geobox, get_lon_coordinate, depth, dtype)

    dset = fid.create_dataset('longitude', data=lon_grid, **kwargs)
    attach_image_attributes(dset, attrs)

    fid.flush()
    return fid


def create_lat_grid(acquisition, out_fname=None, compression='lzf', depth=7,
                    dtype='float64'):
    """Create latitude grid.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset path names will be as follows:

        * latitude

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
    # retrieve the geobox
    geobox = acquisition.gridded_geo_box()

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('latitude.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    # define some base attributes for the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.affine.to_gdal()}
    attrs['Description'] = LAT_DESC
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))

    lat_grid = create_grid(geobox, get_lat_coordinate, depth, dtype)

    dset = fid.create_dataset('latitude', data=lat_grid, **kwargs)
    attach_image_attributes(dset, attrs)

    fid.flush()
    return fid
