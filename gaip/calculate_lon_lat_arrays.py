#!/usr/bin/env python

import numpy

from gaip import interpolate_grid
from gaip import write_img

def calculate_array(shape, eval_func, depth=7, dtype='float64'):
    """
    A general interface to interpolating co-ordinate arrays.

    :param shape:
        A tuple of 2 elements (y,x)  describing the array dimensions.

    :param eval_func:
        Co-ordinate location generation fuction.

    :param depth:
        The number of recursive bisections to be made. Default is 7.

    :param dtype:
        The datatype of the returned array. Default is float64.

    :return:
        A 2D NumPy array of type defined by the dtype keyword.
    """
    ndims = len(shape)
    if len(ndims) != 2:
        err = "Dimensions must be 2D; Received: {dims}".format(dims=ndims)
        raise ValueError(err)

    # Initialise the array to contain the result
    array = numpy.zeros(shape, dtype=dtype)

    interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
                     eval_func=eval_func, grid=array)

    return array

def create_lon_lat_grids(geobox, lon_fname='LON.tif', lat_fname='LAT.tif',
    work_dir='', to_disk=True):
    """
    Creates 2 by 2D NumPy arrays containing longitude and latitude
    co-ordinates for each array element.

    :param geobox:
        A GriddedGeoBox object.

    :param lon_fname:
        If the keyword to_disk is set to True (Default) then the
        longitude array will be written to disk rather than returned.
        The format of the output array is a GeoTiff.
        Default filename is LON.tif.

    :param lat_fname:
        If the keyword to_disk is set to True (Default) then the
        latitude array will be written to disk rather than returned.
        The format of the output file is a GeoTiff.
        Default filename is LAT.tif.

    :param work_dir:
        A string containing the full system filepath to a directory
        that will contain the longitude and latitude files. If unset
        then files will be written to the current working directory.

    :param to_disk:
        If set to True (Default), then the longitude and latitude
        arrays will be written to disk instead of being returned to
        the caller.

    :return:
        If to_disk is set to True (Default), then the longitude and
        latitude arrays are written to disk. If set to False, then
        the longitude and latitude arrays are returned as a tuple
        (longitude, latitude) 2D float64 NumPy arrays.
    """
    # Define the transform funtions
    lon_func = None
    lat_func = None

    lon_arr = calculate_array(shape, lon_func)

    crs = geobox.crs.ExportToWkt()
    transform = geobox.affine.to_gdal()

    if to_disk:
        lon_fname = os.path.join(work_dir, lon_fname)
        write_img(lon_arr, lon_fname, format='GTiff', geotransform=transform,
            projection=crs)
        lon_arr = None

    lat_arr = calculate_array(shape, lat_func)

    if to_disk:
        lat_fname = os.path.join(work_dir, lat_fname)
        write_img(lat_arr, lat_fname, format='GTiff', geotransform=transform,
            projection=crs)
        lat_arr = None
        return
    else:
        return (lon_arr, lat_arr)

