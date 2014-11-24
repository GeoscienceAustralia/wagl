#!/usr/bin/env python

from functools import partial

import numpy
import osr

from gaip import interpolate_grid
from gaip import write_img

GEOCRS = "WGS84"

def get_lon_coordinate(x, y, geobox, geo_crs=None, centre=False):
    """
    Given an image/array x & y co-ordinate return the corresponding
    longitude co-ordinate.

    :param x:
        An integer representing an image/array column coordinate.

    :param y:
        An integer representing an image/array row coordinate.

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
        sr = osr.SpatialReference()
        sr.SetWellKnownGeogCS(GEOCRS)

    xy = (x, y)
    mapx, mapy = geobox.convert_coordinates(xy, to_map=True, centre=centre)
    x = geobox.transform_coordinates((mapx, mapy), geo_crs)[0]

    return x

def get_lat_coordinate(x, y, geobox, geo_crs=None, centre=False):
    """
    Given an image/array x & y co-ordinate return the corresponding
    latitude co-ordinate.

    :param x:
        An integer representing an image/array column coordinate.

    :param y:
        An integer representing an image/array row coordinate.

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
        sr = osr.SpatialReference()
        sr.SetWellKnownGeogCS(GEOCRS)

    xy = (x, y)
    mapx, mapy = geobox.convert_coordinates(xy, to_map=True, centre=centre)
    y = geobox.transform_coordinates((mapx, mapy), geo_crs)[1]

    return y

def interpolate_array(shape, eval_func, depth=7, dtype='float64'):
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

    # Define the lon and lat transform funtions
    lon_func = partial(get_lon_coordinate, geobox=geobox, centre=True)
    lat_func = partial(get_lat_coordinate, geobox=geobox, centre=True)

    crs = geobox.crs.ExportToWkt()
    transform = geobox.affine.to_gdal()
    shape = geobox.getShapeYX()

    lon_arr = interpolate_array(shape, lon_func)

    if to_disk:
        lon_fname = os.path.join(work_dir, lon_fname)
        write_img(lon_arr, lon_fname, format='GTiff', geotransform=transform,
            projection=crs)
        lon_arr = None

    lat_arr = interpolate_array(shape, lat_func)

    if to_disk:
        lat_fname = os.path.join(work_dir, lat_fname)
        write_img(lat_arr, lat_fname, format='GTiff', geotransform=transform,
            projection=crs)
        lat_arr = None
        return
    else:
        return (lon_arr, lat_arr)

