#!/usr/bin/env python

from functools import partial
import os

import numpy
import osr

from EOtools.blrb import interpolate_grid
from gaip import write_img

CRS = "EPSG:4326"

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

def create_lon_lat_grids(geobox, depth=7, dtype='float64',
    lon_fname='LON.tif', lat_fname='LAT.tif', work_dir='', to_disk=True):
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

    # Get some basic info about the image
    crs = geobox.crs.ExportToWkt()
    transform = geobox.affine.to_gdal()
    shape = geobox.getShapeYX()

    # Initialise the array to contain the result
    lon_arr = numpy.zeros(shape, dtype=dtype)
    interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
        eval_func=lon_func, grid=lon_arr)

    if to_disk:
        lon_fname = os.path.join(work_dir, lon_fname)
        write_img(lon_arr, lon_fname, format='GTiff', geotransform=transform,
            projection=crs)
        lon_arr = None

    lat_arr = numpy.zeros(shape, dtype=dtype)
    interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
        eval_func=lat_func, grid=lat_arr)

    if to_disk:
        lat_fname = os.path.join(work_dir, lat_fname)
        write_img(lat_arr, lat_fname, format='GTiff', geotransform=transform,
            projection=crs)
        lat_arr = None
        return
    else:
        return (lon_arr, lat_arr)

