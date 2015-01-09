import numpy
import osr
import os

from functools import partial
from EOtools.blrb import interpolate_grid
from gaip import gridded_geo_box
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


def create_lon_lat_grids(acquisition, depth=7, dtype='float64',
                         lon_fname='LON.bin', lat_fname='LAT.bin',
                         work_dir='', to_disk=True):
    """
    Creates 2 by 2D NumPy arrays containing longitude and latitude
    co-ordinates for each array element.

    :param acquisition:
        An instance of an acquisitions object.

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

    # Compute the geobox
    geobox = gridded_geo_box(acquisition)

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
        write_img(lon_arr, lon_fname, geobox=geobox)
        lon_arr = None

    lat_arr = numpy.zeros(shape, dtype=dtype)
    interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
                     eval_func=lat_func, grid=lat_arr)

    if to_disk:
        lat_fname = os.path.join(work_dir, lat_fname)
        write_img(lat_arr, lat_fname, geobox=geobox)
        lat_arr = None
        return
    else:
        return (lon_arr, lat_arr)


def create_grid(acquisition, coord_fn, fname=None, depth=7, dtype='float64'):
    """
    Creates 2 by 2D NumPy arrays containing co-ordinates for each array element.

    :param acquisition:
        An instance of an acquisitions object.

    :param coord_fn:
        A function that maps coordinates.

    :param fname:
        If set then the array will be written to disk rather than returned.
        The format of the output array is a GeoTiff.

    :return:
        If fname is set, then the longitude and latitude arrays are written
        to disk. If set to False, then the array is returned as a tuple
        (longitude, latitude) 2D float64 NumPy arrays.
    """

    # Compute the geobox
    geobox = gridded_geo_box(acquisition)

    # Define the transform funtions
    func = partial(coord_fn, geobox=geobox, centre=True)

    # Get some basic info about the image
    crs = geobox.crs.ExportToWkt()
    transform = geobox.affine.to_gdal()
    shape = geobox.getShapeYX()

    # Initialise the array to contain the result
    arr = numpy.zeros(shape, dtype=dtype)
    interpolate_grid(depth=depth, origin=(0, 0), shape=shape,
                     eval_func=func, grid=arr)

    if fname is not None:
        write_img(arr, fname, geobox=geobox)
    else:
        return arr


def create_lon_grid(acquisition, fname=None, depth=7, dtype='float64'):
    return create_grid(acquisition, get_lon_coordinate, fname, depth, dtype)

def create_lat_grid(acquisition, fname=None, depth=7, dtype='float64'):
    return create_grid(acquisition, get_lat_coordinate, fname, depth, dtype)
