#!/usr/bin/env python

import numpy
import rasterio
from rasterio import Affine
from rasterio import crs
from rasterio.warp import reproject
from rasterio.warp import RESAMPLING

from GriddedGeoBox import GriddedGeoBox

def write_img(array, filename, format='ENVI', geobox=None):
    """
    Writes a 2D/3D image to disk using rasterio.

    :param array:
        A 2D/3D NumPy array.

    :param filename:
        A string containing the output file name.

    :param format:
        A string containing a GDAL compliant image format. Default is
        'ENVI'.

    :param geobox:
        An instance of a GriddedGeoBox object.
    """

    dtype = array.dtype.name
    # If we have an excluded datatype default to float64 which should cover
    # most data values
    excluded_dtypes = ['int64', 'int8', 'uint64']
    if dtype in excluded_dtypes:
        dtype = 'float64'

    ndims = array.ndim
    dims  = array.shape

    # Get the (z, y, x) dimensions (assuming BSQ interleave)
    if ndims == 2:
        samples = dims[1]
        lines   = dims[0]
        bands   = 1
    elif ndims == 3:
        samples = dims[2]
        lines   = dims[1]
        bands   = dims[0]
    else:
        print 'Input array is not of 2 or 3 dimensions!!!'
        err = 'Array dimensions: {dims}'.format(dims=ndims)
        raise IndexError(err)

    # If we have a geobox, then retrieve the geotransform and projection
    if geobox is not None:
        transform  = geobox.affine
        projection = bytes(geobox.crs.ExportToWkt())
    else:
        transform = None
        projection = None

    kwargs = {'count': bands,
              'width': samples,
              'height': lines,
              'crs': projection,
              'transform': transform,
              'dtype': dtype,
              'driver': format
             }

    with rasterio.open(filename, 'w', **kwargs) as outds:
        if bands == 1:
            outds.write_band(1, array)
        else:
            for i in range(bands):
                outds.write_band(i+1, array[i])


def read_subset(fname, ULxy, URxy, LRxy, LLxy, bands=1):
    """
    Return a 2D or 3D NumPy array subsetted to the given bounding
    extents.

    :param fname:
        A string containing the full file pathname to an image on
        disk.

    :param ULxy:
        A tuple containing the Upper Left (x,y) co-ordinate pair
        in real world (map) co-ordinates.  Co-ordinate pairs can be
        (longitude, latitude) or (eastings, northings), but they must
        be of the same reference as the image of interest.

    :param URxy:
        A tuple containing the Upper Right (x,y) co-ordinate pair
        in real world (map) co-ordinates.  Co-ordinate pairs can be
        (longitude, latitude) or (eastings, northings), but they must
        be of the same reference as the image of interest.

    :param LRxy:
        A tuple containing the Lower Right (x,y) co-ordinate pair
        in real world (map) co-ordinates.  Co-ordinate pairs can be
        (longitude, latitude) or (eastings, northings), but they must
        be of the same reference as the image of interest.

    :param LLxy:
        A tuple containing the Lower Left (x,y) co-ordinate pair
        in real world (map) co-ordinates.  Co-ordinate pairs can be
        (longitude, latitude) or (eastings, northings), but they must
        be of the same reference as the image of interest.

    :param bands:
        Can be an integer of list of integers representing the band(s)
        to be read from disk.  If bands is a list, then the returned
        subset will be 3D, otherwise the subset will be strictly 2D.

    :return:
        A tuple of 3 elements:
            [0] 2D or 3D NumPy array containing the image subset.
            [1] A list of length 6 containing the GDAL geotransform.
            [2] A WKT formatted string representing the co-ordinate
                reference system (projection).

    :additional notes:
        The ending array co-ordinates are increased by +1,
        i.e. xend = 270 + 1
        to account for Python's [inclusive, exclusive) index notation.
    """

    # Open the file
    with rasterio.open(fname) as src:

        # Get the inverse transform of the affine co-ordinate reference
        inv = ~src.affine

        # Convert each map co-ordinate to image/array co-ordinates
        imgULx, imgULy = [int(v) for v in inv*ULxy]
        imgURx, imgURy = [int(v) for v in inv*URxy]
        imgLRx, imgLRy = [int(v) for v in inv*LRxy]
        imgLLx, imgLLy = [int(v) for v in inv*LLxy]

        # Calculate the min and max array extents
        # The ending array extents have +1 to account for Python's
        # [inclusive, exclusive) index notation.
        xstart = min(imgULx, imgLLx)
        ystart = min(imgULy, imgURy)
        xend = max(imgURx, imgLRx) + 1
        yend = max(imgLLy, imgLRy) + 1

        # Read the subset
        subs =  src.read(bands, window=((ystart, yend), (xstart, xend)))

        # Get the projection as WKT
        prj = bytes(src.crs_wkt) # rasterio returns a unicode

        # Get the new UL co-ordinates of the array
        ULx, ULy = src.affine * (xstart, ystart)

        # Get the x & y pixel resolution
        res = src.res

        geobox = GriddedGeoBox(shape=subs.shape, origin=(ULx, ULy),
            pixelsize=res, crs=prj)

    return (subs, geobox)


def read_img(fname):
    """
    A small and simple routine to read a GDAL compliant image.
    This is only intended for reading the raw file into a NumPy memory
    variable.
    Largely used in the unittesting suite.
    """

    ds = gdal.Open(fname)

    img = ds.ReadAsArray()

    ds = None

    return img


def find_file(dir, file):
    """
    A simple routine for checking existance of files on disk.
    No error catching, it'll bail out of the main level program
    as it is designed for the unittests.
    """
    fname = os.path.join(dir, file)
    if os.path.isfile(fname):
        return fname
    else:
        err = "Error! file not found: {filename}".format(filename=fname)
        raise IOError(err)


def reprojectFile2Array(src_filename, src_band=1, dst_geobox=None,
        resampling=RESAMPLING.nearest):
    """
    Given an image on file, reproject to the desired coordinate
    reference system.

    :param src_filename:
        A string containing the full file path name to the source
        image on disk.

    :param src_band:
        An integer representing the band number to be reprojected.
        Default is 1, the 1st band.

    :param dst_geobox:
        An instance of a GriddedGeoBox object containing the
        destination parameters such as origin, affine, projection,
        and array dimensions.

    :param resampling:
        An integer representing the resampling method to be used.
        check rasterio.warp.RESMPLING for more details.
        Default is 0, nearest neighbour resampling.

    :return:
        A NumPy array containing the reprojected result.
    """

    if not isinstance(dst_geobox, GriddedGeoBox):
        msg = 'dst_geobox must be an instance of a GriddedGeoBox! Type: {}'
        msg = msg.format(type(dst_geobox))
        raise TypeError(msg)

    with rasterio.open(src_filename) as src:
        # Define a rasterio band
        rio_band = rasterio.band(src, src_band)

        # Define the output NumPy array
        dst_arr = numpy.zeros(dst_geobox.shape, dtype=src.dtypes[0])

        # Get the rasterio proj4 styled dict
        prj = crs.from_string(dst_geobox.crs.ExportToProj4())

        reproject(rio_band, dst_arr, dst_transform=dst_geobox.affine,
            dst_crs=prj, resampling=resampling)

    return dst_arr


def reprojectImg2Img(src_img, src_geobox, dst_geobox, resampling=RESAMPLING.nearest):
    """
    Reprojects an image/array to the desired co-ordinate reference system.

    :param src_img:
        A NumPy array containing the source image.

    :param src_geobox:
        An instance of a GriddedGeoBox object containing the
        source parameters such as origin, affine, projection.

    :param dst_geobox:
        An instance of a GriddedGeoBox object containing the
        destination parameters such as origin, affine, projection,
        and array dimensions.

    :param resampling:
        An integer representing the resampling method to be used.
        check rasterio.warp.RESMPLING for more details.
        Default is 0, nearest neighbour resampling.

    :return:
        A NumPy array containing the reprojected result.
    """

    if not isinstance(dst_geobox, GriddedGeoBox):
        msg = 'dst_geobox must be an instance of a GriddedGeoBox! Type: {}'
        msg = msg.format(type(dst_geobox))
        raise TypeError(msg)

    if not isinstance(src_geobox, GriddedGeoBox):
        msg = 'src_geobox must be an instance of a GriddedGeoBox! Type: {}'
        msg = msg.format(type(src_geobox))
        raise TypeError(msg)

    # Get the source and destination projections in Proj4 styled dicts
    src_prj = crs.from_string(src_geobox.crs.ExportToProj4())
    dst_prj = crs.from_string(dst_geobox.crs.ExportToProj4())

    # Get the source and destination transforms
    src_trans = src_geobox.affine
    dst_trans = dst_geobox.affine

    # Define the output NumPy array
    dst_arr = numpy.zeros(dst_geobox.shape, dtype=scr_img.dtype)

    reproject(src_img, dst_arr, src_transform=src_trans,
        src_crs=src_prj, dst_transform=dst_trans, dst_crs=dst_prj,
        resampling=resampling)

    return dst_arr
