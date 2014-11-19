#!/usr/bin/env python

import rasterio
from rasterio import Affine

def write_img(array, filename, format='ENVI', geotransform=None, projection=None):
    """
    Writes a 2D/3D image to disk using rasterio.

    :param array:
        A 2D/3D NumPy array.

    :param filename:
        A string containing the output file name.

    :param format:
        A string containing a GDAL compliant image format. Default is
        'ENVI'.

    :param projection:
        A variable containing the projection information of the array.

    :param geotransform:
        A variable containing the geotransform information for the
        array.
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

    # Convert the GDAL geotransform to rasterio Affine transform
    transform = Affine.from_gdal(*geotransform)

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
        A 2D or 3D NumPy array containing the image subset.

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
        prj = str(src.crs_wkt) # rasterio returns a unicode

        # Get the original geotransform
        base_gt = src.get_transform()

        # Get the new UL co-ordinates of the array
        ULx, ULy = src.affine * (xstart, ystart)

        # Setup the new geotransform
        geot = (ULx, base_gt[1], base_gt[2], ULy, base_gt[4], base_gt[5])

    return (subs, geot, prj)


def read_img(fname):
    """
    A small and simple routine to read a GDAL compliant image.
    This is only intended for reading the raw file into a NumPy memory
    variable.
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

