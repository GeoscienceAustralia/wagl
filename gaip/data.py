"""
Data access functions
---------------------
"""

from __future__ import absolute_import
from __future__ import print_function
import subprocess
import numpy as np
import rasterio
import gaip
import os

from osgeo import gdal
from rasterio import crs
from rasterio.warp import reproject
from rasterio.warp import RESAMPLING
from os.path import join as pjoin


def get_pixel(filename, lonlat, band=1):
    """Return a pixel from `filename` at the longitude and latitude given
    by the tuple `lonlat`. Optionally, the `band` can be specified."""
    with rasterio.open(filename) as src:
        x, y = [int(v) for v in ~src.affine * lonlat]
        if isinstance(band, list):
            data = src.read(band, window=((y, y + 1), (x, x + 1))).ravel()
        else:
            data = src.read(band, window=((y, y + 1), (x, x + 1))).flat[0]
        return data


def no_data(acq):
    """
    Read the supplied acquisition's no data value and return to caller.
    The parameter `acq` should behave like a `gaip.Acquisition` object.
    """
    dirname = acq.dir_name
    filename = acq.file_name
    with rasterio.open(pjoin(dirname, filename), 'r') as fo:
        nodata_list = fo.nodatavals
        # currently we don't support multi-band GeoTiffs
        # idx = acq.band_num-1
        #
        idx = 0
        return nodata_list[idx]

def data(acq, out=None, window=None, masked=False, apply_gain_offset=False,
         out_no_data=-999):
    """
    Read the supplied acquisition's data into the `out` array if provided,
    otherwise return a new `numpy.array` containing the data.
    The parameter `acq` should behave like a `gaip.Acquisition` object.
    The parameter `window` defines a subset ((ystart, yend), (xstart, xend))
    in array co-ordinates. Default is None.
    The parameter `masked` indicates whether or not to return a masked array.
    Default is False.
    """
    dirname = acq.dir_name
    filename = acq.file_name
    with rasterio.open(pjoin(dirname, filename), 'r') as fo:
        # convert to at sensor radiance as required
        if apply_gain_offset:
            if out is None:
                out = fo.read(1, window=window, masked=masked)
            else:
                fo.read(1, out=out, window=window, masked=masked)

            # check for no data
            no_data = acq.no_data if acq.no_data is not None else 0
            nulls = out == no_data

            # gain & offset; y = mx + b
            data = acq.gain * out + acq.bias

            # set the out_no_data value inplace of the input no data value
            data[nulls] = out_no_data

            return data
        return fo.read(1, out=out, window=window, masked=masked)


def data_and_box(acq, out=None, window=None, masked=False):
    """
    Return a tuple comprising the `numpy.array` containing the data of
    the acquisition `acq` together with the associated GriddedGeoBox describing
    the data extent. 
    The parameter `acq` should behave like a `gaip.Acquisition` object.
    The `out` parameter, if supplied is a numpy.array into which the
    acquisition data is read.
    The parameter `window` defines a subset ((ystart, yend), (xstart, xend))
    in array co-ordinates. Default is None.
    The parameter `masked` indicates whether or not to return a masked array.
    Default is False.
    """
    dirname = acq.dir_name
    filename = acq.file_name
    with rasterio.open(pjoin(dirname, filename), 'r') as fo:
        box = gaip.GriddedGeoBox.from_dataset(fo)
        if window is not None:
            rows = window[0][1] - window[0][0]
            cols = window[1][1] - window[1][0]
            prj = bytes(fo.crs_wkt)
            res = fo.res
            # Get the new UL co-ordinates of the array
            ul_x, ul_y = fo.affine * (window[1][0], window[0][0])
            box = gaip.GriddedGeoBox(shape=(rows, cols), origin=(ul_x, ul_y),
                                     pixelsize=res, crs=prj)
        return (fo.read(1, out=out, window=window, masked=masked), box)


def gridded_geo_box(acq):
    """Return a GriddedGeoBox instance representing the spatial extent and 
    grid associated with the acquisition `acq`. 
    The parameter `acq` should behave like a `gaip.Acquisition` object."""
    dirname = acq.dir_name
    filename = acq.file_name
    with rasterio.open(pjoin(dirname, filename), 'r') as fo:
        return gaip.GriddedGeoBox.from_dataset(fo)


def select_acquisitions(acqs_list, fn=(lambda acq: True)):
    """
    Given a list of acquisitions, apply the supplied fn to select the
    desired acquisitions.
    """
    acqs = [acq for acq in acqs_list if fn(acq)]
    return acqs


def stack_data(acqs_list, fn=(lambda acq: True), window=None, masked=False):
    """
    Given a list of acquisitions, return the data from each acquisition
    collected in a 3D numpy array (first index is the acquisition number).
    If window is defined, then the subset contained within the window is
    returned along with a GriddedGeoBox instance detailing the 
    spatial information associated with that subset.

    :param acqs_list:
        The list of acquisitions from which to generate a stack of data.

    :param window:
        Defines a subset ((ystart, yend), (xstart, xend)) in array
        co-ordinates. Default is None.

    :param masked:
        Indicates whether or not to return a masked array. Default is False.

    :return:
        A 2-tuple containing:

            * 1. A 3D numpy array (or None) containing the corresponding
                 acquisition data. (None if no data).
            * 2. A GriddedGeoBox instance specifying the spatial context
                 of the 3D numpy array. Note: All Acquisitions share the
                 same GriddedGeoBox.
    """
    # determine data type and dimensions by reading the first band
    acqs = acqs_list
    a, geo_box = acqs[0].data_and_box(window=window, masked=masked)

    # create the result array, setting datatype based on source type

    stack_shape = (len(acqs), a.shape[0], a.shape[1])
    stack = np.empty(stack_shape, a.dtype)
    stack[0] = a
    del a

    # read remaining aquisitions into it

    for i in range(1, stack_shape[0]):
        # can't use this statement because it will cause data to be
        # resampled. But we want an exception thrown if the user
        # tries to stack irreqular aquisitions
        # acqs[i].data(out=stack[i])
        stack[i] = acqs[i].data(window=window, masked=masked)

    return stack, geo_box


def write_img(array, filename, fmt='ENVI', geobox=None, nodata=None,
              compress=None, tags=None, options=None):
    """
    Writes a 2D/3D image to disk using rasterio.

    :param array:
        A 2D/3D NumPy array.

    :param filename:
        A string containing the output file name.

    :param fmt:
        A string containing a GDAL compliant image fmt. Default is
        'ENVI'.

    :param geobox:
        An instance of a GriddedGeoBox object.

    :param nodata:
        A value representing the no data value for the array.

    :param compress:
        A compression algorithm name (e.g. 'lzw').

    :param tags:
        A dictionary of dataset-level metadata.

    :param options:
        A dictionary containing other dataset creation options.
        See creation options for the respective GDAL formats.
    """
    # Get the datatype of the array
    dtype = array.dtype.name

    # Check for excluded datatypes
    excluded_dtypes = ['int64', 'int8', 'uint64']
    if dtype in excluded_dtypes:
        msg = "Datatype not supported: {dt}".format(dt=dtype)
        raise TypeError(msg)

    ndims = array.ndim
    dims = array.shape

    # Get the (z, y, x) dimensions (assuming BSQ interleave)
    if ndims == 2:
        samples = dims[1]
        lines = dims[0]
        bands = 1
    elif ndims == 3:
        samples = dims[2]
        lines = dims[1]
        bands = dims[0]
    else:
        print('Input array is not of 2 or 3 dimensions!!!')
        err = 'Array dimensions: {dims}'.format(dims=ndims)
        raise IndexError(err)

    # If we have a geobox, then retrieve the geotransform and projection
    if geobox is not None:
        transform = geobox.affine
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
              'driver': fmt,
              'nodata': nodata}

    # compression predictor choices
    predictor = {'int8': 2,
                 'uint8': 2,
                 'int16': 2,
                 'uint16': 2,
                 'int32': 2,
                 'uint32': 2,
                 'int64': 2,
                 'uint64': 2,
                 'float32': 3,
                 'float64': 3}

    if fmt == 'GTiff' and compress is not None: 
        kwargs['compress'] = compress
        kwargs['predictor'] = predictor[dtype]

    if options is not None:
        for key in options:
            kwargs[key] = options[key]

    with rasterio.open(filename, 'w', **kwargs) as outds:
        if bands == 1:
            outds.write(array, 1)
        else:
            for i in range(bands):
                outds.write(array[i], i + 1)
        if tags is not None:
            outds.update_tags(**tags)


def read_subset(fname, ul_xy, ur_xy, lr_xy, ll_xy, bands=1):
    """
    Return a 2D or 3D NumPy array subsetted to the given bounding
    extents.

    :param fname:
        A string containing the full file pathname to an image on
        disk.

    :param ul_xy:
        A tuple containing the Upper Left (x,y) co-ordinate pair
        in real world (map) co-ordinates.  Co-ordinate pairs can be
        (longitude, latitude) or (eastings, northings), but they must
        be of the same reference as the image of interest.

    :param ur_xy:
        A tuple containing the Upper Right (x,y) co-ordinate pair
        in real world (map) co-ordinates.  Co-ordinate pairs can be
        (longitude, latitude) or (eastings, northings), but they must
        be of the same reference as the image of interest.

    :param lr_xy:
        A tuple containing the Lower Right (x,y) co-ordinate pair
        in real world (map) co-ordinates.  Co-ordinate pairs can be
        (longitude, latitude) or (eastings, northings), but they must
        be of the same reference as the image of interest.

    :param ll_xy:
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

            * 1. 2D or 3D NumPy array containing the image subset.
            * 2. A list of length 6 containing the GDAL geotransform.
            * 3. A WKT formatted string representing the co-ordinate
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

        # Get the dimensions
        cols = src.width
        rows = src.height

        # Convert each map co-ordinate to image/array co-ordinates
        img_ul_x, img_ul_y = [int(v) for v in inv * ul_xy]
        img_ur_x, img_ur_y = [int(v) for v in inv * ur_xy]
        img_lr_x, img_lr_y = [int(v) for v in inv * lr_xy]
        img_ll_x, img_ll_y = [int(v) for v in inv * ll_xy]

        # Calculate the min and max array extents
        # The ending array extents have +1 to account for Python's
        # [inclusive, exclusive) index notation.
        xstart = min(img_ul_x, img_ll_x)
        ystart = min(img_ul_y, img_ur_y)
        xend = max(img_ur_x, img_lr_x) + 1
        yend = max(img_ll_y, img_lr_y) + 1

        # Check for out of bounds
        if (((xstart < 0) or (ystart < 0)) or
            ((xend -1 > cols) or (yend -1 > rows))):
            msg = ("Error! Attempt to read a subset that is outside of the"
                   "image domain. Index: ({ys}, {ye}), ({xs}, {xe}))")
            msg = msg.format(ys=ystart, ye=yend, xs=xstart, xe=xend)
            raise IndexError(msg)

        # Read the subset
        subs = src.read(bands, window=((ystart, yend), (xstart, xend)))

        # Get the projection as WKT
        prj = bytes(src.crs_wkt)  # rasterio returns a unicode

        # Get the new UL co-ordinates of the array
        ul_x, ul_y = src.affine * (xstart, ystart)

        # Get the x & y pixel resolution
        res = src.res

        geobox = gaip.GriddedGeoBox(shape=subs.shape, origin=(ul_x, ul_y),
                                    pixelsize=res, crs=prj)

    return (subs, geobox)


def read_img(fname):
    """
    A small and simple routine to read a GDAL compliant image.
    This is only intended for reading the raw file into a NumPy memory
    variable.
    Largely used in the unittesting suite.

    :param fname:
        A string containing the full file pathname to an image on
        disk.

    :return:
        A NumPy array containing the full dimensions and specific
        datatype of the image represented on disk.
    """

    ds = gdal.Open(fname)

    img = ds.ReadAsArray()

    ds = None

    return img


def find_file(path, filename):
    """
    A simple routine for checking existance of files on disk.
    No error catching, it'll bail out of the main level program
    as it is designed for the unittests.

    :param path:
        A string containing the directory to search.

    :param filename:
        A string containing the name of the file to search within the
        given directory.

    :return:
        If the file exists, then a full filepath name for the file of
        interest will be returned. If the file is not found then an
        IOError will be raised.
    """
    fname = os.path.join(path, filename)
    if os.path.isfile(fname):
        return fname
    else:
        err = "Error! file not found: {filename}".format(filename=fname)
        raise IOError(err)


def reproject_file_to_array(src_filename, src_band=1, dst_geobox=None,
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

    if not isinstance(dst_geobox, gaip.GriddedGeoBox):
        msg = 'dst_geobox must be an instance of a GriddedGeoBox! Type: {}'
        msg = msg.format(type(dst_geobox))
        raise TypeError(msg)

    with rasterio.open(src_filename) as src:
        # Define a rasterio band
        rio_band = rasterio.band(src, src_band)

        # Define the output NumPy array
        dst_arr = np.zeros(dst_geobox.shape, dtype=src.dtypes[0])

        # Get the rasterio proj4 styled dict
        prj = crs.from_string(dst_geobox.crs.ExportToProj4())

        reproject(rio_band, dst_arr, dst_transform=dst_geobox.affine,
                  dst_crs=prj, resampling=resampling)

    return dst_arr


def reproject_img_to_img(src_img, src_geobox, dst_geobox,
                         resampling=RESAMPLING.nearest):
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

    if not isinstance(dst_geobox, gaip.GriddedGeoBox):
        msg = 'dst_geobox must be an instance of a GriddedGeoBox! Type: {}'
        msg = msg.format(type(dst_geobox))
        raise TypeError(msg)

    if not isinstance(src_geobox, gaip.GriddedGeoBox):
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
    dst_arr = np.zeros(dst_geobox.shape, dtype=src_img.dtype)

    reproject(src_img, dst_arr, src_transform=src_trans,
              src_crs=src_prj, dst_transform=dst_trans, dst_crs=dst_prj,
              resampling=resampling)

    return dst_arr


def as_array(array, dtype, transpose=False):
    """
    Given an array and dtype, array will be converted to dtype if
    and only if array.dtype != dtype. If transpose is set to True
    then array will be transposed before returning.

    :param array:
        A NumPy array.

    :param dtype:
        The type to return the array as.
    :type dtype:
        A NumPy data type (e.g. ``numpy.float32``).

    :param transpose:
        If set then array will be transposed before returning.
        Useful for passing arrays into Fortran routiines. Default is
        False.
    :type transpose:
        Bool.

    :return:
        A :py:class:`numpy.ndarry` of type ``dtype`` with the same
        dimensions as array.
    """
    if array.dtype != dtype:
        if transpose:
            return array.astype(dtype).transpose()
        else:
            return array.astype(dtype)
    else:
        if transpose:
            return array.transpose()
        else:
            return array


def read_meatadata_tags(fname, bands):
    """
    Retrieves the metadata tags for a list of bands from a `GDAL`
    compliant dataset.
    """
    with rasterio.open(fname) as ds:
        tag_data = {k: [] for k in ds.tags(1).keys()}
        for band in bands:
            tags = ds.tags(band)
            for tag in tags:
                tag_data[tag].append(tags[tag])

    return pandas.DataFrame(tag_data)
