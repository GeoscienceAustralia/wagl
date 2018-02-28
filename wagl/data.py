"""
Data access functions
---------------------
"""

from __future__ import absolute_import
from os.path import join as pjoin, basename, dirname
import subprocess
import tempfile
import logging
import numpy as np
import h5py
import rasterio

from rasterio.crs import CRS
from rasterio.warp import reproject
from rasterio.enums import Resampling
from wagl.geobox import GriddedGeoBox
from wagl.tiling import generate_tiles


def get_pixel(filename, lonlat, band=1):
    """Return a pixel from `filename` at the longitude and latitude given
    by the tuple `lonlat`. Optionally, the `band` can be specified."""
    with rasterio.open(filename) as src:
        x, y = [int(v) for v in ~src.transform * lonlat]
        if isinstance(band, list):
            data = src.read(band, window=((y, y + 1), (x, x + 1))).ravel()
        else:
            data = src.read(band, window=((y, y + 1), (x, x + 1))).flat[0]
        return data


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


def write_img(array, filename, driver='GTiff', geobox=None, nodata=None,
              tags=None, options=None, cogtif=False, levels=None,
              resampling=Resampling.nearest):
    """
    Writes a 2D/3D image to disk using rasterio.

    :param array:
        A 2D/3D NumPy array.

    :param filename:
        A string containing the output file name.

    :param driver:
        A string containing a GDAL compliant image driver. Default is
        'GTiff'.

    :param geobox:
        An instance of a GriddedGeoBox object.

    :param nodata:
        A value representing the no data value for the array.

    :param tags:
        A dictionary of dataset-level metadata.

    :param options:
        A dictionary containing other dataset creation options.
        See creation options for the respective GDAL formats.

    :param cogtif:
        If set to True, override the `driver` keyword with `GTiff`
        and create a Cloud Optimised GeoTiff. Default is False.
        See:
        https://trac.osgeo.org/gdal/wiki/CloudOptimizedGeoTIFF

    :param levels:
        If cogtif is set to True, build overviews/pyramids
        according to levels. Default levels are [2, 4, 8, 16, 32].

    :param resampling:
        If cogtif is set to True, build overviews/pyramids using
        a resampling method from `rasterio.enums.Resampling`.
        Default is `Resampling.nearest`.

    :notes:
        If array is an instance of a `h5py.Dataset`, then the output
        file will include blocksizes based on the `h5py.Dataset's`
        chunks. To override the blocksizes, specify them using the
        `options` keyword. Eg {'blockxsize': 512, 'blockysize': 512}.
        If `cogtif` is set to True, the default blocksizes will be
        256x256. To override this behaviour, specify them using the
        `options` keyword. Eg {'blockxsize': 512, 'blockysize': 512}.
    """
    # Get the datatype of the array
    dtype = array.dtype.name

    # Check for excluded datatypes
    excluded_dtypes = ['int64', 'int8', 'uint64']
    if dtype in excluded_dtypes:
        msg = "Datatype not supported: {dt}".format(dt=dtype)
        raise TypeError(msg)

    # convert any bools to uin8
    if dtype == 'bool':
        array = np.uint8(array)
        dtype = 'uint8'

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
        logging.error('Input array is not of 2 or 3 dimensions!!!')
        err = 'Array dimensions: {dims}'.format(dims=ndims)
        raise IndexError(err)

    # If we have a geobox, then retrieve the geotransform and projection
    if geobox is not None:
        transform = geobox.transform
        projection = geobox.crs.ExportToWkt()
    else:
        transform = None
        projection = None

    # override the driver if we are creating a cogtif
    if cogtif:
        driver = 'GTiff'

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

    kwargs = {'count': bands,
              'width': samples,
              'height': lines,
              'crs': projection,
              'transform': transform,
              'dtype': dtype,
              'driver': driver,
              'nodata': nodata,
              'predictor': predictor[dtype]}

    if isinstance(array, h5py.Dataset):
        # TODO: if array is 3D get x & y chunks
        if array.chunks[1] == array.shape[1]:
            # GDAL doesn't like tiled or blocksize options to be set
            # the same length as the columns (probably true for rows as well)
            array = array[:]
        else:
            y_tile, x_tile = array.chunks
            tiles = generate_tiles(samples, lines, x_tile, y_tile)

            # add blocksizes to the creation keywords
            kwargs['tiled'] = 'yes'
            kwargs['blockxsize'] = x_tile
            kwargs['blockysize'] = y_tile

    # the user can override any derived blocksizes by supplying `options`
    if options is not None:
        for key in options:
            kwargs[key] = options[key]

    with tempfile.TemporaryDirectory() as tmpdir:
        out_fname = pjoin(tmpdir, basename(filename)) if cogtif else filename

        with rasterio.open(out_fname, 'w', **kwargs) as outds:
            if bands == 1:
                if isinstance(array, h5py.Dataset):
                    for tile in tiles:
                        idx = (slice(tile[0][0], tile[0][1]),
                               slice(tile[1][0], tile[1][1]))
                        outds.write(array[idx], 1, window=tile)
                else:
                    outds.write(array, 1)
            else:
                if isinstance(array, h5py.Dataset):
                    for tile in tiles:
                        idx = (slice(tile[0][0], tile[0][1]),
                               slice(tile[1][0], tile[1][1]))
                        subs = array[:, idx[0], idx[1]]
                        for i in range(bands):
                            outds.write(subs[i], i + 1, window=tile)
                else:
                    for i in range(bands):
                        outds.write(array[i], i + 1)
            if tags is not None:
                outds.update_tags(**tags)

            # overviews/pyramids
            if cogtif:
                if levels is None:
                    levels = [2, 4, 8, 16, 32]
                outds.build_overviews(levels, resampling)

        if cogtif:
            cmd = ['gdal_translate',
                   '-co',
                   'TILED=YES',
                   '-co',
                   'COPY_SRC_OVERVIEWS=YES',
                   '-co',
                   '{}={}'.format('PREDICTOR', predictor[dtype])]

            for key, value in options.items():
                cmd.extend(['-co', '{}={}'.format(key, value)])

            cmd.extend([out_fname, filename])

            subprocess.check_call(cmd, cwd=dirname(filename))


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
    if isinstance(fname, h5py.Dataset):
        geobox = GriddedGeoBox.from_dataset(fname)
        prj = fname.attrs['crs_wkt']
    else:
        # Open the file
        with rasterio.open(fname) as src:

            # Get the inverse transform of the affine co-ordinate reference
            geobox = GriddedGeoBox.from_dataset(src)
            prj = src.crs.wkt  # rasterio returns a unicode

    inv = ~geobox.transform
    rows, cols = geobox.shape

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

    if isinstance(fname, h5py.Dataset):
        subs = fname[ystart:yend, xstart:xend]
    else:
        with rasterio.open(fname) as src:
            subs = src.read(bands, window=((ystart, yend), (xstart, xend)))

    # Get the new UL co-ordinates of the array
    ul_x, ul_y = geobox.transform * (xstart, ystart)

    geobox_subs = GriddedGeoBox(shape=subs.shape, origin=(ul_x, ul_y),
                                pixelsize=geobox.pixelsize, crs=prj)

    return (subs, geobox_subs)


def reproject_file_to_array(src_filename, src_band=1, dst_geobox=None,
                            resampling=Resampling.nearest):
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
        dst_arr = np.zeros(dst_geobox.shape, dtype=src.dtypes[0])

        # Get the rasterio proj4 styled dict
        prj = CRS.from_string(dst_geobox.crs.ExportToProj4())

        reproject(rio_band, dst_arr, dst_transform=dst_geobox.transform,
                  dst_crs=prj, resampling=resampling)

    return dst_arr


def reproject_img_to_img(src_img, src_geobox, dst_geobox,
                         resampling=Resampling.nearest):
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
    src_prj = CRS.from_string(src_geobox.crs.ExportToProj4())
    dst_prj = CRS.from_string(dst_geobox.crs.ExportToProj4())

    # Get the source and destination transforms
    src_trans = src_geobox.transform
    dst_trans = dst_geobox.transform

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
        return array.astype(dtype)
    if transpose:
        return array.transpose()
    return array
