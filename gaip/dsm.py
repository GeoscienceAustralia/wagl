#!/usr/bin/env python
"""
Digital Surface Model Data extraction and smoothing.
"""

from __future__ import absolute_import
import numpy
from scipy import ndimage
import h5py
from rasterio.warp import RESAMPLING
from gaip import ImageMargins
from gaip import GriddedGeoBox
from gaip import reproject_file_to_array
from gaip import dataset_compression_kwargs
from gaip import attach_image_attributes


def filter_dsm(array):
    """
    Applies a gaussian filter to array.

    :param array:
        A 2D NumPy array.

    :return:
        A 2D NumPy array.
    """
    # Define the kernel
    kernel = [0.009511, 0.078501, 0.009511, 0.078501, 0.647954, 0.078501,
              0.009511, 0.078501, 0.009511]
    kernel = numpy.array(kernel).reshape((3, 3))

    filtered = ndimage.convolve(array, kernel)
    return filtered


def get_dsm(acquisition, national_dsm, margins, out_fname=None,
            compression='lzf'):
    """
    Given an acquisition and a national Digitial Surface Model,
    extract a subset from the DSM based on the acquisition extents
    plus an x & y margins. The subset is then smoothed with a 3x3
    gaussian filter.
    A square margins is applied to the extents.

    :param acquisition:
        An instance of an acquisition object.

    :param national_dsm:
        A string containing the full filepath name to an image on
        disk containing national digital surface model.

    :param margin:
        An integer indictaing the number of pixels to be used as a
        margin around the aqcuisition.
        Eg, a value of 250 indicates that 250 pixels to the top,
        bottom, left and right will be added to the acquisition
        margin/border.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset names will be as follows:

        * dsm
        * dsm-smoothed

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
    # Use the 1st acquisition to setup the geobox
    geobox = acquisition.gridded_geo_box()
    shape = geobox.get_shape_yx()

    # Define Top, Bottom, Left, Right pixel margins
    pixel_buf = ImageMargins(margins)

    # Get the dimensions and geobox of the new image
    dem_cols = shape[1] + pixel_buf.left + pixel_buf.right
    dem_rows = shape[0] + pixel_buf.top + pixel_buf.bottom
    dem_shape = (dem_rows, dem_cols)
    dem_origin = geobox.convert_coordinates((0 - pixel_buf.left,
                                             0 - pixel_buf.top))
    dem_geobox = GriddedGeoBox(dem_shape, origin=dem_origin,
                               pixelsize=geobox.pixelsize,
                               crs=geobox.crs.ExportToWkt())

    # Retrive the DSM data
    dsm_data = reproject_file_to_array(national_dsm, dst_geobox=dem_geobox,
                                       resampling=RESAMPLING.bilinear)

    # Output the reprojected result
    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('dsm-extraction.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))

    grp = fid.create_group('parameters')
    grp.attrs['left_buffer'] = pixel_buf.left
    grp.attrs['right_buffer'] = pixel_buf.right
    grp.attrs['top_buffer'] = pixel_buf.top
    grp.attrs['bottom_buffer'] = pixel_buf.bottom

    # output datasets
    out_dset = fid.create_dataset('dsm', data=dsm_data, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': dem_geobox.affine.to_gdal()}
    desc = "A subset of a Digital Surface Model."
    attrs['Description'] = desc
    attach_image_attributes(out_dset, attrs)

    # Smooth the DSM
    dsm_data = filter_dsm(dsm_data)
    out_sm_dset = fid.create_dataset('dsm-smoothed', data=dsm_data, **kwargs)
    desc = ("A subset of a Digital Surface Model smoothed with a gaussian "
            "kernel.")
    attrs['Description'] = desc
    attach_image_attributes(out_sm_dset, attrs)

    fid.flush()
    return fid
