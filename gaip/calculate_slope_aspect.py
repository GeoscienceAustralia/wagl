#!/usr/bin/env

import gdal
import numpy
import rasterio

from EOtools import tiling
from gaip import ImageMargins
from gaip import setup_spheroid
from gaip import _slope_aspect
from gaip import write_img

X_TILE = None
Y_TILE = 20


def write_header_slope_file(file_name, margins, geobox):
    """Write the header slope file."""
    with open(file_name, 'w') as output:
        # get dimensions, resolution and pixel origin
        rows, cols = geobox.shape
        res = geobox.pixelsize
        origin = geobox.origin

        # Now output the details
        output.write("{nr} {nc}\n".format(nr=rows, nc=cols))
        output.write("{0} {1}\n{2} {3}\n".format(margins.left, margins.right,
                                                 margins.top, margins.bottom))
        output.write("{resy} {resx}\n".format(resx=res[0], resy=res[1]))
        output.write("{yorigin} {xorigin}\n".format(xorigin=origin[0],
                                                    yorigin=origin[1]))


def slope_aspect(acquisition, dsm_fname, margins, slope_out_fname,
                 aspect_out_fname, header_slope_fname=None):
    """
    Calculates slope and aspect.

    :param acquisition:
        An instance of an acquisition object.

    :param dsm_fname:
        A string containing the full file path name to the Digital
        Surface Model to be used in deriving the surface angles.

    :param margins:
        An object with members top, bottom, left and right giving the
        size of the margins (in pixels) which have been added to the
        corresponding sides of dsm.

    :param slope_out_fname:
        A string containing the full file path name to be used for
        writing the slope image on disk.

    :param aspect_out_fname:
        A string containing the full file path name to be used for
        writing the aspect image on disk.

    :param header_slope_fname:
        A string containing the full file path name to be used for
        writing the header slope text file to disk.

    :return:
        None. Outputs are written to disk.
    """

    # Setup the geobox
    geobox = acquisition.gridded_geo_box()

    # Retrive the spheroid parameters
    # (used in calculating pixel size in metres per lat/lon)
    spheroid = setup_spheroid(geobox.crs.ExportToWkt())

    # Are we in projected or geographic space
    is_utm = not geobox.crs.IsGeographic()

    # Define Top, Bottom, Left, Right pixel margins
    pixel_margin = ImageMargins(margins)

    # Get the x and y pixel sizes
    _, y_origin = geobox.origin
    x_res, y_res = geobox.pixelsize
    dresx = x_res + 2
    dresy = y_res + 2

    # Get acquisition dimensions and add 1 pixel top, bottom, left & right
    cols, rows = geobox.get_shape_xy()
    ncol = cols + 2
    nrow = rows + 2

    # Define the index to read the DEM subset
    idx = ((pixel_margin.top - 1, -(margin.bottom - 1)),
           (pixel_margin.left - 1, -(pixel_margin.right - 1)))

    with rasterio.open(dsm_fname) as dsm_ds:
        dsm_subset = as_array(dsm_ds.read_band(1, window=idx, masked=False),
                              dtype=numpy.float32, transpose=True)

    # Define an array of latitudes
    # This will be ignored if is_utm == True
    alat = numpy.array([y_origin - i * dresy for i in range(-1, nrow - 1)],
                       dtype=numpy.float64)  # yes, I did mean float64.

    # Define the output arrays. These will be transposed upon input
    slope = numpy.zeros((rows, cols), dtype='float32').transpose()
    aspect = numpy.zeros((rows, cols), dtype='float32').transpose()

    _slope_aspect(ncol, nrow, cols, rows, dresx, dresy, spheroid, alat, is_utm,
                  dsm_subset, slope.transpose(), aspect.transpose())

    # Output the results
    write_img(slope, slope_out_fname, geobox=geobox, nodata=-999)
    write_img(aspect, aspect_out_fname, geobox=geobox, nodata=-999)

    if header_slope_fname:
        write_header_slope_file(header_slope_fname, pixel_margin, geobox)
