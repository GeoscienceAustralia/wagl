"""
Shadow Calculations
-------------------
"""

import gdal
import numpy
import rasterio

from gaip import ImageMargins
from gaip import setup_spheroid
from gaip import read_img
from gaip import run_slope


X_TILE = None
Y_TILE = 20

def calculate_self_shadow(acquisition, dsm_fname, margins,
                          solar_zenith_fname, solar_azimuth_fname,
                          satellite_view_fname, satellite_azimuth_fname,
                          out_fnames=None, header_slope_fname=None):
    """
    Computes the self shadow mask, slope, aspect, incident, exiting,
    azimuth incident, azimuth exiting and relative slope angles.

    :param acquisition:
        An instance of an acquisition object.

    :param dsm_fname:
        A string containing the full file path name to the Digital
        Surface Model to be used in deriving the surface angles.

    :param margins:
        An object with members top, bottom, left and right giving the
        size of the margins (in pixels) which have been added to the
        corresponding sides of dsm.

    :param solar_zenith_fname:
        A string containing the full file path name to the solar
        zenith angle image.

    :param solar_azimuth_fname:
        A string containing the full file path name to the solar
        azimuth angle image.

    :param satellite_view_fname:
        A string containing the full file path name to the satellite
        view angle image.

    :param satellite_azimuth_fname:
        A string containing the full file path name to the satellite
        azimuth angle image.

    :param out_fnames:
        A list of length 8 containing full file path names for each
        of the outputs in the following order:

            * Self shadow mask
            * Slope
            * Aspect
            * Incident angle
            * Exiting angle
            * Azimuth incident angle
            * Azimuth exiting angle
            * Relative slope

        Default is None, in which case default filenames will be
        used and the results written to the current working
        directory.

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

    # Read the DSM and angle arrays into memory
    dsm = read_img(dsm_fname)
    solar_zenith = read_img(solar_zenith_fname)
    solar_azimuth = read_img(solar_azimuth_fname)
    satellite_view = read_img(satellite_view_fname)
    satellite_azimuth = read_img(satellite_azimuth_fname)

    # Define Top, Bottom, Left, Right pixel margins
    pixel_buf = ImageMargins(margins)

    # Compute self shadow, slope and various other angles
    slope_results = run_slope(acquisition, dsm, solar_zenith, satellite_view,
                              solar_azimuth, satellite_azimuth, pixel_buf,
                              is_utm, spheroid)

    # Output the results
    slope_results.write_arrays(out_fnames=out_fnames, geobox=geobox)

    if header_slope_fname:
        write_header_slope_file(header_slope_fname, pixel_buf, geobox)


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


def self_shadow(incident_fname, exiting_fname, self_shadow_out_fname):
    """
    Computes the self shadow mask.

    :param incident_fname:
        A string containing the full file path name to the incident
        angle image.

    :param exiting_fname:
        A string containing the full file path name to the exiting
        angle image.

    :param self_shadow_out_fname:
        A string containing the full file path name to be used for
        writing the self shadow mask on disk.

    :return:
        None. Output is written to disk.
    """

    with rasterio.open(incident_fname) as inc_ds,\
        rasterio.open(exiting_fname) as exi_ds:

        # Retrieve a geobox and image info
        geobox = GriddedGeoBox.from_rio_dataset(inc_ds)
        cols, rows = geobox.get_shape_xy()
        prj = geobox.crs.ExportToWkt()
        geoT = geobox.affine.to_gdal()

        # Initialise the output file
        drv = gdal.GetDriverByName("ENVI")
        out_dtype = gdal.GDT_Byte
        nbands = 1
        outds = drv.Create(self_shadow_out_fname, cols, rows, nbands,
                           out_dtype)
        outds.SetProjection(prj)
        outds.SetGeoTransform(geoT)
        outband = outds.GetRasterBand(1)

        # Initialise the tiling scheme for processing
        if X_TILE is None:
            X_TILE = cols
        if Y_TILE is None:
            Y_TILE = 1
        tiles = tiling.generate_tiles(cols, rows, X_TILE, Y_TILE,
                                      Generator=False)

        # Loop over each tile
        for tile in tiles:
            # Row and column start and end locations
            ystart = tile[0][0]
            xstart = tile[1][0]
            yend = tile[0][1]
            xend = tile[1[0]

            # Tile size
            ysize = yend - ystart
            xsize = xend - xstart

            # Read the data for the current tile
            inc = numpy.radians(inc_ds.read_band(1, window=tile, masked=False))
            exi = numpy.radians(exi_ds.read_band(1, window=tile, masked=False))

            # Process the tile
            mask = numpy.ones((rows, cols), dtype='uint8')
            mask[inc <= 0.0] = 0
            mask[exi <= 0.0] = 0

            # Write the current tile to disk
            outband.WriteArray(mask, xstart, ystart)
            outband.FlushCache()

        # Close the file
        outband = None
        outds = None
