#!/usr/bin/env python

import gdal
import numpy
import rasterio

from EOtools import tiling
from gaip import GriddedGeoBox
from gaip import _exiting_angle
from gaip import _incident_angle


X_TILE = None
Y_TILE = 100

def incident_angle(solar_zenith_fname, solar_azimuth_fname, slope_fname,
                   aspect_fname, incident_out_fname,
                   azimuth_incident_out_fname):
    """
    Calculates the incident angle and the azimuthal incident angle.

    :param solar_zenith_fname:
        A string containing the full file path name to the solar
        zenith angle image.

    :param solar_azimuth_fname:
        A string containing the full file path name to the solar
        azimuth angle image.

    :param slope_fname:
        A string containing the full file path name to the slope image.

    :param aspect_fname:
        A string containing the full file path name to the aspect image.

    :param incident_out_fname:
        A string containing the full file path name to be used for
        writing the incident angle image on disk.

    :param azimuth_incident_out_fname:
        A string containing the full file path name to be used for
        writing the azimuth incident angle image on disk.

    :return:
        None. Outputs are written to disk.
    """

    # Retrieve a geobox
    with rasterio.open(solar_zenith_fname) as ds:
        geobox = GriddedGeoBox.from_rio_dataset(ds)

    # Initialise the output files
    output_fnames = {'incident': incident_out_fname,
                     'azimuth_incident': azimuth_incident_out_fname}
    output_files = {}
    out_bands = {}
    drv = gdal.GetDriverByName("ENVI")
    out_dtype = gdal.GDT_Float32
    nbands = 1
    prj = geobox.crs.ExportToWkt()
    geoT = geobox.affine.to_gdal()
    for key in output_fnames:
        fname = output_fnames[key]
        output_files[key] = drv.Create(fname, cols, rows, nbands, out_dtype)
        output_files[key].SetProjection(prj)
        output_files[key].SetGeoTransform(geoT)
        out_bands[key] = output_files[key].GetRasterBand(1)
        out_bands[key].SetNoDataValue(-999)

    # Open the files for reading
    with rasterio.open(solar_zenith_fname) as sol_z_ds, \
        rasterio.open(solar_azimuth_fname) as sol_azi_ds, \
        rasterio.open(slope_fname) as slope_ds, \
        rasterio.open(aspect_fname) as aspect_ds:

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
           # Convert to required datatype and transpose
           sol_zen = as_array(sol_z_ds.read_band(1, window=tile, masked=False),
                              dtype=numpy.float32, transpose=True)
           sol_azi = as_array(sol_azi_ds.read_band(1, window=tile,
                                                   masked=False),
                              dtype=numpy.float32, transpose=True)
           slope = as_array(slope_ds.read_band(1, window=tile, masked=False),
                            dtype=numpy.float32, transpose=True)
           aspect = as_array(aspect_ds.read_band(1, window=tile, masked=False),
                             dtype=numpy.float32, transpose=True)

           # Initialise the work arrays
           incident = numpy.zeros((ysize, xsize), dtype='float32').transpose()
           azi_incident = numpy.zeros((ysize, xsize),
                                      dtype='float32').transpose()

           # Process the current tile
           _incident_angle(xsize, ysize, sol_zen, sol_azi, slope, aspect,
                           incident, azi_incident)

           # Write the current tile to disk
           out_bands['incident'].WriteArray(incident, xstart, ystart)
           out_bands['incident'].FlushCache()
           out_bands['azimuth_incident'].WriteArray(azi_incident, xstart,
                                                    ystart)
           out_bands['azimuth_incident'].FlushCache()

    # Close the files to complete the writing
    for key in output_fnames:
        out_bands[key] = None
        output_files[key] = None

    out_bands = None
    output_files = None


def exiting_angle(satellite_view_fname, satellite_azimuth_fname, slope_fname,
                  aspect_fname, exiting_out_fname, azimuth_exiting_out_fname):
    """
    Calculates the exiting angle and the azimuthal exiting angle.

    :param solar_zenith_fname:
        A string containing the full file path name to the solar
        zenith angle image.

    :param solar_azimuth_fname:
        A string containing the full file path name to the solar
        azimuth angle image.

    :param slope_fname:
        A string containing the full file path name to the slope image.

    :param aspect_fname:
        A string containing the full file path name to the aspect image.

    :param exiting_out_fname:
        A string containing the full file path name to be used for
        writing the exiting angle image on disk.

    :param azimuth_exiting_out_fname:
        A string containing the full file path name to be used for
        writing the azimuth exiting angle image on disk.

    :return:
        None. Outputs are written to disk.
    """

    # Retrieve a geobox
    with rasterio.open(satellite_view_fname) as ds:
        geobox = GriddedGeoBox.from_rio_dataset(ds)

    # Initialise the output files
    output_fnames = {'exiting': exiting_out_fname,
                     'azimuth_exiting': azimuth_exiting_out_fname}
    output_files = {}
    out_bands = {}
    drv = gdal.GetDriverByName("ENVI")
    out_dtype = gdal.GDT_Float32
    nbands = 1
    prj = geobox.crs.ExportToWkt()
    geoT = geobox.affine.to_gdal()
    for key in output_fnames:
        fname = output_fnames[key]
        output_files[key] = drv.Create(fname, cols, rows, nbands, out_dtype)
        output_files[key].SetProjection(prj)
        output_files[key].SetGeoTransform(geoT)
        out_bands[key] = output_files[key].GetRasterBand(1)
        out_bands[key].SetNoDataValue(-999)

    # Open the files for reading
    with rasterio.open(satellite_view_fname) as sat_view_ds, \
        rasterio.open(satellite_azimuth_fname) as sat_azi_ds, \
        rasterio.open(slope_fname) as slope_ds, \
        rasterio.open(aspect_fname) as aspect_ds:

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
           # Convert to required datatype and transpose
           sat_view = as_array(sat_view_ds.read_band(1, window=tile,
                                                     masked=False),
                               dtype=numpy.float32, transpose=True)
           sat_azi = as_array(sat_azi_ds.read_band(1, window=tile,
                                                   masked=False),
                              dtype=numpy.float32, transpose=True)
           slope = as_array(slope_ds.read_band(1, window=tile, masked=False),
                            dtype=numpy.float32, transpose=True)
           aspect = as_array(aspect_ds.read_band(1, window=tile, masked=False),
                             dtype=numpy.float32, transpose=True)

           # Initialise the work arrays
           exiting = numpy.zeros((ysize, xsize), dtype='float32').transpose()
           azi_exiting = numpy.zeros((ysize, xsize),
                                     dtype='float32').transpose()

           # Process the current tile
           _exiting_angle(xsize, ysize, sol_zen, sol_azi, slope, aspect,
                           exiting, azi_exiting)

           # Write the current to disk
           out_bands['exiting'].WriteArray(exiting, xstart, ystart)
           out_bands['exiting'].FlushCache()
           out_bands['azimuth_exiting'].WriteArray(azi_exiting, xstart,
                                                    ystart)
           out_bands['azimuth_exiting'].FlushCache()

    # Close the files to complete the writing
    for key in output_fnames:
        out_bands[key] = None
        output_files[key] = None

    out_bands = None
    output_files = None


def relative_azimuth(azimuth_incident_fname, azimuth_exiting_fname,
                     relative_azimuth_out_fname):
    """
    Calculates the relative azimuth angle on the slope surface.

    :param azimuth_incident_fname:
        A string containing the full file path name to the azimuth
        incident angle image.

    :param azimuth_exiting_fname:
        A string containing the full file path name to the azimuth
        exiting angle image.

    :param relative_slope_out_fname:
        A string containing the full file path name to be used for
        writing the relative azimuth angle image on disk.
    """

    with rasterio.open(azimuth_incident_fname) as azi_inc_ds,\
        rasterio.open(azimuth_exiting_fname) as azi_exi_ds:

        # Retrieve a geobox and image info
        geobox = GriddedGeoBox.from_rio_dataset(azi_inc_ds)
        cols, rows = geobox.get_shape_xy()
        prj = geobox.crs.ExportToWkt()
        geoT = geobox.affine.to_gdal()

        # Initialise the output file
        drv = gdal.GetDriverByName("ENVI")
        out_dtype = gdal.GDT_Float32
        nbands = 1
        outds = drv.Create(relative_azimuth_fname, cols, rows, nbands,
                           out_dtype)
        outds.SetProjection(prj)
        outds.SetGeoTransform(geoT)
        outband = outds.GetRasterBand(1)
        outband.SetNoDataValue(-999)

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
            azi_inc = azi_inc_ds.read_band(1, window=tile, masked=False)
            azi_exi = azi_exi_ds.read_band(1, window=tile, masked=False)

            # Process the tile
            rel_azi = azi_inc - azi_exi
            rel_azi[rel_azi <= -180.0] += 360.0
            rel_azi[rel_azi > 180.0] -= 360.0

            # Write the current tile to disk
            outband.WriteArray(rel_azi, xstart, ystart)
            outband.FlushCache()

        # Close the file
        outband = None
        outds = None
