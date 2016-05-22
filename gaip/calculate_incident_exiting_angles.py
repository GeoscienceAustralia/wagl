#!/usr/bin/env python

import numpy
import rasterio

from eotools import tiling
from gaip import as_array
from gaip import GriddedGeoBox
from gaip import exiting_angle
from gaip import incident_angle


def incident_angles(solar_zenith_fname, solar_azimuth_fname, slope_fname,
                   aspect_fname, incident_out_fname,
                   azimuth_incident_out_fname, x_tile=None, y_tile=None):
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

    :param x_tile:
        Defines the tile size along the x-axis. Default is None which
        equates to all elements along the x-axis.

    :param y_tile:
        Defines the tile size along the y-axis. Default is None which
        equates to all elements along the y-axis.

    :return:
        None. Outputs are written to disk.
    """

    # Retrieve a geobox
    with rasterio.open(solar_zenith_fname) as ds:
        geobox = GriddedGeoBox.from_rio_dataset(ds)
        cols, rows = geobox.get_shape_xy()
        crs = geobox.crs.ExportToWkt()

    # Initialise the output files
    out_dtype = 'float32'
    no_data = -999
    kwargs = {'driver': 'GTiff',
              'width': cols,
              'height': rows,
              'count': 1,
              'crs': crs,
              'transform': geobox.affine,
              'dtype': out_dtype,
              'nodata': no_data,
              'compress': 'deflate',
              'zlevel': 1,
              'predictor': 3}
    outds_incident = rasterio.open(incident_out_fname, 'w', **kwargs)
    outds_azi_incident = rasterio.open(azimuth_incident_out_fname, 'w',
                                       **kwargs)

    # Open the files for reading
    with rasterio.open(solar_zenith_fname) as sol_z_ds, \
        rasterio.open(solar_azimuth_fname) as sol_azi_ds, \
        rasterio.open(slope_fname) as slope_ds, \
        rasterio.open(aspect_fname) as aspect_ds:

        # Initialise the tiling scheme for processing
        if x_tile is None:
            x_tile = cols
        if y_tile is None:
            y_tile = rows
        tiles = tiling.generate_tiles(cols, rows, x_tile, y_tile,
                                      generator=False)

        # Loop over each tile
        for tile in tiles:
            # Row and column start and end locations
            ystart = tile[0][0]
            xstart = tile[1][0]
            yend = tile[0][1]
            xend = tile[1][1]

            # Tile size
            ysize = yend - ystart
            xsize = xend - xstart

            # Read the data for the current tile
            # Convert to required datatype and transpose
            sol_zen = as_array(sol_z_ds.read(1, window=tile,
                                                  masked=False),
                               dtype=numpy.float32, transpose=True)
            sol_azi = as_array(sol_azi_ds.read(1, window=tile,
                                                    masked=False),
                               dtype=numpy.float32, transpose=True)
            slope = as_array(slope_ds.read(1, window=tile, masked=False),
                             dtype=numpy.float32, transpose=True)
            aspect = as_array(aspect_ds.read(1, window=tile,
                                                  masked=False),
                              dtype=numpy.float32, transpose=True)

            # Initialise the work arrays
            incident = numpy.zeros((ysize, xsize), dtype='float32')
            azi_incident = numpy.zeros((ysize, xsize), dtype='float32')

            # Process the current tile
            incident_angle(xsize, ysize, sol_zen, sol_azi, slope, aspect,
                           incident.transpose(), azi_incident.transpose())

            # Write the current tile to disk
            outds_incident.write(incident, 1, window=tile)
            outds_azi_incident.write(azi_incident, 1, window=tile)

    # Close the files to complete the writing
    outds_incident.close()
    outds_azi_incident.close()


def exiting_angles(satellite_view_fname, satellite_azimuth_fname, slope_fname,
                  aspect_fname, exiting_out_fname, azimuth_exiting_out_fname,
                  x_tile=None, y_tile=None):
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

    :param x_tile:
        Defines the tile size along the x-axis. Default is None which
        equates to all elements along the x-axis.

    :param y_tile:
        Defines the tile size along the y-axis. Default is None which
        equates to all elements along the y-axis.

    :return:
        None. Outputs are written to disk.
    """

    # Retrieve a geobox
    with rasterio.open(satellite_view_fname) as ds:
        geobox = GriddedGeoBox.from_rio_dataset(ds)
        cols, rows = geobox.get_shape_xy()
        crs = geobox.crs.ExportToWkt()

    # Initialise the output files
    out_dtype = 'float32'
    no_data = -999
    kwargs = {'driver': 'GTiff',
              'width': cols,
              'height': rows,
              'count': 1,
              'crs': crs,
              'transform': geobox.affine,
              'dtype': out_dtype,
              'nodata': no_data,
              'compress': 'deflate',
              'zlevel': 1,
              'predictor': 3}
    outds_exiting = rasterio.open(exiting_out_fname, 'w', **kwargs)
    outds_azi_exiting = rasterio.open(azimuth_exiting_out_fname, 'w', **kwargs)

    # Open the files for reading
    with rasterio.open(satellite_view_fname) as sat_view_ds, \
        rasterio.open(satellite_azimuth_fname) as sat_azi_ds, \
        rasterio.open(slope_fname) as slope_ds, \
        rasterio.open(aspect_fname) as aspect_ds:

        # Initialise the tiling scheme for processing
        if x_tile is None:
            x_tile = cols
        if y_tile is None:
            y_tile = rows
        tiles = tiling.generate_tiles(cols, rows, x_tile, y_tile,
                                      generator=False)

        # Loop over each tile
        for tile in tiles:
            # Row and column start and end locations
            ystart = tile[0][0]
            xstart = tile[1][0]
            yend = tile[0][1]
            xend = tile[1][1]

            # Tile size
            ysize = yend - ystart
            xsize = xend - xstart

            # Read the data for the current tile
            # Convert to required datatype and transpose
            sat_view = as_array(sat_view_ds.read(1, window=tile,
                                                      masked=False),
                                dtype=numpy.float32, transpose=True)
            sat_azi = as_array(sat_azi_ds.read(1, window=tile,
                                                    masked=False),
                               dtype=numpy.float32, transpose=True)
            slope = as_array(slope_ds.read(1, window=tile, masked=False),
                             dtype=numpy.float32, transpose=True)
            aspect = as_array(aspect_ds.read(1, window=tile,
                                                  masked=False),
                              dtype=numpy.float32, transpose=True)

            # Initialise the work arrays
            exiting = numpy.zeros((ysize, xsize), dtype='float32')
            azi_exiting = numpy.zeros((ysize, xsize), dtype='float32')

            # Process the current tile
            exiting_angle(xsize, ysize, sat_view, sat_azi, slope, aspect,
                          exiting.transpose(), azi_exiting.transpose())

            # Write the current to disk
            outds_exiting.write(exiting, 1, window=tile)
            outds_azi_exiting.write(azi_exiting, 1, window=tile)

    # Close the files to complete the writing
    outds_exiting.close()
    outds_azi_exiting.close()


def relative_azimuth_slope(azimuth_incident_fname, azimuth_exiting_fname,
                           relative_azimuth_slope_out_fname,
                           x_tile=None, y_tile=None):
    """
    Calculates the relative azimuth angle on the slope surface.

    :param azimuth_incident_fname:
        A string containing the full file path name to the azimuth
        incident angle image.

    :param azimuth_exiting_fname:
        A string containing the full file path name to the azimuth
        exiting angle image.

    :param relative_azimuth_slope_out_fname:
        A string containing the full file path name to be used for
        writing the relative azimuth angle image on disk.

    :param x_tile:
        Defines the tile size along the x-axis. Default is None which
        equates to all elements along the x-axis.

    :param y_tile:
        Defines the tile size along the y-axis. Default is None which
        equates to all elements along the y-axis.

    :return:
        None. Output is written to disk.
    """

    with rasterio.open(azimuth_incident_fname) as azi_inc_ds,\
        rasterio.open(azimuth_exiting_fname) as azi_exi_ds:

        # Retrieve a geobox and image info
        geobox = GriddedGeoBox.from_rio_dataset(azi_inc_ds)
        cols, rows = geobox.get_shape_xy()
        crs = geobox.crs.ExportToWkt()

        # Initialise the output file
        out_dtype = 'float32'
        no_data = -999
        kwargs = {'driver': 'GTiff',
                  'width': cols,
                  'height': rows,
                  'count': 1,
                  'crs': crs,
                  'transform': geobox.affine,
                  'dtype': out_dtype,
                  'nodata': no_data,
                  'compress': 'deflate',
                  'zlevel': 1,
                  'predictor': 3}
        outds = rasterio.open(relative_azimuth_slope_out_fname, 'w', **kwargs)

        # Initialise the tiling scheme for processing
        if x_tile is None:
            x_tile = cols
        if y_tile is None:
            y_tile = rows
        tiles = tiling.generate_tiles(cols, rows, x_tile, y_tile,
                                      generator=False)

        # Loop over each tile
        for tile in tiles:
            # Read the data for the current tile
            azi_inc = azi_inc_ds.read(1, window=tile, masked=False)
            azi_exi = azi_exi_ds.read(1, window=tile, masked=False)

            # Process the tile
            rel_azi = azi_inc - azi_exi
            rel_azi[rel_azi <= -180.0] += 360.0
            rel_azi[rel_azi > 180.0] -= 360.0

            # Write the current tile to disk
            outds.write(rel_azi, 1, window=tile)

        # Close the file
        outds.close()
