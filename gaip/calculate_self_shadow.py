"""
Shadow Calculations
-------------------
"""

import numpy
import rasterio

from eotools import tiling
from gaip import GriddedGeoBox


def self_shadow(incident_fname, exiting_fname, self_shadow_out_fname,
                x_tile=None, y_tile=None):
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

    :param x_tile:
        Defines the tile size along the x-axis. Default is None which
        equates to all elements along the x-axis.

    :param y_tile:
        Defines the tile size along the y-axis. Default is None which
        equates to all elements along the y-axis.

    :return:
        None. Output is written to disk.
    """

    with rasterio.open(incident_fname) as inc_ds,\
        rasterio.open(exiting_fname) as exi_ds:

        # Retrieve a geobox and image info
        geobox = GriddedGeoBox.from_rio_dataset(inc_ds)
        cols, rows = geobox.get_shape_xy()

        # Initialise the output file
        outds = tiling.TiledOutput(self_shadow_out_fname, cols, rows,
                                   geobox=geobox)

        # Initialise the tiling scheme for processing
        if x_tile is None:
            x_tile = cols
        if y_tile is None:
            y_tile = rows
        tiles = tiling.generate_tiles(cols, rows, x_tile, y_tile,
                                      generator=False)

        # Loop over each tile
        for tile in tiles:
            # Row and column start locations
            ystart = tile[0][0]
            xstart = tile[1][0]
            yend = tile[0][1]
            xend = tile[1][1]

            # Tile size
            ysize = yend - ystart
            xsize = xend - xstart

            # Read the data for the current tile
            inc = numpy.radians(inc_ds.read_band(1, window=tile, masked=False))
            exi = numpy.radians(exi_ds.read_band(1, window=tile, masked=False))

            # Process the tile
            mask = numpy.ones((ysize, xsize), dtype='uint8')
            mask[numpy.cos(inc) <= 0.0] = 0
            mask[numpy.cos(exi) <= 0.0] = 0

            # Write the current tile to disk
            outds.write_tile(mask, tile)

        # Close the file
        outds.close()
        #outband = None
        outds = None
