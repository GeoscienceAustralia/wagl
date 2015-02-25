"""
Shadow Calculations
-------------------
"""

import gdal
import numpy
import rasterio

from EOtools import tiling
from gaip import GriddedGeoBox

X_TILE = None
Y_TILE = 100


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
            x_tile = cols
        else:
            x_tile = X_TILE
        if Y_TILE is None:
            y_tile = 1
        else:
            y_tile = Y_TILE
        tiles = tiling.generate_tiles(cols, rows, x_tile, y_tile,
                                      Generator=False)

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
            mask[inc <= 0.0] = 0
            mask[exi <= 0.0] = 0

            # Write the current tile to disk
            outband.WriteArray(mask, xstart, ystart)
            outband.FlushCache()

        # Close the file
        outband = None
        outds = None
