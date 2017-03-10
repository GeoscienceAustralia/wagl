from __future__ import absolute_import
import errno
import os
import numpy
import numexpr
import gdal

from gaip import tiling
from gaip import endmembers
from gaip import GriddedGeoBox
from gaip import LandsatAcquisition
from gaip import stack_data

#from datacube.api.model import BANDS
#from datacube.api.model import DatasetType, Satellite

"""
Utility functions used in fractional cover. These should not be used outside of this package as they
may be moved, renamed or deleted in the future. If you need to add more functions that only for internal
use, this is the place to put them.
"""

def fractional_cover(acquisitions, x_tile, y_tile, out_fnames):
    """
    Given a list of spectral acquisitions compute the fractional
    components.

    :param acquisitions:
        A list of acquisition class objects that will be run through
        the terrain correction workflow.

    :param x_tile:
        Defines the tile size along the x-axis. Default is None which
        equates to all elements along the x-axis.

    :param y_tile:
        Defines the tile size along the y-axis. Default is None which
        equates to all elements along the y-axis.

    :param out_fnames:
        A list of strings with each containing the full file path name
        to a location on disk.
        The order should be:

            * Photosynthetic vegetation
            * Non-photosynthetic vegetation
            * Bare soil
            * Unmixing error

    :return:
        None. Outputs are written to disk.

    :Notes:
        The number of outputs, the order and the type of fractional
        component can change with the algorithm used.
    """
    # Define the output datatype and format
    out_dtype = gdal.GDT_Byte
    fmt = "GTiff"
    out_no_data = 0

    # Compute the geobox and get the array dimensions from the 1st
    # acquisition
    if isinstance(acquisitions, LandsatAcquisition):
        geobox = acquisitions[0].gridded_geo_box()
        cols, rows = geobox.get_shape_xy()

        # Define int16 zero and no data value
        zero = numpy.int16(0)
        no_data = acquisitions[0].no_data 
        if no_data is None:
            no_data = -999
        
        # Initialise the output files
        outds_pv = tiling.TiledOutput(out_fnames[0], cols, rows, geobox=geobox,
                                      dtype=out_dtype, nodata=out_no_data,
                                      fmt=fmt)
        outds_npv = tiling.TiledOutput(out_fnames[1], cols, rows,
                                       geobox=geobox, dtype=out_dtype,
                                       nodata=out_no_data, fmt=fmt)
        outds_bs = tiling.TiledOutput(out_fnames[2], cols, rows, geobox=geobox,
                                      dtype=out_dtype, nodata=out_no_data,
                                      fmt=fmt)
        outds_ue = tiling.TiledOutput(out_fnames[3], cols, rows, geobox=geobox,
                                      dtype=out_dtype, nodata=out_no_data,
                                      fmt=fmt)
    else:
        # To avoid confusion between the Acquisition and StackedDataset
        # objects we'll alias acquisitions to sd (StackedDataset)
        sd = acquisitions

        geobox = GriddedGeoBox.from_gdal_dataset(gdal.Open(sd.fname))
        cols = sd.samples
        rows = sd.lines
        zero = numpy.int16(0)
        no_data = sd.no_data
        if no_data is None:
            no_data = -999

        # Define the output file
        outds = tiling.TiledOutput(out_fnames, cols, rows, geobox=geobox,
                                   dtype=out_dtype, nodata=out_no_data,
                                   fmt=fmt, bands=4)

    # Initialise the tiling scheme for processing
    if x_tile is None:
        x_tile = cols
    if y_tile is None:
        y_tile = rows
    tiles = generate_tiles(cols, rows, x_tile, y_tile)


    # Scarth 20090810 14:06:35 CEST
    # Define the weight of the sum to one constraint
    # This value determined how well the resulting fractions will sum to 100%
    # I typically determine this by running the unmixing against field data for
    # a number of values, picking the best one
    sum_to_one_weight = endmembers.sum_weight('2014_07_23')


    # Scarth 20090810 14:06:35 CEST
    # 2009v gives green, dead, bare1 and bare2
    # 2012v gives green, dead, bare1 and bare2
    # 2013v gives green, dead1, dead2 and bare fractions
    # Note the last row is the sum to one constraint value
    endmembers_array = endmembers.endmember_version('2014_07_23')

    # Define separate loops for the different data sources
    if isinstance(acquisitions, LandsatAcquisition):
        # Loop over each tile
        for tile in tiles:
            # Read the data for the current tile from acquisitions
            stack, _ = stack_data(acquisitions, window=tile)

            # set no data values to zero
            #numpy.maximum(stack, zero, out=stack) # Only works if negative
            wh_any = numpy.any(numexpr.evaluate("stack == no_data"), axis=0)
            stack[:, wh_any] = zero

            # Compute the fractions
            (green, dead1, dead2, bare,
             err) = unmix(stack[0], stack[1], stack[2], stack[3], stack[4],
                          sum_to_one_weight, endmembers_array)

            # Find any occurences of an unmixing error
            # If an pixel is in error then all pixels for that location are
            # counted as an error
            wh_unmix_err = numexpr.evaluate("(green == -10) |"
                                            "(dead1 == -10) |"
                                            "(dead2 == -10) |"
                                            "(bare == -10)")

            # scale the results and clip the range to (0, 255)
            green = numexpr.evaluate("green / 0.01 + 100")
            green[wh_any] = out_no_data
            green[wh_unmix_err] = out_no_data
            numpy.clip(green, a_min=0, a_max=255, out=green)

            dead = numexpr.evaluate("(dead1 + dead2) / 0.01 + 100")
            dead[wh_any] = out_no_data
            dead[wh_unmix_err] = out_no_data
            numpy.clip(dead, a_min=0, a_max=255, out=dead)

            bare = numexpr.evaluate("bare / 0.01 + 100")
            bare[wh_any] = out_no_data
            bare[wh_unmix_err] = out_no_data
            numpy.clip(bare, a_min=0, a_max=255, out=bare)

            err = numexpr.evaluate("err + 100")
            err[wh_any] = out_no_data
            err[wh_unmix_err] = out_no_data
            numpy.clip(err, a_min=0, a_max=255, out=err)

            # Output to disk
            outds_pv.write_tile(green.astype('uint8'), tile)
            outds_npv.write_tile(dead.astype('uint8'), tile)
            outds_bs.write_tile(bare.astype('uint8'), tile)
            outds_ue.write_tile(err.astype('uint8'), tile)
    else:
        # Get the sensor and construct a Satellite object
        satellite = Satellite[os.path.basename(sd.fname).split('_')[0]]

        # Get the required bands for the unmixing algorithm
        bands = BANDS[DatasetType.ARG25, satellite]
        bands = [bands.GREEN.value, bands.RED.value,
                 bands.NEAR_INFRARED.value,
                 bands.SHORT_WAVE_INFRARED_1.value,
                 bands.SHORT_WAVE_INFRARED_2.value]

        # ******* QDERM TEST **********
        #bands = [2,3,4,5,6]

        # Loop over each tile
        for tile in tiles:
            # Read the data for the current tile from acquisitions
            stack = sd.read_tile(tile, raster_bands=bands)

            # For some reason we might get an un-writeable array
            stack.flags['WRITEABLE'] = True

            # set no data values to zero
            #numpy.maximum(stack, zero, out=stack)
            wh_any = numpy.any(numexpr.evaluate("stack == no_data"), axis=0)
            stack[:, wh_any] = zero

            # Compute the fractions
            (green, dead1, dead2, bare,
             err) = unmix(stack[0], stack[1], stack[2], stack[3], stack[4],
                          sum_to_one_weight, endmembers_array)

            # Find any occurences of an unmixing error
            # If an pixel is in error then all pixels for that location are
            # counted as an error
            wh_unmix_err = numexpr.evaluate("(green == -10) |"
                                            "(dead1 == -10) |"
                                            "(dead2 == -10) |"
                                            "(bare == -10)")

            # scale the results and clip the range to (0, 255)
            green = numexpr.evaluate("green / 0.01 + 100")
            green[wh_any] = out_no_data
            green[wh_unmix_err] = out_no_data
            numpy.clip(green, a_min=0, a_max=255, out=green)

            dead = numexpr.evaluate("(dead1 + dead2) / 0.01 + 100")
            dead[wh_any] = out_no_data
            dead[wh_unmix_err] = out_no_data
            numpy.clip(dead, a_min=0, a_max=255, out=dead)

            bare = numexpr.evaluate("bare / 0.01 + 100")
            bare[wh_any] = out_no_data
            bare[wh_unmix_err] = out_no_data
            numpy.clip(bare, a_min=0, a_max=255, out=bare)

            err = numexpr.evaluate("err + 100")
            err[wh_any] = out_no_data
            err[wh_unmix_err] = out_no_data
            numpy.clip(err, a_min=0, a_max=255, out=err)

            # Output to disk
            outds.write_tile(bare.astype('uint8'), tile, raster_band=1)
            outds.write_tile(green.astype('uint8'), tile, raster_band=2)
            outds.write_tile(dead.astype('uint8'), tile, raster_band=3)
            outds.write_tile(err.astype('uint8'), tile, raster_band=4)

    # Close the files to complete the writing
    if isinstance(acquisitions, LandsatAcquisition):
        outds_pv.close()
        outds_npv.close()
        outds_bs.close()
        outds_ue.close()
    else:
        outds.close()


def unmix(green, red, nir, swir1, swir2, sum_to_one_weight, endmembers_array):
    # NNLS Unmixing v1.0
    # Scarth 20090810 14:06:35 CEST
    # This implements a constrained unmixing process to recover the fraction images from
    # a synthetic reflectance generated from a large number of interactive
    # terms produced from the original and log-transformed landsat bands

    # GA wrapped and modified version of Scarth 20090810 14:06:35 CEST

    band2 = numexpr.evaluate("(1.0 + green) * 0.0001")
    band3 = numexpr.evaluate("(1.0 + red) * 0.0001")
    band4 = numexpr.evaluate("(1.0 + nir) * 0.0001")
    band5 = numexpr.evaluate("(1.0 + swir1) * 0.0001")
    band7 = numexpr.evaluate("(1.0 + swir2) * 0.0001")

    #b_logs = numexpr.evaluate("log(subset)")
    logb2 = numexpr.evaluate("log(band2)")
    logb3 = numexpr.evaluate("log(band3)")
    logb4 = numexpr.evaluate("log(band4)")
    logb5 = numexpr.evaluate("log(band5)")
    logb7 = numexpr.evaluate("log(band7)")

    b2b3  = numexpr.evaluate("band2 * band3")
    b2b4  = numexpr.evaluate("band2 * band4")
    b2b5  = numexpr.evaluate("band2 * band5")
    b2b7  = numexpr.evaluate("band2 * band7")
    b2lb2 = numexpr.evaluate("band2 * logb2")
    b2lb3 = numexpr.evaluate("band2 * logb3")
    b2lb4 = numexpr.evaluate("band2 * logb4")
    b2lb5 = numexpr.evaluate("band2 * logb5")
    b2lb7 = numexpr.evaluate("band2 * logb7")

    b3b4  = numexpr.evaluate("band3 * band4")
    b3b5  = numexpr.evaluate("band3 * band5")
    b3b7  = numexpr.evaluate("band3 * band7")
    b3lb2 = numexpr.evaluate("band3 * logb2")
    b3lb3 = numexpr.evaluate("band3 * logb3")
    b3lb4 = numexpr.evaluate("band3 * logb4")
    b3lb5 = numexpr.evaluate("band3 * logb5")
    b3lb7 = numexpr.evaluate("band3 * logb7")

    b4b5  = numexpr.evaluate("band4 * band5")
    b4b7  = numexpr.evaluate("band4 * band7")
    b4lb2 = numexpr.evaluate("band4 * logb2")
    b4lb3 = numexpr.evaluate("band4 * logb3")
    b4lb4 = numexpr.evaluate("band4 * logb4")
    b4lb5 = numexpr.evaluate("band4 * logb5")
    b4lb7 = numexpr.evaluate("band4 * logb7")

    b5b7  = numexpr.evaluate("band5 * band7")
    b5lb2 = numexpr.evaluate("band5 * logb2")
    b5lb3 = numexpr.evaluate("band5 * logb3")
    b5lb4 = numexpr.evaluate("band5 * logb4")
    b5lb5 = numexpr.evaluate("band5 * logb5")
    b5lb7 = numexpr.evaluate("band5 * logb7")

    b7lb2 = numexpr.evaluate("band7 * logb2")
    b7lb3 = numexpr.evaluate("band7 * logb3")
    b7lb4 = numexpr.evaluate("band7 * logb4")
    b7lb5 = numexpr.evaluate("band7 * logb5")
    b7lb7 = numexpr.evaluate("band7 * logb7")

    lb2lb3 = numexpr.evaluate("logb2 * logb3")
    lb2lb4 = numexpr.evaluate("logb2 * logb4")
    lb2lb5 = numexpr.evaluate("logb2 * logb5")
    lb2lb7 = numexpr.evaluate("logb2 * logb7")

    lb3lb4 = numexpr.evaluate("logb3 * logb4")
    lb3lb5 = numexpr.evaluate("logb3 * logb5")
    lb3lb7 = numexpr.evaluate("logb3 * logb7")

    lb4lb5 = numexpr.evaluate("logb4 * logb5")
    lb4lb7 = numexpr.evaluate("logb4 * logb7")

    lb5lb7 = numexpr.evaluate("logb5 * logb7")

    band_ratio1 = numexpr.evaluate("(band4 - band3) / (band4 + band3)")
    band_ratio2 = numexpr.evaluate("(band4 - band5) / (band4 + band5)")
    band_ratio3 = numexpr.evaluate("(band5 - band3) / (band5 + band3)")
    band_ratio4 = numexpr.evaluate("(band3 - band2) / (band3 + band2)")

    # The 2009_08_10 and 2012_12_07 versions use a different interactive
    # terms array compared to the 2013_01_08 version
    # 2013_01_08 uses 59 endmebers
    # 2009_08_10 uses 56 endmebers
    # 2012_12_07 uses 56 endmebers
    # 2014_07_23 uses 60 endmembers
    # TODO write an interface that can retrieve the correct
    # interactiveTerms array according to the specified version.

    interactive_terms = numpy.array([b2b3, b2b4, b2b5, b2b7, b2lb2, b2lb3,
                                     b2lb4, b2lb5, b2lb7, b3b4, b3b5, b3b7,
                                     b3lb2, b3lb3, b3lb4, b3lb5, b3lb7, b4b5,
                                     b4b7, b4lb2, b4lb3, b4lb4, b4lb5, b4lb7,
                                     b5b7, b5lb2, b5lb3, b5lb4, b5lb5, b5lb7,
                                     b7lb2, b7lb3, b7lb4, b7lb5, b7lb7, lb2lb3,
                                     lb2lb4, lb2lb5, lb2lb7, lb3lb4, lb3lb5,
                                     lb3lb7, lb4lb5, lb4lb7, lb5lb7, band2,
                                     band3, band4, band5, band7, logb2, logb3,
                                     logb4, logb5, logb7, band_ratio1,
                                     band_ratio2, band_ratio3, band_ratio4])

    # Now add the sum to one constraint to the interactive terms
    # First make a zero array of the right shape
    weighted_spectra = numpy.zeros((interactive_terms.shape[0] + 1,) +
                                   interactive_terms.shape[1:])
    # Insert the interactive terms
    weighted_spectra[:-1, ...] = interactive_terms
    # Last element is special weighting
    weighted_spectra[-1] = sum_to_one_weight

    in_null = 0.0001
    out_unmix_null = -10.0

    from gaip import unmiximage
    fractions = unmiximage.unmiximage(weighted_spectra, endmembers_array,
                                      in_null, out_unmix_null)

    # 2013v gives green, dead1, dead2 and bare fractions
    # the last band should be the unmixing error
    return fractions
