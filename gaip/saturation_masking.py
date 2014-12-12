import os
import logging
import numpy
import numexpr
from pqa_result import PQAResult
import constants


def SaturationMask(band_array, under_sat=1, over_sat=255, use_numexpr=True):
    """
    Creates a saturation mask for a single band.

    Under and over saturated pixels are masked.

    :param band_array:
        A 2D numpy array containing a single Landsat band.

    :param under_sat:
        Value used to detect under saturation. Default is 1.

    :param over_sat:
        Value used to detect over saturation. Default is 255.

    :param use_numexpr:
        If True (default) numexpression will be used to evalute the mask.
        Otherwise standard NumPy will be used to evaluate the mask.

    :return:
        A 2D Numpy Boolean array with True == Unsaturated.
    """

    if len(band_array) == 0: return None
    assert type(band_array) == numpy.ndarray, 'Input is not valid'

    if use_numexpr:
        logging.debug('numexpr used: numexpr.evaluate("(band_array != %i) & (band_array != %i)")'%(under_sat,over_sat))
        mask = numexpr.evaluate("(band_array != under_sat) & (band_array != over_sat)")
    else:
        logging.debug('numpy used: (band_array != %i) & (band_array != %i)'%(under_sat,over_sat))
        mask = (band_array != under_sat) & (band_array != over_sat)
    
    logging.debug('saturation mask computed')
    return mask


def setSaturationBits(l1t_stack, pq_const, result):
    logging.debug('setSaturationBits() called')
    band_list = pq_const.saturation_bands
    full_band_list = pq_const.available_bands
    band_index_list = pq_const.getArrayBandLookup(band_list)

    bit_index_list = pq_const.saturation_bits
    logging.debug('bit_index_list = %s', bit_index_list)

    for band in band_list:
        if band not in full_band_list:
            logging.warning('Ignoring invalid band number: %d', band)
            continue

        band_index = full_band_list.index(band)
        bit_index = bit_index_list[band_list.index(band)]
        logging.debug('Processing saturation for band = %d, band_index = %d, bit_index = %d', band, band_index, bit_index)

        band_array = l1t_stack[band_index]

        # Use numpy for single bands run in parallel - numexpr is not thread safe
        mask = SaturationMask(band_array, (len(band_list) > 1 or CONFIG.sequential))

        result.set_mask(mask, bit_index)

        # *** This will need to change. Tests not run will not be set. ***
        # Copy results for first thermal band to second one if there is only one available
        if bit_index == 5 and 6 not in bit_index_list: # If this is thermal band 1 and there is no thermal band 2 in this dataset
            bit_index = 6
            logging.debug('Copying thermal band mask to bit %d', bit_index)
            result.set_mask(mask, bit_index)
