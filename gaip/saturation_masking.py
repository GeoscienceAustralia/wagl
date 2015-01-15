import os
import logging
import numpy
import numexpr
from pqa_result import PQAResult
import constants


def saturation_mask(band_array, under_sat=1, over_sat=255, use_numexpr=True):
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
        msg = ('numexpr used: numexpr.evaluate("(band_array != {under_sat}) & '
               '(band_array != {over_sat})")')
        msg = msg.format(under_sat=under_sat, over_sat=over_sat)
        logging.debug(msg)
        mask = numexpr.evaluate("(band_array != under_sat)"
                                "& (band_array != over_sat)")
    else:
        msg = 'numpy used: (band_array != {0}) & (band_array != {1})'
        msg = msg.format(under_sat, over_sat)
        logging.debug(msg)
        mask = (band_array != under_sat) & (band_array != over_sat)
    
    logging.debug('saturation mask computed')
    return mask


def set_saturation_bits(l1t_stack, pq_const, result):
    logging.debug('set_saturation_bits() called')
    band_list = pq_const.saturation_bands
    full_band_list = pq_const.available_bands
    band_index_list = pq_const.getArrayBandLookup(band_list)

    bit_index_list = pq_const.saturation_bits
    logging.debug('bit_index_list = %s', bit_index_list)

    for band in band_list:
        if band not in full_band_list:
            logging.warning('Ignoring invalid band number: {}'.format(band))
            continue

        band_index = full_band_list.index(band)
        bit_index = bit_index_list[band_list.index(band)]
        msg = ('Processing saturation for band = {band},'
               'band_index = {band_index}, bit_index = {bit_index}')
        msg = msg.format(band=band, band_index=band_index, bit_index=bit_index)
        logging.debug(msg)

        band_array = l1t_stack[band_index]

        # Use numpy for single bands run in parallel - numexpr is not thread
        # safe
        mask = saturation_mask(band_array, (len(band_list) > 1
                                            or CONFIG.sequential))

        result.set_mask(mask, bit_index)

        # *** This will need to change. Tests not run will not be set. ***
        # Copy results for first thermal band to second one if there is only
        # one available
        # If this is thermal band 1 and there is no thermal band 2 in this
        # dataset
        if bit_index == 5 and 6 not in bit_index_list:
            bit_index = 6
            logging.debug('Copying thermal band mask to bit %d', bit_index)
            result.set_mask(mask, bit_index)
