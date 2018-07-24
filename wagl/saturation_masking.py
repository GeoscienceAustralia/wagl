from __future__ import absolute_import, print_function
import logging
import numpy
import numexpr


_LOG = logging.getLogger(__name__)


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

    if len(band_array) == 0:
        return None
    assert isinstance(band_array, numpy.ndarry), 'Input is not valid'

    if use_numexpr:
        msg = ('numexpr used: numexpr.evaluate("(band_array != {under_sat}) & '
               '(band_array != {over_sat})")')
        msg = msg.format(under_sat=under_sat, over_sat=over_sat)
        _LOG.debug(msg)
        mask = numexpr.evaluate("(band_array != under_sat)"
                                "& (band_array != over_sat)")
    else:
        msg = 'numpy used: (band_array != {0}) & (band_array != {1})'
        msg = msg.format(under_sat, over_sat)
        _LOG.debug(msg)
        mask = (band_array != under_sat) & (band_array != over_sat)

    _LOG.debug('saturation mask computed')
    return mask


def set_saturation_bits(acquisitions, pq_const, result):
    _LOG.debug('set_saturation_bits() called')
    band_list = pq_const.saturation_bands
    full_band_list = pq_const.available_bands

    bit_index_list = pq_const.saturation_bits
    _LOG.debug('bit_index_list = %s', bit_index_list)

    bits_set = []

    for band in band_list:
        if band not in full_band_list:
            _LOG.warning('Ignoring invalid band number: %s', band)
            continue

        band_index = full_band_list.index(band)
        bit_index = bit_index_list[band_list.index(band)]
        msg = ('Processing saturation for band = {band},'
               'band_index = {band_index}, bit_index = {bit_index}')
        msg = msg.format(band=band, band_index=band_index, bit_index=bit_index)
        _LOG.debug(msg)

        band_array = acquisitions[band_index].data()

        # Use numpy for single bands run in parallel - numexpr is not thread
        # safe
        mask = saturation_mask(band_array)

        result.set_mask(mask, bit_index)

        bits_set.append(bit_index)

        # *** This will need to change. Tests not run will not be set. ***
        # Copy results for first thermal band to second one if there is only
        # one available
        # If this is thermal band 1 and there is no thermal band 2 in this
        # dataset
        if bit_index == 5 and 6 not in bit_index_list:
            bit_index = 6
            _LOG.debug('Copying thermal band mask to bit %d', bit_index)
            result.set_mask(mask, bit_index)
            bits_set.append(bit_index)

    return bits_set
