"""
Contiguity Mask
---------------
"""
from __future__ import absolute_import, print_function
import logging
import numpy

from scipy import ndimage
from idl_functions import histogram
from gaip.data import stack_data
from gaip.tiling import generate_tiles


def calc_contiguity_mask(acquisitions, platform_id):
    """
    Determines locations of null values.

    Null values for every band are located in order to create band
    contiguity.

    :param acquisitions:
        A `list` of `acquisition` objects.

    :param platform_id:
        A `str` containing the platform id as given by
        `acquisition.platform_id`.

    :return:
        A single ndarray determining band/pixel contiguity. 1 for
        contiguous, 0 for non-contiguous.

    :notes:
        Attempts to flag thermal anomolies for Landsat 5TM as well.
    """
    cols = acquisitions[0].samples
    rows = acquisitions[0].lines
    tiles = list(generate_tiles(cols, rows, cols))

    logging.debug('Determining pixel contiguity')
    # Create mask array with True for all pixels which are non-zero in all
    # bands
    mask = numpy.zeros((rows, cols), dtype='bool')

    for tile in tiles:
        idx = (slice(tile[0][0], tile[0][1]), slice(tile[1][0], tile[1][1]))
        stack, _ = stack_data(acquisitions, window=tile)
        mask[idx] = stack.all(0)

    # The following is only valid for Landsat 5 images
    logging.debug('calc_contiguity_mask: platform_id=%s', platform_id)
    if platform_id == 'LANDSAT_5':
        logging.debug('Finding thermal edge anomalies')
        # Apply thermal edge anomalies
        struct = numpy.ones((7, 7), dtype='bool')
        erode = ndimage.binary_erosion(mask, structure=struct)

        dims = mask.shape
        th_anom = numpy.zeros(dims, dtype='bool').flatten()

        pix_3buff_mask = mask - erode
        pix_3buff_mask[pix_3buff_mask > 0] = 1
        edge = pix_3buff_mask == 1

        low_sat = acquisitions[5].data() == 1
        low_sat_buff = ndimage.binary_dilation(low_sat, structure=struct)

        s = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
        low_sat, _ = ndimage.label(low_sat_buff, structure=s)

        labels = low_sat[edge]
        ulabels = numpy.unique(labels[labels > 0])

        # Histogram method, a lot faster
        mx = numpy.max(ulabels)
        h = histogram(low_sat, minv=0, maxv=mx, reverse_indices='ri')
        hist = h['histogram']
        ri = h['ri']

        for i in numpy.arange(ulabels.shape[0]):
            if hist[ulabels[i]] == 0:
                continue
            th_anom[ri[ri[ulabels[i]]:ri[ulabels[i] + 1]]] = True

        th_anom = ~(th_anom.reshape(dims))
        mask &= th_anom

    return mask


def set_contiguity_bit(l1t_acqs, platform_id, pq_const, pqa_result):
    """Set the contiguity bit."""
    mask = calc_contiguity_mask(l1t_acqs, platform_id)
    bit_index = pq_const.contiguity
    pqa_result.set_mask(mask, bit_index)
