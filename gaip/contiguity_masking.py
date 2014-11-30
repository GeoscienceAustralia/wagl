import os
import numpy
import logging
import numexpr
from scipy import ndimage
from IDL_functions import histogram


def calc_contiguity_mask(image_stack, spacecraft_id):
    """
    Determines locations of null values.

    Null values for every band are located in order to create band
    contiguity.

    :param image:
        An nD Numpy array of all bands (ordered).

    :param mask:
        Output array.

    :param slc_off:
        Whether to perform image scaling to flag potential
        saturated/non-contiguous pixels. Should only be applied to Non
        landsat 7 products. (L7 products prior to slc-off could be run).
        Default is False.

    :return:
        A single ndarray determining band/pixel contiguity. 1 for
        contiguous, 0 for non-contiguous.

    :notes:
        Attempts to flag thermal anomolies for Landsat 5TM as well.
    """

    if len(image_stack) == 0:
        return None

    assert type(image_stack[0]) == numpy.ndarray, 'Input is not valid'
    assert len(image_stack.shape) == 3, 'Input array must contain 3 dimensions!'

    logging.debug('Determining pixel contiguity')
    # Create mask array with True for all pixels which are non-zero in all bands
    mask = image_stack.all(0)
    #mask = numexpr.evaluate('prod(image_stack, 0)') != 0 # ***This has the potential to overflow, roll back to array.all(0)***

    # The following is only valid for Landsat 5 images
    if spacecraft_id == 'LS5':
        logging.debug('Finding thermal edge anomalies')
        # Apply thermal edge anomalies
        struct = numpy.ones((7,7), dtype='bool')
        erode  = ndimage.binary_erosion(mask, structure=struct)

        dims = mask.shape
        th_anom = numpy.zeros(dims, dtype='bool').flatten()

        pix_3buff_mask = mask - erode
        pix_3buff_mask[pix_3buff_mask > 0] = 1
        edge = pix_3buff_mask == 1

        low_sat = image_stack[5,:,:] == 1
        low_sat_buff = ndimage.binary_dilation(low_sat, structure=struct)

        s = [[1,1,1],[1,1,1],[1,1,1]]
        low_sat, num_labels = ndimage.label(low_sat_buff, structure=s)

        labels = low_sat[edge]
        ulabels = numpy.unique(labels[labels > 0])

        # Testing a new method, more code but might be quicker
        #find_lab = numpy.in1d(low_sat.flatten(), ulabels)
        #th_anom |= find_lab

        # Histogram method, a lot faster
        mx = numpy.max(ulabels)
        h = histogram(low_sat.flatten(), min=0, max=mx, reverse_indices='ri')
        hist = h['histogram']
        ri = h['ri']

        for i in numpy.arange(ulabels.shape[0]):
            if hist[ulabels[i]] == 0:
                continue
            th_anom[ri[ri[ulabels[i]]:ri[ulabels[i]+1]]] = True

        th_anom = ~(th_anom.reshape(dims))
        mask &= th_anom

    return mask

def setContiguityBit(l1t_data, spacecraft_id, pq_const, pqaResult):
    mask = calc_contiguity_mask(l1t_data, spacecraft_id)
    bit_index = pq_const.contiguity
    pqaResult.set_mask(mask, bit_index)
