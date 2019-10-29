#!/usr/bin/env python

"""
Module geared towards applying convolution specifically to enable
-----------------------------------------------------------------

atmospheric adjacency correction
--------------------------------

The module contains various methods for filtering, handling null
data, and padding, in order to avoid the extreme biased behaviour
that can occur during convolution.
"""

from os.path import join as pjoin
import tempfile
import numpy
from scipy import ndimage
import h5py
from astropy.convolution import convolve_fft


def _sequential_valid_rows(mask):
    """
    Check that the mask of null data contains sequential rows.
    i.e. non-sequential means an invalid row is surrounded by valid
    rows.
    """
    nrows = mask.shape[0]
    row_ids = numpy.arange(nrows)

    # identify any rows that are completely null
    invalid_rows_mask = numpy.all(mask, axis=1)

    valid_rows = row_ids[~invalid_rows_mask]
    start_idx = valid_rows[0]
    end_idx = valid_rows[-1] + 1

    # check that we are ascending by 1
    # i.e. an invalid row has valid rows before and after
    # TODO;
    # don't raise, simply return and we use another method to fill nulls
    sequential = numpy.all(numpy.diff(valid_rows) == 1)
    if not sequential:
        msg = "Rows with valid data are non-sequential."
        raise Exception(msg)

    return start_idx, end_idx


def _fill_nulls(data, mask, outds):
    """
    Calculate run-length averages and insert (inplace) at null pixels.
    In future this could also be an averaging kernel.
    We are assuming 2D arrays only (y, x).
    """
    # TODO: how to account for a row of data that is all null?
    # potentially use alternate methods if row length average doesn't satisfy
    for row in range(data.shape[0]):
        data[row][mask[row]] = numpy.mean(data[row][~mask[row]])

    outds.write_direct(data)


def convolve(data, kernel, data_mask, fourier=False):
    """
    Apply convolution.
    Internally uses temp files to store intermediates in order to conserve
    memory.

    :param data:
        A 2D numpy.array or h5py.Dataset containing the data that is
        to be convolved.

    :param kernel:
        The kernel to be used in convolving over the data array.

    :param fourier:
        A bool indicating whether to undertake convolution via fourier.
        Useful when dealing with large kernels.

    :return:
        A 2D numpy.array of type float32.

    :notes:
        The data input fills null pixels with an average value
        calculated by a row/run length average.
    """
    # can we correctly apply row-length averages?
    start_idx, end_idx = _sequential_valid_rows(data_mask)

    with tempfile.TemporaryDirectory(suffix='.tmp', prefix='convol-') as tmpd:
        with h5py.File(pjoin(tmpd, 'replace-nulls.h5', 'w')) as fid:

            # chunksize doesn't matter as we simply load all chunks
            outds = fid.create_dataset('fill-null', shape=data.shape,
                                       dtype=data.dtype, compression='lzf',
                                       shuffle=True)

            # fill nulls with run-length averages
            _fill_nulls(data[:], data_mask, outds)

            # apply convolution
            if fourier:
                # determine required buffering/padding
                # half kernel size (+1 for good measure :) )
                pad = [i // 2 + 1 for i in kernel.shape]

                # index to unpad the array
                unpad_idx = (slice(pad[0], -pad[0]), slice(pad[1], -pad[1]))

                # pad/buffer and convolve
                buffered = numpy.pad(outds[start_idx:end_idx], pad,
                                     mode='reflect')
                convolved = convolve_fft(buffered, kernel, allow_huge=True)

                # copy convolved into result taking into account both
                # the potential row subset and pad subset
                result = numpy.full(data.shape, fill_value=numpy.nan,
                                    dtype='float32')
                result[start_idx:end_idx] = convolved[unpad_idx]
            else:
                result = numpy.full(data.shape, fill_value=numpy.nan,
                                    dtype='float32')
                ndimage.convolve(outds[start_idx:end_idx], kernel,
                                 output=result[start_idx:end_idx])

    return result
