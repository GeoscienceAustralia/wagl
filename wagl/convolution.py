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

from collections import namedtuple
from os.path import join as pjoin
import tempfile
import numpy
from scipy import ndimage
import h5py
from astropy.convolution import convolve_fft
import pyfftw

from wagl.tiling import generate_tiles


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


def _determine_pad(data, kernel):
    """
    Determines the required padding for a data array given a kernel.
    Ideally half the kernel size is added to each edge of the data array,
    plus and an additional amount to get to the nearest power of 2.
    The power of 2 enables faster fourier transforms to take place.

    :param data:
        A 2D numpy.ndarray or h5py.Dataset.

    :param kernel:
        A 2D numpy.ndarray or h5py.Dataset.

    :return:
        An instance of a namedtuple:

        Pad(dims, width, idx)

        where:

        dims -> array dimensions/shape
        width -> padding width
        idx -> index location of the original array
    """
    # define a container to hold specific return info
    Pad = namedtuple('Pad', 'dims width idx')

    # data and kernel y & x dimensions (shape)
    data_dims = data.shape
    kern_dims = kernel.shape

    # minimum pad required for kernel (half +1 for good measure)
    min_pad = [i // 2 + 1 for i in kern_dims]

    # additional pad to enable FFT (power of 2)
    out_dims = tuple([
        pyfftw.next_fast_len(data_dims[0] + min_pad[0]*2),
        pyfftw.next_fast_len(data_dims[1] + min_pad[1]*2)
    ])

    data_offset = (
        out_dims[0] - data_dims[0],
        out_dims[1] - data_dims[1]
    )
    kern_offset = (
        out_dims[0] - kern_dims[0],
        out_dims[1] - kern_dims[0]
    )

    # additional padding on each 2D edge
    d_pad = tuple([
        ((data_offset[0] + 1) // 2, data_offset[0] // 2),
        ((data_offset[1] + 1) // 2, data_offset[1] // 2),
    ])
    k_pad = tuple([
        ((kern_offset[0] + 1) // 2, kern_offset[0] // 2),
        ((kern_offset[1] + 1) // 2, kern_offset[1] // 2),
    ])

    # indices for original locations of data and kernel
    d_unpad_idx = (
        slice(d_pad[0][0], -d_pad[0][1]),
        slice(d_pad[1][0], -d_pad[1][1])
    )
    k_unpad_idx = (
        slice(k_pad[0][0], -k_pad[0][1]),
        slice(k_pad[1][0], -k_pad[1][1])
    )

    data_pad = Pad(dims=out_dims, width=d_pad, idx=d_unpad_idx)
    kern_pad = Pad(dims=out_dims, width=k_pad, idx=k_unpad_idx)

    return data_pad, kern_pad


def _pad(data, pad, mode, out_group, dname):
    """
    Pad/buffer an array according to the pad_width.
    The output dataset will be cast to type complex.
    The reason is that internal to astropy it will be recast as complex,
    so this is simply a way of conserving memory.

    :param data:
        A 2D numpy.ndarray or h5py.Dataset.

    :param pad:
        An instance of a namedtuple Pad(dims, width, idx).

    :param mode:
        See numpy.pad for mode options.

    :param out_group:
        A h5py.Group at which to write the h5py.Dataset.
        The padded data will be output as type complex.

    :param dname:
        A name to use for the creation of the h5py.Dataset.

    :return:
        A h5py.Dataset containing the padded array.
    """
    dtype = 'complex'
    ds = out_group.create_dataset(dname, shape=pad.dims, compression='lzf',
                                  shuffle=True, dtype=dtype)
    chunks = ds.chunks

    buffered = numpy.pad(data, pad.width, mode=mode)

    for tile in generate_tiles(pad.dims[1], pad.dims[0], chunks[1], chunks[0]):
        ds[tile] = buffered.astype(dtype)

    return ds


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
        We'll try to handle the required padding before passing to
        astropy.
        A temporary HDF5 file will be created to store intermediates
        in order to conserve memory.
        As we're expecting a normalised kernel, we set astropy's
        normalization_zero_tol equal to 2 and set nan_treatment to fill,
        in order to for astropy not to do any normalisation.
        One would think setting to False would mean it wouldn't do any
        normalisation.
        The data will be padded using the 'reflect' method, whereas
        the kernel will be padded using the 'constant' method which
        uses a constant value of 0.
    """
    # can we correctly apply row-length averages?
    start_idx, end_idx = _sequential_valid_rows(data_mask)
    row_idx = slice(start_idx, end_idx)

    with tempfile.TemporaryDirectory(suffix='.tmp', prefix='convol-') as tmpd:
        with h5py.File(pjoin(tmpd, 'replace-nulls.h5'), 'w') as fid:

            # chunksize doesn't matter as we simply load all chunks
            outds = fid.create_dataset('fill-null', shape=data.shape,
                                       dtype=data.dtype, compression='lzf',
                                       shuffle=True)

            # fill nulls with run-length averages
            _fill_nulls(data[:], data_mask, outds)

            # apply convolution
            if fourier:
                # determine required buffering/padding
                data_pad, kern_pad = _determine_pad(data[row_idx], kernel)

                # pad/buffer
                buffered_data = _pad(outds, data_pad, 'reflect', fid,
                                     'buffered-data')
                buffered_kern = _pad(kernel, kern_pad, 'constant', fid,
                                     'buffered-kernel')

                # convolve
                convolved = convolve_fft(buffered_data, buffered_kern,
                                         allow_huge=True, fft_pad=False,
                                         psf_pad=False, normalize_kernel=False,
                                         normalization_zero_tol=2,
                                         nan_treatment='fill')

                # copy convolved into result taking into account both
                # the potential row subset and pad subset
                result = numpy.full(data.shape, fill_value=numpy.nan,
                                    dtype='float32')
                result[row_idx] = convolved[data_pad.idx]
            else:
                result = numpy.full(data.shape, fill_value=numpy.nan,
                                    dtype='float32')
                ndimage.convolve(outds[row_idx], kernel,
                                 output=result[row_idx])

    return result
