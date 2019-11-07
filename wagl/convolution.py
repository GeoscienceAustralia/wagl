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
import pyfftw
import numexpr

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

    return slice(start_idx, end_idx)


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

    if outds is not None:
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

    # create compressed, chunk shape doesn't matter too much at this point
    # as when we need the data, the whole lot is read at once
    ds = out_group.create_dataset(dname, shape=pad.dims, compression='lzf',
                                  shuffle=True, dtype=dtype)
    chunks = ds.chunks

    buffered = numpy.pad(data, pad.width, mode=mode)

    for tile in generate_tiles(pad.dims[1], pad.dims[0], chunks[1], chunks[0]):
        ds[tile] = buffered[tile].astype(dtype)

    return ds


def conv_fft(data, kernel, out_group):
    """
    Convolve an array with a kernel via fourier.
    We are using are own version here as we can minimise the
    memory use as best we can. The version from astropy while
    convienient, was very memory hungry.
    In this function, we reuse arrays and storing temporary results
    on disk as much as possible.

    :param data:
        A 2D numpy.ndarray or h5py.Dataset that will be convolved with
        the kernel.

    :param kernel:
        A 2D numpy.ndarray or h5py.Dataset of the same shape/dimensions
        as data.

    :param out_group:
        A h5py.Group object to be used for storing intermediate results
        on disk.

    :return:
        A 2D numpy.ndarray of type float64 containin the result of the
        convolution process.
    """
    dims = data.shape
    indata = pyfftw.empty_aligned(dims, dtype='complex')
    outdata = pyfftw.empty_aligned(dims, dtype='complex')

    # define forward and reverse transforms
    flags = ['FFTW_DESTROY_INPUT', 'FFTW_MEASURE']
    fft_object = pyfftw.FFTW(indata, outdata, axes=(0, 1), flags=flags)
    ifft_object = pyfftw.FFTW(indata, outdata, axes=(0, 1),
                              direction='FFTW_BACKWARD', flags=flags)

    # read data into the input
    indata[:] = data[:]

    # forward transform
    result = fft_object()

    # temporary storage to hold the result
    ft_data_ds = out_group.create_dataset('data-ft', data=result, shuffle=True,
                                          compression='lzf')
    chunks = ft_data_ds.chunks

    # resuse existing input and output arrays for the kernel
    indata[:] = numpy.fft.ifftshift(kernel[:])
    result = fft_object()

    # temporary storage to hold the result
    ft_kern_ds = out_group.create_dataset('kern-ft', data=result, shuffle=True,
                                          compression='lzf', chunks=chunks)

    # apply convolution; resusing the input array to store the result
    for tile in generate_tiles(dims[1], dims[0], chunks[1], chunks[0]):
        numexpr.evaluate(
            "a * b",
            local_dict={
                'a': ft_data_ds[tile],
                'b': ft_kern_ds[tile]
            },
            out=indata[tile]
        )

    # now for the inverse
    result = ifft_object()

    return result.real


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
        The data will be padded using the 'reflect' method, whereas
        the kernel will be padded using the 'constant' method which
        uses a constant value of 0.
        The astropy library whilst convienient, proved very memory
        hungry. For example Landsat 8's panchromtatic band of
        dimensions (15461, 15241) would be killed on a laptop running
        8GB of memory when trying to convolve with a kernel of
        dimensions (77, 77). As such, a more memory efficient version
        has been implemented here.
    """
    # can we correctly apply row-length averages?
    row_idx = _sequential_valid_rows(data_mask)
    dims = (row_idx.stop - row_idx.start, data.shape[1])
    dtype = data.dtype
    fill = numpy.nan

    with tempfile.TemporaryDirectory(suffix='.tmp', prefix='convol-') as tmpd:
        with h5py.File(pjoin(tmpd, 'replace-nulls.h5'), 'w') as fid:

            # chunksize doesn't matter as we simply load all chunks
            outds = fid.create_dataset('fill-null', shape=dims, dtype=dtype,
                                       compression='lzf', shuffle=True)

            # fill nulls with run-length averages
            _fill_nulls(data[row_idx], data_mask, outds)

            # apply convolution
            if fourier:
                # determine required buffering/padding
                data_pad, kern_pad = _determine_pad(outds, kernel)

                # pad/buffer
                buffered_data = _pad(outds, data_pad, 'symmetric', fid,
                                     'buffered-data')
                buffered_kern = _pad(kernel, kern_pad, 'constant', fid,
                                     'buffered-kernel')

                # convolve
                convolved = conv_fft(buffered_data, buffered_kern, fid)

                # copy convolved into result taking into account both
                # the potential row subset and pad subset
                result = numpy.full(data.shape, fill_value=fill, dtype=dtype)
                result[row_idx] = convolved[data_pad.idx]
            else:
                # we could allocate the output array once, but we can conserve
                # memory by allocating right at the end where we need it
                result = numpy.full(data.shape, fill_value=fill, dtype=dtype)
                ndimage.convolve(outds, kernel, output=result[row_idx])

    return result
