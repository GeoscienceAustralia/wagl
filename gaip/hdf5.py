#!/usr/bin/env python

import numpy
import h5py

DEFAULT_IMAGE_CLASS = {'CLASS': 'IMAGE',
                       'IMAGE_VERSION': '1.2',
                       'DISPLAY_ORIGIN': 'UL'}

DEFAULT_TABLE_CLASS = {'CLASS': 'TABLE',
                       'VERSION': '0.2'}


def dataset_compression_kwargs(compression='lzf', shuffle=True,
                               chunks=(512, 512), compression_opts=None):
    """
    A helper function that aids in setting up the use of different
    compression filters for HDF5.
    The default is 'lzf', and other filters include:
    * 'mafisc' which uses LZMA and shuffle filters for heavy compression
    * 'bitshuffle' which uses LZ4 and shuffle filters for fast compression

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param shuffle:
        A `bool` indicating whether to apply the shuffle filter
        prior to compression. Can improve compression ratios.
        If set in conjuction with either 'lz4' or 'mafisc',
        then `shuffle` will be set to False as shuffle filters
        will be handled internally by the `bitshuffle` (for 'lz4')
        and `mafisc` filters.
        Default is `True`.

    :param chunks:
        A `tuple` containing the desired chunks sizes for each
        dimension axis of the dataset to be written to disk.
        Default is (512, 512).

    :param compression_opts:
        Finer grained control over compressio filters.
        Default is `None`.

    :return:
        A `dict` of key/value pairs of compression options for use
        with h5py's 'create_dataset' function.
    """
    if compression == 'mafisc':
        compression = 32002
        shuffle = False
        compression_opts = (1, 0)

    if compression == 'bitshuffle':
        compression = 32008
        shuffle = False
        compression_opts = (0, 2)

    kwargs = {'compression': compression,
              'shuffle': shuffle,
              'chunks': chunks,
              'compression_opts': compression_opts}

    return kwargs


def attach_image_attributes(dataset, attrs=None):
    """
    Attaches attributes to an HDF5 `IMAGE` Class dataset.

    :param dataset:
        A `NumPy` compound dataset.

    :param attrs:
        A `dict` containing any additional attributes to attach to
        the dataset.
    """
    for key in DEFAULT_IMAGE_CLASS:
        dataset.attrs[key] = DEFAULT_IMAGE_CLASS[key]

    if attrs is not None:
        for key in attrs:
            dataset.attrs[key] = attrs[key]


def attach_table_attributes(dataset, title='Table', attrs=None):
    """
    Attaches attributes to an HDF5 `TABLE` Class dataset.

    :param dataset:
        A `NumPy` compound dataset.

    :param title:
        A `string` containing the title of the TABLE dataset.
        Default is `Table`.

    :param attrs:
        A `dict` containing any additional attributes to attach to
        the dataset.
    """
    for key in DEFAULT_TABLE_CLASS:
        dataset.attrs[key] = DEFAULT_TABLE_CLASS[key]

    dataset.attrs['TITLE'] = title

    # column names
    col_fmt = 'FIELD_{}_NAME'
    columns = dataset.dtype.names
    for i, col in enumerate(columns):
        dataset.attrs[col_fmt.format(i)] = col

    if attrs is not None:
        for key in attrs:
            dataset.attrs[key] = attrs[key]


def create_image_dataset(fid, dataset_name, shape, dtype, compression='lzf',
                         shuffle=True, chunks=(512, 512), attrs=None):
    """
    Initialises a HDF5 dataset populated with some basic attributes
    that detail the Image specification.
    https://support.hdfgroup.org/HDF5/doc/ADGuide/ImageSpec.html

    :param fid:
        The file id opened for writing.

    :param dataset_name:
        The name used identify the dataset. The name can contain '/'
        used to identify the full path name. Any 'GROUP' names in the
        path will be created if they don't already exist.

    :param shape:
        The shape/dimensions of the dataset to create.

    :param compression:
        The compression algorithm to use. Default is 'lzf'.
        If 'mafisc' is used, then the `shuffle` filter will be turned
        off as alternate interal shuffling filter will be used instead.

    :shuffle:
        A `bool`, default is True, indicating whether or not to use
        HDF5's internal shuffle filter.

    :chunks:
        A `tuple` indicating the chunksize to use for each dimension
        of the dataset.

    :attrs:
        A `dict` by which the keys will be the attribute name, and the
        values will be the attribute value.
    """
    compression_opts = None
    if compression == 'mafisc':
        compression = 32002
        shuffle = False
        compression_opts = (1, 0)

    dset = fid.create_dataset(dataset_name, shape=shape, dtype=dtype,
                              compression=compression, chunks=chunks,
                              compression_opts=compression_opts,
                              shuffle=shuffle)

    if attrs is not None:
        for key in attrs:
            dset.attrs[key] = attrs[key]

    return dset


def write_h5_image(data, dset_name, group, attrs=None, **kwargs):

    dset = group.create_dataset(dset_name, data=data, **kwargs)

    minv = data.min()
    maxv = data.max()

    dset.attrs['CLASS'] = 'IMAGE'
    dset.attrs['IMAGE_VERSION'] = '1.2'
    dset.attrs['DISPLAY_ORIGIN'] = 'UL'
    dset.attrs['IMAGE_MINMAXRANGE'] = [minv, maxv]

    if attrs is not None:
        for key in attrs:
            dset.attrs[key] = attrs[key]


def write_h5_table(data, dset_name, group, compression='lzf', title='Table',
                   attrs=None):
    """
    Writes a `NumPy` structured array to a HDF5 `compound` type
    dataset.

    :param data:
        The `NumPy` structured array containing the data to be
        written to disk.

    :param dset_name:
        A `str` containing the name and location of the dataset
        to write to.

    :param group:
        A h5py `Group` or `File` object from which to write the
        dataset to.

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param title:
        A `str` containing the title name of the `Table` dataset.
        Default is 'Table'.

    :param attrs:
        A `dict` of key, value items to be attached as attributes
        to the `Table` dataset.
    """
    kwargs = dataset_compression_kwargs(compression=compression, chunks=True)
    dset = group.create_dataset(dset_name, data=data, **kwargs)

    dset.attrs['CLASS'] = 'TABLE'
    dset.attrs['VERSION'] = '0.2'
    dset.attrs['TITLE'] = title

    # column names
    col_fmt = 'FIELD_{}_NAME'
    columns = data.dtype.names
    for i, col in enumerate(columns):
        dset.attrs[col_fmt.format(i)] = col

    if attrs is not None:
        for key in attrs:
            dset.attrs[key] = attrs[key]


def write_dataframe(df, dset_name, group, compression='lzf', title='Table',
                    attrs=None):
    """
    Converts a `pandas.DataFrame` to a HDF5 `Table`, stored
    internall as a compound datatype.

    :param df:
        A `pandas.DataFrame` object.

    :param dset_name:
        A `str` containing the name and location of the dataset
        to write to.

    :param group:
        A h5py `Group` or `File` object from which to write the
        dataset to.

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param title:
        A `str` containing the title name of the `Table` dataset.
        Default is 'Table'.

    :param attrs:
        A `dict` of key, value items to be attached as attributes
        to the `Table` dataset.
    """
    # check for object types, for now write fixed length strings
    dtype = []
    for i, val in enumerate(df.dtypes):
        if val == object:
            str_sz = len(df[df.columns[i]].max())
            dtype.append((df.columns[i], '|S{}'.format(str_sz)))
        else:
            dtype.append((df.columns[i], val))

    dtype = numpy.dtype(dtype)

    kwargs = dataset_compression_kwargs(compression=compression, chunks=True)
    kwargs['shape'] = df.shape[0]
    kwargs['dtype'] = dtype
    dset = group.create_dataset(dset_name, **kwargs)

    for col in df.columns:
        dset[col] = df[col].values

    attach_table_attributes(dset, title=title, attrs=attrs)


def attach_attributes(dataset, attrs=None):
    """
    A small utility for attaching attributes to a h5py `Dataset` or
    `Group` object.

    :param dataset:
        The h5py `Dataset` or `Group` object on which to attach
        attributes to.

    :param attrs:
        A `dict` of key, value pairs used to attach the atrrbutes
        onto the h5py `Dataset` or `Group` object.

    :return:
        None.
    """
    if attrs is not None:
        for key in attrs:
            dataset.attrs[key] = attrs[key]

    return


def create_external_link(fname, dataset_path, out_fname, new_dataset_path):
    """
    Creates an external link of `fname:dataset_path` into
    `out_fname:new_dataset_path`.

    :param fname:
        The filename of the file containing the dataset that is to
        be linked to.
        `type:str`

    :param dataset_path:
        A `str` containg the path to the dataset contained within `fname`.

    :param out_fname:
        A `str` for the output filename that will contain the link to
        `fname:dataset_path`.

    :param new_dataset_path:
        A `str` containing the dataset path within `out_fname` that will
        link to `fname:dataset_path`.
    """
    with h5py.File(out_fname) as fid:
        fid[new_dataset_path] = h5py.ExternalLink(fname, dataset_path)

    return
