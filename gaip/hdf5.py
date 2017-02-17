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


def attach_image_attributes(dataset, attrs={}):
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

    for key in attrs:
        dataset.attrs[key] = attrs[key]


def attach_table_attributes(dataset, title='Table', attrs={}):
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

    for key in attrs:
        dataset.attrs[key] = attrs[key]


def create_image_dataset(fid, dataset_name, shape, dtype, compression='lzf',
                         shuffle=True, chunks=(512, 512), attrs={}):
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

    for key in attrs:
        dset.attrs[key] = attrs[key]

    return dset


def write_h5_image(data, dset_name, group, kwargs, attrs={}):

    dset = group.create_dataset(dset_name, data=data, **kwargs)

    minv = data.min()
    maxv = data.max()

    dset.attrs['CLASS'] = 'IMAGE'
    dset.attrs['IMAGE_VERSION'] = '1.2'
    dset.attrs['DISPLAY_ORIGIN'] = 'UL'
    dset.attrs['IMAGE_MINMAXRANGE'] = [minv, maxv]

    for key in attrs:
        dset.attrs[key] = attrs[key]


def write_h5_table(data, dset_name, fid, compression='lzf', shuffle=True,
                   chunks=True, title='Table', attrs={}):

    compression_opts = None
    if compression == 'mafisc':
        compression = 32002
        shuffle = False
        compression_opts = (1, 0)

    dset = fid.create_dataset(dset_name, data=data, chunks=chunks,
                              shuffle=shuffle, compression=compression,
                              compression_opts=compression_opts)

    dset.attrs['CLASS'] = 'TABLE'
    dset.attrs['VERSION'] = '0.2'
    dset.attrs['TITLE'] = title

    # column names
    col_fmt = 'FIELD_{}_NAME'
    columns = data.dtype.names
    for i, col in enumerate(columns):
        dset.attrs[col_fmt.format(i)] = col

    for key in attrs:
        dset.attrs[key] = attrs[key]

    fid.flush()


def write_dataframe(df, dset_name, group, compression='lzf', title='Table',
                    attrs={}):
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


def attach_attributes(dataset, attrs):
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
    for key in attrs:
        dset.attrs[key] = attrs[key]

    return
