#!/usr/bin/env python

"""
Helper classes for accessing some of the compression filters
supported by HDF5. See https://support.hdfgroup.org/services/contributions.html
for more details on the compression filters supported here, as well as
others that could be implemented here.

This method of compression filter access is based around HDF5's
dynamically loaded filters see
https://support.hdfgroup.org/HDF5/doc/Advanced/DynamicallyLoadedFilters/
https://support.hdfgroup.org/HDF5/doc/Advanced/DynamicallyLoadedFilters/HDF5DynamicallyLoadedFilters.pdf
for more details.

The class implementation utilises python-attrs
https://github.com/python-attrs/attrs for it's simplicity in
writing classes that can make use of default values after initialisation,
as well as the validators to help ensure that incorrect user inputs are
captured prior to dataset creation.
Also, all sub-classes of H5CompressionConfig are implemented as frozen instances
(immuteable). This is purely so each instance can act as a configuration,
such that a different setting will be a different configuration.
See http://www.attrs.org/en/stable/examples.html#immutability for more details.
"""

from enum import Enum, IntEnum
import attr


class H5CompressionFilter(IntEnum):

    """
    The available HDF5 compression filters.

    :example:
        >>> import numpy
        >>> import h5py
        >>> from wagl.hdf5 import H5CompressionFilter
        >>> compressor = H5CompressionFilter.LZF
        >>> config = compressor.config(chunks=(200, 300), shuffle=False)
        >>> kwargs = config.dataset_compression_kwargs()
        >>> with h5py.File('tmp.h5', 'w') as fid:
        >>>     data = numpy.random.randint(0, 256, (2000, 3000))
        >>>     fid.create_dataset('random-ints', data=data, **kwargs)
    """

    BLOSC_LZ = 0
    BLOSC_LZ4 = 1
    BLOSC_LZ4HC = 2
    BLOSC_SNAPPY = 3
    BLOSC_ZLIB = 4
    BLOSC_ZSTANDARD = 5
    LZF = 6
    GZIP = 7
    BITSHUFFLE = 8
    MAFISC = 9
    ZSTANDARD = 10

    def config(self, **kwargs):
        """
        Return the appropriate compression filter configuration
        class.
        """
        switch = {H5CompressionFilter.LZF: H5lzf,
                  H5CompressionFilter.GZIP: H5gzip,
                  H5CompressionFilter.BITSHUFFLE: H5bitshuffle,
                  H5CompressionFilter.MAFISC: H5mafisc,
                  H5CompressionFilter.ZSTANDARD: H5zstandard,
                  H5CompressionFilter.BLOSC_LZ: H5blosc,
                  H5CompressionFilter.BLOSC_LZ4: H5blosc,
                  H5CompressionFilter.BLOSC_LZ4HC: H5blosc,
                  H5CompressionFilter.BLOSC_SNAPPY: H5blosc,
                  H5CompressionFilter.BLOSC_ZLIB: H5blosc,
                  H5CompressionFilter.BLOSC_ZSTANDARD: H5blosc}

        return switch.get(self)(compression_filter=self, **kwargs)
    

class BloscCompression(IntEnum):

    """
    Comression filter ids for the blosc metafamily of compressors.
    """

    BLOSC_LZ = 0
    BLOSC_LZ4 = 1
    BLOSC_LZ4HC = 2
    BLOSC_SNAPPY = 3
    BLOSC_ZLIB = 4
    BLOSC_ZSTANDARD = 5


class BloscShuffle(IntEnum):

    """
    The shuffle filters available to the blosc metafamily of
    compressors.
    """

    NO_SHUFFLE = 0
    SHUFFLE = 1
    BITSHUFFLE = 2


@attr.s()
class H5CompressionConfig(object):

    """
    Base class for defining HDF5 compression filters.
    Users that know what they're doing can customise the
    dataset compression filter configuration here.
    But mostly, this class acts as a way for sub-classes to
    access the dataset_compression_kwargs method.
    """

    compression_filter = attr.ib()
    compression = attr.ib()
    compression_opts = attr.ib()
    chunks = attr.ib(default=True)
    shuffle = attr.ib(default=True)

    def dataset_compression_kwargs(self):
        """
        Return a dict with the required keywords for h5py's
        create_dataset() function.
        """
        include = (attr.fields(self.__class__).compression,
                   attr.fields(self.__class__).compression_opts,
                   attr.fields(self.__class__).chunks,
                   attr.fields(self.__class__).shuffle)
        kwargs = attr.asdict(self, filter=attr.filters.include(*include),
                             retain_collection_types=True)
        return kwargs


@attr.s(frozen=True)
class H5lzf(H5CompressionConfig):

    """
    Provides access to the LZF compression filter in-built with
    h5py, and optionally a shuffle filter.
    This filter requires that compression_opts be None.
    """

    compression_filter = attr.ib(default=H5CompressionFilter.LZF)
    compression = attr.ib(default='lzf')
    compression_opts = attr.ib(default=None)

    @compression_opts.validator
    def check_compression_opts(self, attribute, value):
        """
        Validate the compression_opts parameter.
        """
        if value:
            raise ValueError("compression_opts must be None")


@attr.s(frozen=True)
class H5gzip(H5CompressionConfig):

    """
    Provides access to the gzip/zlib/deflate compression filter.
    """

    compression_filter = attr.ib(default=H5CompressionFilter.GZIP)
    compression: str = attr.ib(default='gzip')
    aggression: int = attr.ib(default=4)
    compression_opts: int = attr.ib()

    @compression_opts.default
    def c_opts_default(self):
        """
        Set the default value for the compression_opts parameter.
        """
        return self.aggression

    @aggression.validator
    def check_aggression(self, attribute, value):
        """
        Validate the aggression parameter.
        """
        if value < 0 or value > 9:
            msg = "compression_opts must be in the interval [0, 9]"
            raise ValueError(msg)


@attr.s(frozen=True)
class H5zstandard(H5CompressionConfig):

    """
    See https://github.com/aparamon/HDF5Plugin-Zstandard for more
    details.
    Accesses the ZStandard compression filter, and optionally
    apply a shuffle filter.
    """

    compression_filter = attr.ib(default=H5CompressionFilter.ZSTANDARD)
    compression: int = attr.ib(default=32015)
    aggression: int = attr.ib(default=6)
    compression_opts = attr.ib()

    @compression_opts.default
    def c_opts_default(self):
        """
        Set the default value for the compression_opts parameter.
        """
        return (self.aggression,)

    @aggression.validator
    def check_aggression(self, attribute, value):
        """
        Validate the aggression parameter.
        """
        if value < 0 or value > 22:
            msg = "aggression must be in the interval [0, 22]"
            raise ValueError(msg)


@attr.s(frozen=True)
class H5bitshuffle(H5CompressionConfig):

    """
    See https://github.com/kiyo-masui/bitshuffle for more details.
    Uses a bitshuffle filter on top of the LZ4 compression filter.
    It can also be applied to the LZF compression filter, but it has
    not been integrated here.
    This implementation will utilise the LZ4 compression, and the
    bitshuffle filter.
    """

    compression_filter = attr.ib(default=H5CompressionFilter.BITSHUFFLE)
    compression: int = attr.ib(default=32008)
    shuffle: bool = attr.ib(default=False)
    compression_opts = attr.ib((0, 2))

    @shuffle.validator
    def check_shuffle(self, attribute, value):
        """
        Validate the shuffle parameter.
        """
        if value:
            cfilter = self.compression_filter.name
            msg = "shuffle keyword must be false for {}".format(cfilter)
            raise ValueError(msg)


@attr.s(frozen=True)
class H5mafisc(H5bitshuffle):

    """
    See https://wr.informatik.uni-hamburg.de/research/projects/icomex/mafisc
    for more details.
    Essentially it uses the LZMA comression filter, combined with
    data shuffling filters similar to bitshuffle.
    The goal of mafisc is to achieve the best compression possible,
    regardless of the time it takes.
    """

    compression_filter = attr.ib(default=H5CompressionFilter.MAFISC)
    compression: int = attr.ib(default=32002)
    compression_opts = attr.ib((1, 0))


@attr.s(frozen=True)
class H5blosc(H5bitshuffle):

    """
    See https://github.com/Blosc/hdf5-blosc for more details.
    Provides access to the blosc metafamily of compressors.
    Available compression filters are given by BloscCompression,
    and the shuffle options are given by BloscShuffle.
    The results of filters such as ZStandard, lz4+bitshuffle, zlib
    accessed via blosc may not provide the same result as if accessing
    those filters directly.
    """

    compression_filter = attr.ib(default=H5CompressionFilter.BLOSC_LZ)
    shuffle_id = attr.ib(default=BloscShuffle.SHUFFLE,
                         validator=attr.validators.in_(BloscShuffle))
    compression: int = attr.ib(default=32001)
    aggression: int = attr.ib(default=4)
    compression_opts = attr.ib()

    @compression_opts.default
    def c_opts_default(self):
        """
        Set the default value for the compression_opts parameter.
        """
        blosc_id = self.compression_filter.value
        return (0, 0, 0, 0, self.aggression, self.shuffle_id.value, blosc_id)

    @compression_filter.validator
    def check_compression_filter(self, attribute, value):
        """
        Validate the compression_filter parameter.
        """
        if value < 0 or value > 5:
            msg = "compression_filter: '{}' not supported by blosc".format(value)
            raise ValueError(msg)

    @aggression.validator
    def check_aggression(self, attribute, value):
        """
        Validate the aggression parameter.
        """
        if value < 0 or value > 9:
            msg = "aggression must be in the interval [0, 9]"
            raise ValueError(msg)
