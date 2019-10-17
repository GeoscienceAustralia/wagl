

# ===============================================================================
# Copyright 2015 Geoscience Australia
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ===============================================================================

from __future__ import absolute_import, print_function
import numpy


def generate_tiles(samples, lines, xtile=None, ytile=None):
    """
    Generates a list of tile indices for a 2D array.

    :param samples:
        An integer expressing the total number of samples in an array.

    :param lines:
        An integer expressing the total number of lines in an array.

    :param xtile:
        (Optional) The desired size of the tile in the x-direction.
        Default is all samples

    :param ytile:
        (Optional) The desired size of the tile in the y-direction.
        Default is min(100, lines) lines.

    :return:
        Each tuple in the generator contains
        ((ystart,yend),(xstart,xend)).

    Example:

        >>> from wagl.tiling import generate_tiles
        >>> tiles = generate_tiles(8624, 7567, xtile=1000, ytile=400)
        >>> for tile in tiles:
        >>>     # A rasterio dataset
        >>>     subset = rio_ds.read([4, 3, 2], window=tile)
        >>>     # Or simply move the tile window across an array
        >>>     subset = array[tile] # 2D
        >>>     subset = array[:,tile[0],tile[1]] # 3D
    """
    def create_tiles(samples, lines, xstart, ystart):
        """
        Creates a generator object for the tiles.
        """
        for ystep in ystart:
            if ystep + ytile < lines:
                yend = ystep + ytile
            else:
                yend = lines
            for xstep in xstart:
                if xstep + xtile < samples:
                    xend = xstep + xtile
                else:
                    xend = samples
                yield (slice(ystep, yend), slice(xstep, xend))

    # check for default or out of bounds
    if xtile is None or xtile < 0:
        xtile = samples
    if ytile is None or ytile < 0:
        ytile = min(100, lines)

    xstart = numpy.arange(0, samples, xtile)
    ystart = numpy.arange(0, lines, ytile)

    tiles = create_tiles(samples, lines, xstart, ystart)

    return tiles


def scatter(iterable, n):
    """
    Evenly scatters an interable by `n` blocks.
    Sourced from:
    http://stackoverflow.com/questions/2130016/splitting-a-list-of-arbitrary-size-into-only-roughly-n-equal-parts

    :param iterable:
        An iterable or preferably a 1D list or array.

    :param n:
        An integer indicating how many blocks to create.

    :return:
        A `list` consisting of `n` blocks of roughly equal size, each
        containing elements from `iterable`.
    """

    q, r = len(iterable) // n, len(iterable) % n
    res = (iterable[i * q + min(i, r):(i + 1) * q + min(i + 1, r)]
           for i in range(n))
    return list(res)
