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
from osgeo import gdal
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
        >>>     ystart = int(tile[0][0])
        >>>     yend = int(tile[0][1])
        >>>     xstart = int(tile[1][0])
        >>>     xend = int(tile[1][1])
        >>>     xsize = int(xend - xstart)
        >>>     ysize = int(yend - ystart)
        >>>     # When used to read data from disk
        >>>     subset = gdal_indataset.ReadAsArray(xstart, ystart, xsize, ysize)
        >>>     # The same method can be used to write to disk.
        >>>     gdal_outdataset.WriteArray(array, xstart, ystart)
        >>>     # A rasterio dataset
        >>>     subset = rio_ds.read([4, 3, 2], window=tile)
        >>>     # Or simply move the tile window across an array
        >>>     subset = array[ystart:yend,xstart:xend] # 2D
        >>>     subset = array[:,ystart:yend,xstart:xend] # 3D
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
                yield ((ystep, yend), (xstep, xend))

    # check for default or out of bounds
    if xtile is None or xtile < 0:
        xtile = samples
    if ytile is None or ytile < 0:
        ytile = min(100, lines)

    xstart = numpy.arange(0, samples, xtile)
    ystart = numpy.arange(0, lines, ytile)

    tiles = create_tiles(samples, lines, xstart, ystart)

    return tiles


class TiledOutput:
    def __init__(
        self,
        out_fname,
        samples=None,
        lines=None,
        bands=1,
        geobox=None,
        fmt="ENVI",
        nodata=None,
        dtype=gdal.GDT_Byte,
    ):
        """
        A class to aid in data processing using a tiling scheme.
        The `TiledOutput` class takes care of writing each tile/chunk
        of data to disk.

        :param out_fname:
            A string containing the full filepath name used for
            creating the image on disk.

        :param samples:
            An integer indicating the number of samples/columns contained
            in the entire image.

        :param lines:
            An integer indicating the number of lines/rows contained in
            the entire image.

        :param bands:
            The number of bands contained in the image. Default is 1.

        :param geobox:
            An instance of a GriddedGeoBox object.

        :param fmt:
            A GDAL compliant file format for the output image.

        :param nodata:
            The no data value for the image.

        :param dtype:
            An integer indicating datatype for the output image.
            Default is gdal.GDT_Byte which corresponds to 1.

        :example:
            >>> a = numpy.random.randint(0, 256, (1000, 1000)).astype('uint8')
            >>> tiles = generate_tiles(a.shape[1], a.shape[0], 100, 100)
            (a - img).sum() == 0
            >>> len(tiles)
            100
            >>> outds = TiledOutput('test_tiled_output_2D', samples=a.shape[1], lines=a.shape[0])  # noqa: E501
            >>> outds.closed
            False
            >>> for tile in tiles:
            ...     ys, ye = tile[0]
            ...     xs, xe = tile[1]
            ...     subset = a[ys:ye, xs:xe]
            ...     outds.write_tile(subset, tile)
            ...
            >>> outds.close()
            >>> outds.closed
            True
            >>>
            >>> inds = gdal.Open('test_tiled_output_2D')
            >>> img = inds.ReadAsArray()
            >>> a.shape == img.shape
            True
            >>> (a - img).sum() == 0
            True
            >>> a = numpy.random.randint(0, 256, (10, 100, 100)).astype('uint8')
            outds.closed
            >>> tiles = generate_tiles(a.shape[2], a.shape[1], 10, 10)
            >>> outds = TiledOutput('test_tiled_output_3D', samples=a.shape[2], lines=a.shape[1], bands=a.shape[0])  # noqa: E501
            for tile in tiles:
                ys, ye = tile[0]
                xs, xe = tile[1]
                subset = a[:, ys:ye, xs:xe]
                outds.write_tile(subset, tile)
            >>> outds.closed
            False
            >>> for tile in tiles:
            ...     ys, ye = tile[0]
            ...     xs, xe = tile[1]
            ...     subset = a[:, ys:ye, xs:xe]
            ...     outds.write_tile(subset, tile)
            ...
            outds.close()
            outds.closed
            inds = gdal.Open('test_tiled_output_3D')
            img = inds.ReadAsArray()
            a.shape == img.shape
            (a - img).sum() == 0
            >>> outds.close()
            >>> outds.closed
            True
            >>> inds = gdal.Open('test_tiled_output_3D')
            >>> img = inds.ReadAsArray()
            >>> a.shape == img.shape
            True
            >>> (a - img).sum() == 0
            True
        """

        # Check we have the correct dimensions to create the file
        if samples is None or lines is None:
            msg = (
                "Samples and lines are required inputs! Samples: {ns} " "Lines: {nl}"
            ).format(ns=samples, nl=lines)
            raise TypeError(msg)

        driver = gdal.GetDriverByName(fmt)
        self.outds = driver.Create(out_fname, samples, lines, bands, dtype)

        self.nodata = nodata
        self.geobox = geobox
        self.bands = bands
        self.out_bands = None

        self._set_geobox()
        self._set_bands_lookup()
        self._set_nodata()
        self.closed = False

    def _set_geobox(self):
        """
        Assign the georeference information (if we have any) to the output
        image.
        """
        if self.geobox is not None:
            transform = self.geobox.affine.to_gdal()
            projection = bytes(self.geobox.crs.ExportToWkt())

            self.outds.SetGeoTransform(transform)
            self.outds.SetProjection(projection)

    def _set_bands_lookup(self):
        """
        Define a dictionary that points to each raster band object on disk.
        Used for writing the actual data and setting items such as no data
        values.
        """
        out_bands = {}
        for i in range(1, self.bands + 1):
            out_bands[i] = self.outds.GetRasterBand(i)

        self.out_bands = out_bands

    def _set_nodata(self):
        """
        If we have a no data value, then assign it to the output file.
        """
        if self.nodata is not None:
            for i in range(1, self.bands + 1):
                self.out_bands[i].SetNoDataValue(self.nodata)

    def write_tile(self, array, tile, raster_band=None):
        """
        Given an array and tile index in the form:

            ((ystart, yend), (xstart, xend))

        and optionally a raster band, write the current tile to disk.
        The input `array` can be 2D or 3D. For the 2D case it is assumed
        that `array` will be written to the 1st raster bands unless
        specified by the parameter `raster_band`. For the 3D case it is
        assumed that `array` will have the dimensions (bands, rows, cols),
        whereby the number of bands will written to disk seqentially.
        The `raster_band` parameter will be ignored if `array` is 3D.
        """

        dims = array.ndim
        if array.ndim not in [2, 3]:
            msg = (
                "Input array is not 2 or 3 dimensions. " "Array dimensions: {dims}"
            ).format(dims=dims)
            raise TypeError(msg)

        ystart = tile[0][0]
        xstart = tile[1][0]
        if dims == 3:
            for i in range(self.bands):
                band = i + 1
                self.out_bands[band].WriteArray(array[i], xstart, ystart)
                self.out_bands[band].FlushCache()
        else:
            band = 1 if raster_band is None else raster_band
            self.out_bands[band].WriteArray(array, xstart, ystart)
            self.out_bands[band].FlushCache()

    def close(self):
        """
        Close the output image and flush everything still in cache to disk.
        """

        for band in self.out_bands:
            self.out_bands[band] = None

        self.out_bands = None
        self.outds = None
        self.closed = True


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
    res = (iterable[i * q + min(i, r) : (i + 1) * q + min(i + 1, r)] for i in range(n))
    return list(res)
