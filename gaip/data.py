"""
Data access functions.
"""

import numpy as np
import rasterio
import os

from os.path import join as pjoin


def get_pixel(filename, lonlat, band=1):
    """Return a pixel from `filename` at the longitude and latitude given
    by the tuple `lonlat`. Optionally, the `band` can be specified."""
    with rasterio.open(filename) as src:
        x, y = [int(v) for v in ~src.affine * lonlat]
        return src.read_band(band, window=((y, y + 1), (x, x + 1))).flat[0]


def data(acq):
    """Return a `numpy.array` containing the data of the acquisition `acq`. 
    The parameter `acq` should behave like a `gaip.Acquisition` object."""
    dirname = acq.dir_name
    filename = acq.file_name
    with rasterio.open(pjoin(dirname, filename), 'r') as fo:
        return fo.read_band(1)
