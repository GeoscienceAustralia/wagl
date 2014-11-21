"""
Data access functions.
"""

import numpy as np
import rasterio
import os
import gaip

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

def data_and_box(acq):
    """Return a tuple comprising the `numpy.array` containing the data of
    the acquisition `acq` together with the associated GriddedGeoBox describing
    the data extent. 
    The parameter `acq` should behave like a `gaip.Acquisition` object."""
    dirname = acq.dir_name
    filename = acq.file_name
    with rasterio.open(pjoin(dirname, filename), 'r') as fo:
        box = gaip.GriddedGeoBox.from_dataset(fo)
        return (fo.read_band(1), box)

def gridded_geo_box(acq):
    """Return a GriddedGeoBox instance representing the spatial extent and 
    grid associated with the acquisition `acq`. 
    The parameter `acq` should behave like a `gaip.Acquisition` object."""
    dirname = acq.dir_name
    filename = acq.file_name
    with rasterio.open(pjoin(dirname, filename), 'r') as fo:
        return gaip.GriddedGeoBox.from_dataset(fo)

def vstack_data(acqs):
    # determine result shape (all acquistions must have the same shape)
    stack_shape = (len(acqs), acqs[0].width, acqs[0].height)
    stack = np.empty(stack_shape)

    # determine data type by reading a byte

    a = acqs[0].data()

    # create the result array, setting datatype based on source type

    stack = np.empty(stack_shape, type(a[0][0]))
    
    # read aquisitions into it

    for i in range(0, stack_shape[0]):
        print acqs[i].dir_name, acqs[i].file_name
        stack[i] = acqs[i].data()

    return stack
    

