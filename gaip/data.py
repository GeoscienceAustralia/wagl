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

def stack_data(acqs_list, filter=(lambda acq: True)):
    """
    Given a list of acquisitions, apply the supplied filter to select the
    desired acquisitions and return the data from each acquisition
    collected in a 3D numpy array (first index is the acquisition number).

    :param acqs_list:
        The list of acquisitions to consider

    :param filter:
        A function that takes a single acquisition and returns True if the
        acquisition is to be selected for inclusion in the output

    :return:
        A tuple containing the list of selected acquisitions (possibly empty)
        and a 3D numpy array (or None) containing the corresponding
        acquisition data.
    """

    # get the subset of acquisitions required

    acqs = [acq for acq in acqs_list if filter(acq)]
    if len(acqs) == 0:
       return acqs, None

    # determine data type by reading the first band

    a = acqs[0].data()

    # create the result array, setting datatype based on source type

    stack_shape = (len(acqs), acqs[0].width, acqs[0].height)
    stack = np.empty(stack_shape, type(a[0][0]))
    stack[0] = a

    # read remaining aquisitions into it

    for i in range(1, stack_shape[0]):
        stack[i] = acqs[i].data()

    return acqs, stack

