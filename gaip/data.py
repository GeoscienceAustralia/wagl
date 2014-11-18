import numpy as np
import rasterio
import os

from os.path import join as pjoin

def data(acq):
    dirname = acq.dir_name
    filename = acq.file_name
    with rasterio.open(pjoin(dirname, filename), 'r') as fo:
        return fo.read_band(1)
