from __future__ import absolute_import
import os
from glob import glob
import logging
from scipy import ndimage
import numexpr
from . import fmask_cloud_masking as _fmask
from gaip.acca_cloud_masking import majority_filter

def fmask_cloud_mask(mtl, null_mask=None, cloud_prob=None, wclr_max=None,
                   sat_tag=None, aux_data={}):
    Lnum=int(sat_tag[-1:])
    (zen, azi, ptm, Temp, t_templ,
     t_temph, WT, Snow, fmask_byte,
     Shadow, dim, ul, resolu, zc,
     geoT, prj) = _fmask.plcloud(filename=mtl, mask=null_mask, num_Lst=Lnum,
                                 aux_data=aux_data)

    # Convert to bool, True = Cloud, False not Cloud
    fmask_byte = fmask_byte == 1
    # Use a majority filter to fill holes, 2 iterations works well to smoothe
    # things over
    fmask_byte = majority_filter(array=fmask_byte, iterations=2)

    return ~fmask_byte # Invert
