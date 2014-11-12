'''
Created on 08/01/2013

@author: u76345
'''
import os
from glob import glob
import logging
from scipy import ndimage
import numexpr
import fmask_cloud_masking as _fmask

def majority_filter(array, iterations=1):
    weights_array = [[1,1,1],[1,1,1],[1,1,1]]
    for i in range(iterations):
        array = ndimage.convolve(array, weights_array)
        array = numexpr.evaluate("array > 4")
    return array.astype('uint8')

def FMaskCloudMask(mtl, null_mask=None, cloud_prob=None, wclr_max=None, sat_tag=None, aux_data={}):
    Lnum=int(sat_tag[-1:])
    zen,azi,ptm,Temp,t_templ,t_temph,WT,Snow,fmask_byte,Shadow,dim,ul,resolu,zc, \
        geoT,prj = _fmask.plcloud(filename=mtl, mask=null_mask, num_Lst=Lnum, aux_data=aux_data)

    fmask_byte = fmask_byte == 1 # Convert to bool, True = Cloud, False not Cloud
    # Use a majority filter to fill holes, 2 iterations works well to smoothe things over
    fmask_byte = majority_filter(array=fmask_byte, iterations=2)

    return (fmask_byte != 1).astype('bool') # Invert to a 'land mask'

def calc_fmask_cloud_mask(l1t_data, l1t_sd, pq_const, contiguity_mask, aux_data={}):
