'''
Created on 08/01/2013

@author: u76345
'''
import os
from glob import glob
import logging
from scipy import ndimage
import numexpr

from ULA3.dataset import SceneDataset
from ULA3.common.pqa_result import PQAResult
from ULA3.image_processor import ProcessorConfig
from ULA3.image_processor import constants
from ULA3 import DataManager
from ULA3.utils import dump_array
import fmask_cloud_masking as _fmask

def majority_filter(array, iterations=1):
    weights_array = [[1,1,1],[1,1,1],[1,1,1]]
    for i in range(iterations):
        array = ndimage.convolve(array, weights_array)
        array = numexpr.evaluate("array > 4")
    return array.astype('uint8')

def FMaskCloudMask(mtl, null_mask=None, cloud_prob=None, wclr_max=None):
    #cloud_prob = cloud_prob or CONFIG.pqa_param['fmask_cloudprob']
    #wclr_max = wclr_max or CONFIG.pqa_param['fmask_wclr_max']

    # Open the MTL file.
    # The original MATLAB code opens the file twice to retrieve the Landsat number.
    # It would be better to open it once and restructure the function parameters.
    data = {}
    fl = open(mtl,'r')
    file_lines = fl.readlines()
    for line in file_lines:
        values = line.split(' = ')
        if len(values) != 2:
            continue

        data[values[0].strip()] = values[1].strip().strip('"')

    fl.close()

    # Identify Landsat Number (Lnum = 4, 5 or 7)
    LID=data['SPACECRAFT_ID']
    Lnum=int(LID[len(LID)-1])

    try: # Need to change directory to ensure log files are written to the right place
        current_dir = os.curdir
        os.chdir(os.path.join(pqa_temp_output, 'scene01'))

        zen,azi,ptm,Temp,t_templ,t_temph,WT,Snow,fmask_byte,Shadow,dim,ul,resolu,zc, \
            geoT,prj = _fmask.plcloud(filename=mtl, mask=null_mask, num_Lst=Lnum)

    finally:
        os.chdir(current_dir)

    fmask_byte = fmask_byte == 1 # Convert to bool, True = Cloud, False not Cloud
    # Use a majority filter to fill holes, 2 iterations works well to smoothe things over
    fmask_byte = majority_filter(array=fmask_byte, iterations=2)

    return (fmask_byte != 1).astype('bool') # Invert to a 'land mask'

def calc_fmask_cloud_mask(l1t_data, l1t_sd, pq_const, contiguity_mask, aux_data={}):
