#! /usr/bin/env python
'''
image_reader.py

Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)

Creates a multi-band stack for each input dataset and saves a reference to each one
in the global DataManager object
'''
import logging
from ULA3.dataset import SceneDataset
from ULA3.common.pqa_result import PQAResult
from ULA3.image_processor import ProcessorConfig
from ULA3 import DataManager

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    nbar_input_dataset = DATA.get_item(CONFIG.input['nbar']['path'], SceneDataset)
    assert nbar_input_dataset, 'Unable to retrieve SceneDataset object for NBAR input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', nbar_input_dataset.pathname)

    # SceneDataset.ReadAsArray() only returns reflective & thermal bands (no panchromatic)
    l1t_stack = l1t_input_dataset.ReadAsArray()
    assert l1t_stack.shape[0] == len(l1t_input_dataset.satellite.BAND_TYPES['REFLECTIVE']) + len(l1t_input_dataset.satellite.BAND_TYPES['THERMAL']), "Unexpected number of L1T bands"

    nbar_stack = nbar_input_dataset.ReadAsArray()
    assert nbar_stack.shape[0] == len(l1t_input_dataset.satellite.BAND_TYPES['REFLECTIVE']), "Unexpected number of nbar bands"

    result = PQAResult(image_array=l1t_stack, dtype=CONFIG.pqa_dtype)

    DATA.set_item('l1t_stack', l1t_stack)
    DATA.set_item('nbar_stack', nbar_stack)
    DATA.set_item('result.tif', result)
