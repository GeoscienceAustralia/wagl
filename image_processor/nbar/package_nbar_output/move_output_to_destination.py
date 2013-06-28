'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
Code adapted from NbarProcessor.py and main.py

Create composite preview thumbnail from RGB product files
'''

import logging
from ULA3 import DataManager
from ULA3.image_processor.utils import move_outputs

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    DATA = DataManager()

    nbar_temp_output = DATA.get_item('nbar_temp_output.dat', str)
    assert nbar_temp_output, 'Unable to retrieve nbar_temp_output string'
    logger.debug('string for nbar_temp_output retrieved')

    nbar_output_path = DATA.get_item('nbar_output_path.dat', str)
    assert nbar_output_path, 'Unable to retrieve nbar_output_path string'
    logger.debug('string for nbar_output_path retrieved')

    move_outputs(nbar_temp_output, nbar_output_path)
