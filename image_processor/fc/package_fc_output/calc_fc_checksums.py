'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
Code adapted from NbarProcessor.py and main.py

Create composite preview thumbnail from RGB product files
'''
# As at 26/7/12, old code has been imported and partially cleaned up, but not tested

import logging
from ULA3 import DataManager
from ULA3.utils import generate_md5

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    DATA = DataManager()

    fc_temp_output = DATA.get_item('fc_temp_output.dat', str)
    assert fc_temp_output, 'Unable to retrieve fc_temp_output string'
    logger.debug('string for fc_temp_output retrieved')

    generate_md5(fc_temp_output)
