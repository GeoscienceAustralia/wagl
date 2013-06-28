'''
calc_lat_long_grids.py

Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)

Calculates latitude & longitude grids, writes them to binary files and saves them in global DataManager object
'''

import logging
from ULA3 import DataManager
from ULA3.image_processor import ProcessorConfig
from ULA3.image_processor.utils import calc_lat_long_grids

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)
    calc_lat_long_grids(DataManager(), ProcessorConfig())
