'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
Incorporates original code from ULA main.py and sol_grids.py written by Roger Edberg
'''

import logging
from ULA3 import DataManager
from ULA3.image_processor import ProcessorConfig
from ULA3.image_processor.utils import calc_solar_grids

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    """Generate all solar grids
    """
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    calc_solar_grids(DATA, CONFIG)
