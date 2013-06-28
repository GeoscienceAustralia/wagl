'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''
import logging
from ULA3.image_processor import ProcessorConfig
from ULA3.utils import execute

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()

    if not CONFIG.DEBUG:
        execute('rm -rf %s/mod' % CONFIG.work_path)

