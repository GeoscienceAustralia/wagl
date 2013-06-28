'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
Code adapted from NbarProcessor.py and main.py

Create composite preview thumbnail from RGB product files
'''
# As at 26/7/12, old code has been imported and partially cleaned up, but not tested

import logging
from ULA3 import DataManager
from ULA3.image_processor import ProcessorConfig
from ULA3.image_processor.utils import create_output_image

logger = logging.getLogger('root.' + __name__)

THUMBNAIL_WIDTH = 1024

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    res = create_output_image(DATA, CONFIG, THUMBNAIL_WIDTH, 'thumbnail')

    # Store thumbnail info for downstream metadata output
    try:
        DATA.set_item('thumbnail.dat', {
            'filename': res['browse_image_path'],
            'size': res['image_dim'],
            'rgb_bands': [band_file_number // 10 for band_file_number in res['l1t_input_dataset'].satellite.rgb_bands]})
    except Exception, e:
        logger.error('ERROR: thumbnail.create_thumbnail RAISED EXCEPTION (ignored)')
        logger.error('[ %s ]' % e)
