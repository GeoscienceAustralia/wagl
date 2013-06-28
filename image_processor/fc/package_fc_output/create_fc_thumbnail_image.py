'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
Code adapted from NbarProcessor.py and main.py

Create composite preview thumbnail from RGB product files
'''

#TODO: Make this work for single multi-band TIF files as well as multiple single-band TIFs

import os, logging
from ULA3 import DataManager
from ULA3.image_processor import ProcessorConfig
#from ULA3.image_processor.utils import create_output_image
from ULA3.utils import create_output_image
from ULA3.fc.utils import create_dir

logger = logging.getLogger('root.' + __name__)

THUMBNAIL_WIDTH = 1024

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    fc_temp_output = DATA.get_item('fc_temp_output.dat', str)
    assert fc_temp_output, 'Unable to retrieve fc_temp_output string'
    logger.debug('string for fc_temp_output retrieved')

    fc_dataset_id = DATA.get_item('fc_dataset_id.dat', str)
    assert fc_dataset_id, 'Unable to retrieve fc_dataset_id string'
    logger.debug('string for fc_dataset_id retrieved')

    rgb_files = [os.path.join(fc_temp_output, 'scene01', fc_dataset_id + '_' + file_suffix + '.tif') for file_suffix in ['BS', 'PV', 'NPV']]
    thumbnail_path = os.path.join(fc_temp_output, fc_dataset_id + '.jpg')
    thumbnail_temp_dir = os.path.join(CONFIG.work_path, 'fc_thumbnail')
    create_dir(thumbnail_temp_dir)

    create_output_image(rgb_files,
                        thumbnail_path,
                        thumbnail_temp_dir,
                        1024)
