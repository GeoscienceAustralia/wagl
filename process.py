#!/usr/bin/env python
'''
Created on Jun 20, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
Main script to be invoked from the command line
Replaces ULA nbar.py
'''

import logging, sys, os, traceback
import image_processor
import ULA3.image_processor as process_manager
from ULA3.image_processor import ProcessorConfig
from ULA3 import DataManager
from ULA3.utils import log_multiline

CONFIG = ProcessorConfig()
DATA = DataManager()

#if __name__ == '__main__':
# Set top logging level settings explicitly
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.DEBUG)
console_formatter = logging.Formatter('%(name)s : %(message)s')
console_handler.setFormatter(console_formatter)

root_logger = logging.getLogger('root')
if not root_logger.level:
    root_logger.setLevel(logging.INFO)
    root_logger.addHandler(console_handler)

root_logger.info('process started with args %s', CONFIG._args)

try:
    assert CONFIG.input, 'Input path for scene(s) must be specified'

    # Subprocess list specified (in order) by:
    #    1. Individual subprocess
    #    2. Custom subprocess list file
    #    3. Specified process level (e.g. 'nbar', 'pqa')
    #    4. Default (nbar) - uses __init__.txt
    subprocess = (process_manager.get_process_list_from_text(CONFIG.subprocess) or
        process_manager.get_process_list_from_file(CONFIG.subprocess_file) or
        process_manager.get_process_list_from_file(
            os.path.join(CONFIG.code_root, 'image_processor', CONFIG.process_level + '.txt')) or
        [])

    image_processor.process(subprocess, CONFIG.resume)
    root_logger.info('process completed successfully')

    if CONFIG.debug:
        DATA.save_all() # Save all in-memory data for debugging

    sys.exit(0)

except (Exception), e:
    root_logger.info('process failed')
    log_multiline(root_logger.error, traceback.format_exc(), 'process failed: ' + e.message, '\t')
    DATA.save_all() # Save all in-memory data for debugging or resuming
    sys.exit(1)

