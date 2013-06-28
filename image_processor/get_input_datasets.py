#! /usr/bin/env python
'''
get_input_datasets.py

Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)

Instantiates SceneDataset objects for each input dataset and saves a reference in the global
DataManager object
'''
import logging, os, re
from ULA3 import DataManager
from ULA3.dataset import SceneDataset
from ULA3.image_processor import ProcessorConfig
from ULA3.utils import log_multiline

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    for _level in CONFIG.input.keys():

        input_re_pattern = getattr(CONFIG, '%s_RE_PATTERN' % _level.upper())
        logger.debug('CONFIG.%s_RE_PATTERN = %s', _level.upper(), input_re_pattern)
        # Warn if non-standard name
        if not re.search(input_re_pattern, CONFIG.input[_level]['name']):
            logger.warning('Non-standard input pathname: %s', CONFIG.input[_level]['path'])

        _input_path = os.path.abspath(CONFIG.input[_level]['path'])
        _input_dataset = SceneDataset(_input_path)
        assert _input_dataset, 'Unable to open scene ' + CONFIG.input[_level]['path']

        log_multiline(logger.info, _input_dataset.__dict__, 'SceneDataset object for %s' % _input_dataset._pathname, '\t')
        if CONFIG.debug:
            log_multiline(logger.info, _input_dataset._metadata._metadata_dict, 'Full metadata %s' % _input_dataset._pathname, '\t')
            log_multiline(logger.info, _input_dataset.satellite.__dict__, 'Satellite object for %s' % _input_dataset._pathname, '\t')

        # Check whether processing level is valid
        logger.info('Input dataset processor level = %s', _input_dataset.processor_level)
        if _level == 'l1t':
            assert _input_dataset.processor_level.lower() in ['l1t', 'ortho'], 'Dataset is ' + _input_dataset.processor_level
        else:
            assert _input_dataset.processor_level.lower() == _level, 'Dataset is ' + _input_dataset.processor_level

        DATA.set_item(_input_path, _input_dataset)


