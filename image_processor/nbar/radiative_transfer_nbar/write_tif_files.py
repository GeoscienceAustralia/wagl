'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''

import logging, numpy
from datetime import datetime
from osgeo import gdal, gdalconst, osr
from ULA3 import DataManager
from ULA3.dataset import SceneDataset
from ULA3.utils import write_tif_file
from ULA3.image_processor import ProcessorConfig

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    ref_dtype = numpy.dtype('int16')

    # Change string argument to list if necessary
    if type(subprocess_list) == str:
        subprocess_list = subprocess_list.split(',')

    # Create sub-sublogger to allow interleaved multithreaded logs to be sorted out afterwards
    if len(subprocess_list) == 1:
        if type(subprocess_list[0]) == str:
            sublogger_name = logger.name + '.' + subprocess_list[0]
        elif type(subprocess_list[0]) == list:
            sublogger_name = logger.name + '.' + '.'.join(subprocess_list[0])
        else:
            sublogger_name = logger.name
        logger.debug('sublogger_name = %s', sublogger_name)
        sublogger = logging.getLogger(sublogger_name)
        sublogger.debug('%s.process(%s, %s) called', __name__, subprocess_list, resume)
    else:
        sublogger = logger

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    sublogger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    # List containing valid reflective band numbers as strings
    band_strings = [str(band_number) for band_number in l1t_input_dataset.bands('REFLECTIVE')]

    _subprocess_list = subprocess_list or band_strings

#    sublogger.debug( '_subprocess_list = %s', _subprocess_list)

    nbar_dataset_id = DATA.get_item('nbar_dataset_id.dat', str)
    assert nbar_dataset_id, 'Unable to retrieve nbar_dataset_id string'
    logger.debug('string for %s retrieved', nbar_dataset_id)

    nbar_temp_output = DATA.get_item('nbar_temp_output.dat', str)
    assert nbar_temp_output, 'Unable to retrieve nbar_temp_output string'
    logger.debug('string for %s retrieved', nbar_temp_output)

    for subprocess in _subprocess_list:
        if subprocess:
            if type(subprocess) == str:
                subprocess = subprocess.split('.')
            assert len(subprocess) == 1, 'Only band can be specified'
            if subprocess[0] in band_strings:
                write_tif_file(l1t_input_dataset, int(subprocess[0]), nbar_dataset_id, nbar_temp_output, ref_dtype, sublogger, CONFIG.work_path, CONFIG.debug)
            else:
                sublogger.warning('Ignoring invalid band: %s' % subprocess[0])

    # Set earliest creation time for metadata
    create_datetime = DATA.get_item('create_datetime', datetime) or datetime.utcnow()
    DATA.set_item('create_datetime', create_datetime)

    sublogger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

