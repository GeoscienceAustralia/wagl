'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''

import logging
from ULA3 import DataManager
from ULA3.dataset import SceneDataset
from ULA3.utils import log_multiline
from ULA3.image_processor import ProcessorConfig
from ULA3.ancillary.elevation import get_elevation_data

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)
    CONFIG = ProcessorConfig()
    DATA = DataManager()
    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)
    elevation_data = get_elevation_data(l1t_input_dataset.lonlats['CENTRE'], CONFIG.DIR_DEM)
    log_multiline(logger.info, elevation_data, 'Elevation result', '\t')
    DATA.set_item('elevation.dat', elevation_data)
