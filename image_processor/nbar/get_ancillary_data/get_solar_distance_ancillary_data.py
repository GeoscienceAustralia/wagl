'''
Created on 10/07/2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''

import logging, os
from ULA3 import DataManager
from ULA3.dataset import SceneDataset
from ULA3.utils import log_multiline
from ULA3.ancillary.solar import get_solar_distance
from ULA3.image_processor import ProcessorConfig

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)
    solar_dist_file = os.path.join(CONFIG.DIR_EarthSun_LUT, "EarthSun_distanceLUT.txt")
    solar_dist_data = {
        'data_source': 'Solar distance',
        'data_file': solar_dist_file,
        'value': get_solar_distance(solar_dist_file, l1t_input_dataset.DOY)
        }
    log_multiline(logger.info, solar_dist_data, 'solar_dist result', '\t')

    DATA.set_item('solar_dist.dat', solar_dist_data)
