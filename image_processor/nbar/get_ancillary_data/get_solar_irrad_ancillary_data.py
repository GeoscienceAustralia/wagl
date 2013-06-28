'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''

import logging, os
from ULA3 import DataManager
from ULA3.dataset import SceneDataset
from ULA3.utils import log_multiline
from ULA3.ancillary.solar import get_solar_irrad
from ULA3.image_processor import ProcessorConfig

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    satellite = l1t_input_dataset.satellite
    assert satellite, 'Unable to find Satellite object'

    solar_irrad_file = os.path.join(CONFIG.DIR_SolarIrradianceLUT, satellite.SOLAR_IRRAD_FILE)

    solar_irrad_data = {} # Dict of irradiance values for each band
    for band, value in get_solar_irrad(l1t_input_dataset, solar_irrad_file).items():
        solar_irrad_data[band] = {
            'data_source': 'Solar Irradiance',
            'data_file': solar_irrad_file,
            'value': value
            }
    log_multiline(logger.info, solar_irrad_data, 'solar_irrad result', '\t')

    DATA.set_item('solar_irrad.dat', solar_irrad_data)
