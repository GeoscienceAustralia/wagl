'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''

import logging
from datetime import datetime
from ULA3 import DataManager
from ULA3.dataset import SceneDataset
from ULA3.ancillary.aerosol import get_aerosol_value_for_region
from ULA3.image_processor import ProcessorConfig
from ULA3.utils import log_multiline
from ULA3.meta import print_call

logger = logging.getLogger('root.' + __name__)

@print_call(logger.info)
def get_aerosol_data(aerosol_dir, enable_aeronet, default_aerosol_value, l1t_input_dataset, bin_dir):
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug('SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    dt = datetime(
        l1t_input_dataset.scene_centre_date.year,
        l1t_input_dataset.scene_centre_date.month,
        l1t_input_dataset.scene_centre_date.day,
        l1t_input_dataset.scene_centre_time.hour,
        l1t_input_dataset.scene_centre_time.minute,
        l1t_input_dataset.scene_centre_time.second)

    return get_aerosol_value_for_region(
        dt,
        l1t_input_dataset.ll_lat,
        l1t_input_dataset.ll_lon,
        l1t_input_dataset.ur_lat,
        l1t_input_dataset.ur_lon,
        default_aerosol_value,
        aerosol_dir,
        bin_dir,
        enable_aeronet,
        l1t_input_dataset.scene_centre_lat,
        l1t_input_dataset.scene_centre_long)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)
    
    CONFIG = ProcessorConfig()
    DATA = DataManager()
    
    aod_result = get_aerosol_data(CONFIG.DIR_Aerosol, CONFIG.ENABLE_AERONET, CONFIG.DEFAULT_AEROSOL_VALUE, DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset), CONFIG.BIN_DIR)
    log_multiline(logger.info, aod_result, 'Aerosol result', '\t')
    DATA.set_item('aerosol.dat', aod_result)
