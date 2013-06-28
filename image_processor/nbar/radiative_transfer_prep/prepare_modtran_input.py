'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''
import os, logging
from ULA3 import DataManager
from ULA3.image_processor import ProcessorConfig
from ULA3.dataset import SceneDataset
from ULA3.modtran import prepare_modtran_input

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    # Retrieve all required ancillary data values
    ozone_data = DATA.get_item('ozone.dat', dict)
    assert ozone_data, 'Unable to retrieve dict for ozone.dat'
    water_vapour_data = DATA.get_item('water_vapour.dat', dict)
    assert water_vapour_data, 'Unable to retrieve dict for water_vapour.dat'
    aerosol_data = DATA.get_item('aerosol.dat', dict)
    assert aerosol_data, 'Unable to retrieve dict for aerosol.dat'
    elevation_data = DATA.get_item('elevation.dat', dict)
    assert elevation_data, 'Unable to retrieve dict for elevation.dat'
    brdf_data = DATA.get_item('brdf.dat', dict)
    assert brdf_data, 'Unable to retrieve dict for brdf.dat'
    solar_dist_data = DATA.get_item('solar_dist.dat', dict)
    assert solar_dist_data, 'Unable to retrieve dict for solar_dist.dat'
    solar_irrad_data = DATA.get_item('solar_irrad.dat', dict)
    assert solar_irrad_data, 'Unable to retrieve dict for solar_irrad.dat'


    ha_file = os.path.join(CONFIG.work_path, 'HEADERANGLE')
    modtran_input_file = os.path.join(CONFIG.work_path, 'MODTRANINPUT')

    prepare_modtran_input(
        l1t_input_dataset,
        ha_file,
        modtran_input_file,
        brdf_data,
        solar_irrad_data,
        solar_dist_data,
        ozone_data,
        water_vapour_data,
        aerosol_data,
        elevation_data,
        max_view_angle=None,
        CONFIG=CONFIG
        )
