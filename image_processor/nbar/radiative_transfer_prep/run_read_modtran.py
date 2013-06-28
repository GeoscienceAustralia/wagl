'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''
import logging
from ULA3 import DataManager
from ULA3.dataset import SceneDataset
from ULA3.modtran import read_modtran
from ULA3.image_processor import ProcessorConfig

logger = logging.getLogger('root.' + __name__)

#TODO: Implement resume
def process(subprocess_list=[], resume=False):
    """Run external modtran application.
    subprocess_list can be empty to run all combinations of coordinator and albedo, or
    its elements can consist of any combination of coordinator and/or albedo specifiers
    where modtran_config.FULL_COORD_LIST = ['UL','UM','UR','ML','MM','MR','BL','BM','BR']
    and modtran_config.FULL_ALBEDO_LIST=['0','1','T']
    """
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    DATA = DataManager()
    CONFIG = ProcessorConfig()

    assert not subprocess_list, 'Subprocesses must not be specified'

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    # List containing valid reflective band numbers as strings
    band_strings = [str(band_number) for band_number in l1t_input_dataset.bands('REFLECTIVE')]

    read_modtran(band_strings, logger, CONFIG.work_path, CONFIG.BIN_DIR)
