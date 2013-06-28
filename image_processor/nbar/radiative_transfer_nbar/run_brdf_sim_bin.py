'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''
import logging
from ULA3 import DataManager
from ULA3.dataset import SceneDataset
from ULA3.modtran import run_brdf_sim_bin
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

    for subprocess in _subprocess_list:
        if subprocess:
            if type(subprocess) == str:
                subprocess = subprocess.split('.')
            assert len(subprocess) == 1, 'Only band can be specified'
            if subprocess[0] in band_strings:
                run_brdf_sim_bin(l1t_input_dataset, int(subprocess[0]), sublogger, CONFIG.work_path, CONFIG.BIN_DIR, CONFIG.DEBUG, CONFIG.L7_SLC_DATE)
            else:
                logger.warning('Ignoring invalid band: %s' % subprocess[0])

