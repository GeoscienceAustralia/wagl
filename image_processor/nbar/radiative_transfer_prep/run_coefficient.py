'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''
import logging
from ULA3.image_processor import ProcessorConfig
from ULA3.modtran import modtran_config, run_coefficient

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

    _subprocess_list = subprocess_list or modtran_config.FULL_COORD_LIST

#    sublogger.debug( '_subprocess_list = %s', _subprocess_list)

    for subprocess in _subprocess_list:
        if subprocess:
            if type(subprocess) == str:
                subprocess = subprocess.split('.')
            coordinator_list = [coordinator.upper() for coordinator in subprocess]
            assert len(coordinator_list) == 1, 'Only coordinator can be specified'
            if coordinator in modtran_config.FULL_COORD_LIST: # Coordinator only specified
                run_coefficient(coordinator, sublogger, CONFIG.work_path, CONFIG.BIN_DIR)
            else:
                raise Exception('Invalid coordinator specified: %s' % repr(coordinator))

