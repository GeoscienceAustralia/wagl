'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''
import logging
from ULA3.image_processor import ProcessorConfig
from ULA3.modtran import run_modtran, modtran_config

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
    _FULL_SUBPROCESS_LIST = [coord + '.' + albedo for coord in modtran_config.FULL_COORD_LIST for albedo in modtran_config.FULL_ALBEDO_LIST]

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

    _subprocess_list = subprocess_list or _FULL_SUBPROCESS_LIST

#    sublogger.debug( '_subprocess_list = %s', _subprocess_list)

    for subprocess in _subprocess_list:
        if subprocess:
            if type(subprocess) == str:
                subprocess = subprocess.split('.')
            module_list = [module.upper() for module in subprocess]
            assert len(module_list) <= 2, 'Only coordinator and/or albedo can be specified'
            if len(module_list) == 1: # Only one sub-module specified
                if module in modtran_config.FULL_COORD_LIST: # Coordinator only specified
                    _subsubprocess_list = [[module, albedo] for albedo in modtran_config.FULL_ALBEDO_LIST]
                elif module in modtran_config.FULL_ALBEDO_LIST: # Albedo only specified
                    _subsubprocess_list = [[coord, module] for coord in modtran_config.FULL_COORD_LIST]
                else:
                    raise Exception('run_modtran error: Invalid coordinator or albedo specified: %s' % repr(module))
                for subsubprocess in _subsubprocess_list:
                    run_modtran(subsubprocess[0], subsubprocess[1], sublogger, CONFIG.work_path, CONFIG.MODTRAN_EXE)
            else: # Single run specified by both coord and albedo
                # Swap parameters if they are in the wrong order - must have coord followed by albedo
                if (module_list[1] in modtran_config.FULL_COORD_LIST and module_list[0] in modtran_config.FULL_ALBEDO_LIST):
                    module_list = (module_list[1], module_list[0])
                assert (module_list[0] in modtran_config.FULL_COORD_LIST and
                    module_list[1] in modtran_config.FULL_ALBEDO_LIST), 'Invalid coordinator or albedo specified'

                run_modtran(module_list[0], module_list[1], sublogger, CONFIG.work_path, CONFIG.MODTRAN_EXE)

