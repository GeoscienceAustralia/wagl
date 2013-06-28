'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''
import logging, threading
from ULA3 import DataManager
from ULA3.image_processor import ProcessorConfig
from ULA3.dataset import SceneDataset
from ULA3.modtran import run_bilinear_ortho, modtran_config
#from sub_process import execute
#from log_multiline import log_multiline
#import modtran_config

logger = logging.getLogger('root.' + __name__)

_output_lock = threading.Lock()

#TODO: Implement resume
def process(subprocess_list=[], resume=False):
    """Run external modtran application.
    subprocess_list can be empty to run all combinations of coordinator and albedo, or
    its elements can consist of any combination of coordinator and/or albedo specifiers
    where modtran_config.FULL_COORD_LIST = ['UL','UM','UR','ML','MM','MR','BL','BM','BR']
    and modtran_config.FULLmodtran_config.FULL_FACTOR_LIST=['0','1','T']
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

    satellite = l1t_input_dataset.satellite
    assert satellite, 'Unable to find Satellite object'

    # List containing valid band numbers as strings
    band_string_list = [str(band_number) for band_number in l1t_input_dataset.bands('REFLECTIVE')]

    # Full list of subprocesses to be run if none are specified
    full_subprocess_list = [[band_number, factor] for band_number in band_string_list for factor in modtran_config.FULL_FACTOR_LIST]
    _subprocess_list = subprocess_list or full_subprocess_list
    sublogger.debug('_subprocess_list = %s', _subprocess_list)

    _expanded_subprocess_list = []
    for subprocess in _subprocess_list:
        if subprocess:
            if type(subprocess) == str:
                subprocess = subprocess.split('.')
            # Inspect individual sub-process and expand to both band_index and factor if required
            module_list = [module.lower() for module in subprocess]
            assert len(module_list) <= 2, 'Only band_index and/or factor can be specified'
            if len(module_list) == 1: # Only one sub-module specified
                module = module_list[0]
                if module in band_string_list: # band_index only specified
                    _expanded_subprocess_list += [[module, factor] for factor in modtran_config.FULL_FACTOR_LIST]
                elif module in modtran_config.FULL_FACTOR_LIST: # Factor only specified
                    _expanded_subprocess_list += [[band_number, module] for band_number in band_string_list]
                else:
                    sublogger.warning('Ignoring invalid band number or factor: %s', module_list)
                    continue
            else: # Single process specified by both band_index and factor
                # Swap parameters if they are in the wrong order - must have band no followed by factor
                if (module_list[0] in modtran_config.FULL_FACTOR_LIST and module_list[1] in band_string_list):
                    module_list = (module_list[1], module_list[0])
                else:
                    if not (module_list[0] in band_string_list and module_list[1] in modtran_config.FULL_FACTOR_LIST):
                        sublogger.warning('Ignoring invalid band number or factor: %s', module_list)
                        continue
                _expanded_subprocess_list.append(module_list)
            sublogger.debug('module_list = %s', module_list)

    # Remove duplicates (if any)
    _expanded_subprocess_list = [_expanded_subprocess_list[i] for i,x in enumerate(_expanded_subprocess_list) if x not in _expanded_subprocess_list[i+1:]]
    sublogger.debug('_expanded_subprocess_list = %s', _expanded_subprocess_list)

    if not _expanded_subprocess_list:
        sublogger.warning('No bands or factors to process')
        return

    with _output_lock:
        outputs = DATA.get_item("bilinear_ortho_outputs", item_type=dict)
        if outputs is None:
            outputs = DATA.save_item(item_path="bilinear_ortho_outputs", item_instance={})

    for subprocess in _expanded_subprocess_list:
        band_number = int(subprocess[0])
        factor = subprocess[1]
        outputs[(band_number, factor)] = run_bilinear_ortho(band_number, factor, sublogger, CONFIG.work_path, CONFIG.BIN_DIR)

# What this script should be more like...
#    [threading.Thread(target=run_bilinear_ortho, args=(b, f, sublogger, CONFIG.work_path, CONFIG.BIN_DIR)).run() for b in l1t_input_dataset.bands('REFLECTIVE') for f in factors]
#    [run_bilinear_ortho(b, f, sublogger, CONFIG.work_path, CONFIG.BIN_DIR) for b in l1t_input_dataset.bands('REFLECTIVE') for f in factors]

    DATA.save_item(item_path="bilinear_ortho_outputs", item_type=dict)
