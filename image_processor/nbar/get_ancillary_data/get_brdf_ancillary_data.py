'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''

import logging, os,datetime, threading
from osgeo import gdal, gdalconst
from ULA3 import DataManager
from ULA3.utils import log_multiline
from ULA3.dataset import SceneDataset
from ULA3.ancillary.brdf import average_brdf_value, get_brdf_dirs_modis, get_brdf_dirs_pre_modis
from ULA3.image_processor import ProcessorConfig
from ULA3.meta import print_call

#TODO: Implement resume

logger = logging.getLogger('root.' + __name__)

# keywords for BRDF parameters f0, f1, f2
_FACTOR_LIST = ['geo', 'iso', 'vol']

# Lock object to prevent multiple instances of brdf_dict
_brdf_lock = threading.Lock()



@print_call(logger.info)
def get_average_brdf_value(
    l1t_input_dataset,
    brdf_root,
    brdf_dir,
    band_string,
    band_strings,
    factor,
    work_path,
    sublogger):
        if band_string not in band_strings:
            return None

        band_number = int(band_string)
        assert band_number > 0 and band_number <= l1t_input_dataset.RasterCount, 'Invalid band number %d' % band_number

        wavelength_range = l1t_input_dataset.sensor_band_info(band_number)['WAVELENGTH']
        assert wavelength_range, 'No info for band index %d found in configuration' % band_number

        return average_brdf_value(
            #l1t_input_dataset.ll_lat, l1t_input_dataset.ll_lon,
            #l1t_input_dataset.ur_lat, l1t_input_dataset.ur_lon,
            l1t_input_dataset.lonlats['LL'][1], 
            l1t_input_dataset.lonlats['LL'][0],
            l1t_input_dataset.lonlats['UR'][1],
            l1t_input_dataset.lonlats['UR'][0],
            band_string,
            wavelength_range,
            factor,
            work_path,
            brdf_root, brdf_dir)





def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

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

    with _brdf_lock:
        brdf_dict = DATA.get_item('brdf.dat', dict)
        if brdf_dict is None:
            # Create and register new result dict
            brdf_dict = {}
            DATA.set_item('brdf.dat', brdf_dict)

    sublogger.debug( 'initial dict object for brdf.dat retrieved : %s', brdf_dict)

    # List containing valid reflective band numbers as strings
    band_strings = [str(band_number) for band_number in l1t_input_dataset.bands('REFLECTIVE')]

    # Full list of subprocesses to be run if none are specified
    full_subprocess_list = [[band_string, factor] for band_string in band_strings for factor in _FACTOR_LIST]
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
                if module in band_strings: # band_index only specified
                    _expanded_subprocess_list += [[module, factor] for factor in _FACTOR_LIST]
                elif module in _FACTOR_LIST: # Factor only specified
                    _expanded_subprocess_list += [[band_number, module] for band_number in band_strings]
                else:
                    sublogger.warning('Ignoring invalid band number or factor: %s', module_list)
                    continue
            else: # Single process specified by both band_index and factor
                # Swap parameters if they are in the wrong order - must have band no followed by factor
                if (module_list[0] in _FACTOR_LIST and module_list[1] in band_strings):
                    module_list = (module_list[1], module_list[0])
                else:
                    if not (module_list[0] in band_strings and module_list[1] in _FACTOR_LIST):
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

    # Compare the scene date and MODIS BRDF start date to select the BRDF data
    # root directory. Scene dates outside the range of the CSIRO mosaic data
    # (currently 2000-02-18 through 2013-01-09) should use the pre-MODIS,
    # Jupp-Li BRDF. Configuration setting CONFIG.FORCE_PREMODIS_BRDF=True
    # forces the use of Jupp-Li BRDF.

    modis_brdf_dir_list = sorted(os.listdir(CONFIG.DIR_BRDF))
    modis_brdf_dir_range = [modis_brdf_dir_list[0], modis_brdf_dir_list[-1]]
    modis_brdf_range = [datetime.date(*[int(x) for x in y.split('.')]) for y in modis_brdf_dir_range]

    use_jupp_li_brdf = (l1t_input_dataset.scene_centre_date < modis_brdf_range[0] or
                        l1t_input_dataset.scene_centre_date > modis_brdf_range[1])

    sublogger.info('use_jupp_li_brdf = %s' % use_jupp_li_brdf)
    sublogger.info('FORCE_PREMODIS_BRDF = %s' % CONFIG.FORCE_PREMODIS_BRDF)

    # Get the list of BRDF data directories.

    brdf_root = CONFIG.DIR_BRDF
    __func = get_brdf_dirs_modis

    if CONFIG.FORCE_PREMODIS_BRDF or use_jupp_li_brdf:
        brdf_root = CONFIG.DIR_BRDF_PREMODIS
        __func = get_brdf_dirs_pre_modis

    brdf_dirs = __func(brdf_root, l1t_input_dataset.scene_centre_date)

    sublogger.info('brdf_root = %s', brdf_root)
    sublogger.info('brdf_dirs = %s', brdf_dirs)

    assert os.path.exists(brdf_root), 'brdf_root [%s] not found. Aborting NBAR...' % brdf_root

    # Calculate average BRDF values over the scene footprint.

    for brdf_dir in brdf_dirs: # Try primary then secondary BRDF directories
        success = True
        for subprocess in _expanded_subprocess_list:
            band_string = subprocess[0]
            if band_string not in band_strings:
                sublogger.warning('Ignored invalid band number %s', band_string)
                continue

            factor = subprocess[1]
            hdf_file, average_brdf = get_average_brdf_value(
                l1t_input_dataset,
                brdf_root, brdf_dir,
                band_string, band_strings,
                factor,
                CONFIG.work_path, sublogger)

            success = bool(average_brdf)
            sublogger.info('  Average BRDF from %s = %s', hdf_file, average_brdf)
            if success: # Lookup succeeded
                brdf_dict[(int(band_string), factor)] = {
                                                         'data_source': 'BRDF',
                                                         'data_file': hdf_file,
                                                         'value': average_brdf
                                                         }
            else: # Lookup failed - try next directory (if any)
                sublogger.info('BRDF extraction failed from %s', brdf_dir)
                break # Stop lookups on this directory
        if success: # If all primary lookups succeeded
            break # Don't try secondary

    assert success, 'Primary and secondary BRDF extraction both failed. Aborting NBAR calculation'

    log_multiline(sublogger.info, brdf_dict, 'BRDF results', '\t')

    #===========================================================================
    # if resume and len(subprocess_list) == 1:
    #    # Remove bands at before the nominated band
    #    _band_list = list(_BAND_INDEX_LIST)
    #    while _band_list and _band_list[0] != int(subprocess_list[0]):
    #        _band_list.pop(0)
    # else:
    #    _band_list = [int(subprocess) for subprocess in subprocess_list] or _BAND_INDEX_LIST
    #
    #
    # for band_no in _band_list:
    #    process_band(band_no)
    #===========================================================================
