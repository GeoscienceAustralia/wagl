import os, logging
import numpy, numexpr
from osgeo import gdal
from ULA3.dataset import SceneDataset
from ULA3.common.pqa_result import PQAResult
from ULA3.image_processor import ProcessorConfig
from ULA3 import DataManager
from ULA3.utils import dump_array

logger = logging.getLogger('root.' + __name__)

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

    l1t_stack = DATA.get_item('l1t_stack', numpy.ndarray)
    assert l1t_stack is not None, 'Unable to retrieve ndarray object for l1t_stack'
    sublogger.debug( 'ndarray object for l1t_stack retrieved')

    result = DATA.get_item('result.tif', PQAResult)
    assert result, 'Unable to retrieve PQAResult object for result'
    sublogger.debug( 'PQAResult object for result retrieved')

    def SaturationMask(band_array, use_numexpr=True):
        """
        Creates a saturation mask for a single band.

        Values of 1 and 255 are considered under and over saturated respectively.

        :param band_array:
            A 2D numpy array containing a single Landsat band.

        :return:
            A 2D Numpy Boolean array with True == Unsaturated.
        """

        if len(band_array) == 0: return None
        assert type(band_array) == numpy.ndarray, 'Input is not valid'

        # Work-around for non-thread-safe numexpr expression
        if use_numexpr:
            sublogger.debug('numexpr used: numexpr.evaluate("(band_array > 1) & (band_array < 255)")')
            mask = numexpr.evaluate("(band_array > 1) & (band_array < 255)")
        else:
            sublogger.debug('numpy used: (band_array > 1) & (band_array < 255)')
            mask = (band_array > 1) & (band_array < 255)

        return mask

    # List of all band numbers
    full_band_list = sorted(set(l1t_input_dataset.bands('REFLECTIVE')) |
                        set(l1t_input_dataset.bands('THERMAL')))
    assert len(full_band_list) == l1t_stack.shape[0], 'Mismatch between array depth and band count'
    sublogger.debug('full_band_list = %s', full_band_list)

    # Default to processing all reflective & thermal bands
    band_list = [int(subprocess) for subprocess in subprocess_list] or full_band_list

    # Determine bit indices of all available bands from file numbers
    bit_index_list = [CONFIG.pqa_test_index['SATURATION'][band_file_no] for band_file_no in sorted(
        set(l1t_input_dataset.satellite.BAND_TYPES['REFLECTIVE']) |
        set(l1t_input_dataset.satellite.BAND_TYPES['THERMAL'])
        )]
    assert len(bit_index_list) == l1t_stack.shape[0], 'Unable to determine bit indices for all bands'
    sublogger.debug('bit_index_list = %s', bit_index_list)

    for band in band_list:
        if band not in full_band_list:
            sublogger.warning('Ignoring invalid band number: %d', band)
            continue

        band_index = full_band_list.index(band)
        sublogger.debug('Processing saturation for band = %d, band_index = %d, bit_index = %d', band, band_index, bit_index_list[band_index])

        band_array = l1t_stack[band_index]

        # Use numpy for single bands run in parallel - numexpr is not thread safe
        mask = SaturationMask(band_array, (len(band_list) > 1 or CONFIG.sequential))

        bit_index = bit_index_list[band_index]
        result.set_mask(mask, bit_index)
        if CONFIG.debug:
            dump_array(mask,
                       os.path.join(CONFIG.work_path, 'mask_%02d.tif' % bit_index),
                       l1t_input_dataset)

        # Copy results for first thermal band to second one if there is only one available
        if bit_index == 5 and 6 not in bit_index_list: # If this is thermal band 1 and there is no thermal band 2 in this dataset
            bit_index = 6
            sublogger.debug('Copying thermal band mask to bit %d', bit_index)
            result.set_mask(mask, bit_index)


