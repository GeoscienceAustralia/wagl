import os, numpy, logging, osr
from numpy.core.numeric import nan
from osgeo import gdal
from ULA3 import DataManager
from ULA3.dataset import SceneDataset
from ULA3.common.pqa_result import PQAResult
from ULA3.image_processor import ProcessorConfig

logger = logging.getLogger('root.' + __name__)

NaN = numpy.float32(nan)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    result = DATA.get_item('result.tif', PQAResult)
    assert result, 'Unable to retrieve PQAResult object for result'
    logger.debug( 'PQAResult object for result retrieved')

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    pqa_temp_output = DATA.get_item('pqa_temp_output.dat', str)
    assert pqa_temp_output, 'Unable to retrieve pqa_temp_output string'
    logger.debug('string for %s retrieved', pqa_temp_output)

    pqa_dataset_id = DATA.get_item('pqa_dataset_id.dat', str)
    assert pqa_temp_output, 'Unable to retrieve pqa_dataset_id string'
    logger.debug('string %s for pqa_dataset_id retrieved', pqa_dataset_id)

    tif_filename=os.path.join(pqa_temp_output,
            'scene01',
            pqa_dataset_id + '_' + result.test_string + '.tif'
            )
    result.save_image(
        filename=tif_filename,
        image_type='GTiff',
        dtype=None,
        projection_ref=l1t_input_dataset.GetProjection(),
        geotransform=l1t_input_dataset.GetGeoTransform(),
        convert_angles=False)

    logger.info('PQA output saved to %s', tif_filename)
