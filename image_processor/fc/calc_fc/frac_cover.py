import os, logging
from ULA3.fc import fractional_cover
from ULA3.image_processor import ProcessorConfig
from ULA3 import DataManager
from ULA3.dataset import SceneDataset

logger = logging.getLogger('root.' + __name__)

def process(*args, **kwargs):
    CONFIG = ProcessorConfig()
    DATA = DataManager()

    #TODO This would need to change if output format changed.
    asfloat32 = False

    nbar_input_dataset = DATA.get_item(CONFIG.input['nbar']['path'], SceneDataset)
    assert nbar_input_dataset, 'Unable to retrieve input scene dataset'
    logger.debug('SceneDataset object for %s retrieved', nbar_input_dataset.pathname)

    fc_temp_output = DATA.get_item('fc_temp_output.dat', str)
    assert fc_temp_output, 'Unable to retrieve fc_temp_output string'
    logger.debug('string for fc_temp_output retrieved')

    fc_dataset_id = DATA.get_item('fc_dataset_id.dat', str)
    assert fc_dataset_id, 'Unable to retrieve fc_dataset_id string'
    logger.debug('string for fc_dataset_id retrieved')

    fc_result = fractional_cover(nbar_data_path=CONFIG.input["nbar"]["path"],
                                 asfloat32=asfloat32,
                                 fc_data_path=fc_temp_output,
                                 single_tif=False)
