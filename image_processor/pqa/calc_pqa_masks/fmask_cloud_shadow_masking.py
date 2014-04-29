#! /usr/bin/env python


import os, logging, numpy, osr
from osgeo import gdal
from ULA3.common.pqa_result import PQAResult
from ULA3.dataset import SceneDataset
from ULA3.image_processor import ProcessorConfig
from ULA3.image_processor import constants
from ULA3 import DataManager, DataGrid
from ULA3.utils import dump_array

from cloud_shadow_masking import Cloud_Shadow

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    nbar_stack = DATA.get_item('nbar_stack', numpy.ndarray)
    assert nbar_stack is not None, 'Unable to retrieve ndarray object for nbar_stack'
    logger.debug( 'ndarray object for nbar_stack retrieved')

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    result = DATA.get_item('result.tif', PQAResult)
    assert result, 'Unable to retrieve PQAResult object for result'
    logger.debug( 'PQAResult object for result retrieved')

    kelvin_grid = DATA.get_item('kelvin.tif', DataGrid)
    assert kelvin_grid is not None, 'Unable to retrieve DataGrid object for kelvin_grid'
    logger.debug( 'DataGrid object for kelvin_grid retrieved')
    kelvin_array = kelvin_grid.array

    pqa_temp_output = DATA.get_item('pqa_temp_output.dat', str)
    assert pqa_temp_output, 'Unable to retrieve string object for pqa_temp_output'
    logger.debug( 'string object for pqa_temp_output retrieved')

    #===========================================================================
    # fmask_cloud_mask = DATA.get_item('fmask_cloud_mask', numpy.ndarray)
    # assert fmask_cloud_mask is not None, 'Unable to retrieve ndarray object for fmask_cloud_mask'
    # logger.debug( 'ndarray object for fmask_cloud_mask retrieved')
    #
    # contiguity_mask = DATA.get_item('contiguity_mask', numpy.ndarray)
    # assert contiguity_mask is not None, 'Unable to retrieve ndarray object for contiguity_mask'
    # logger.debug( 'ndarray object for contiguity_mask retrieved')
    #
    # land_sea_mask = DATA.get_item('land_sea_mask', numpy.ndarray)
    # assert land_sea_mask is not None, 'Unable to retrieve ndarray object for land_sea_mask'
    # logger.debug( 'ndarray object for land_sea_mask retrieved')
    #===========================================================================

    # Initialise the PQA constants
    pq_const = constants.pqaContants(l1t_input_dataset.sensor)

    if pq_const.run_cloud_shadow: # TM/ETM/OLI_TIRS
        #contiguity_mask = result.get_mask(CONFIG.pqa_test_index['CONTIGUITY'])
        #fmask_cloud_mask = result.get_mask(CONFIG.pqa_test_index['FMASK'])
        #land_sea_mask = result.get_mask(CONFIG.pqa_test_index['LAND_SEA'])
        contiguity_mask = result.get_mask(pq_const.contiguity)
        fmask_cloud_mask = result.get_mask(pq_const.fmask)
        land_sea_mask = result.get_mask(pq_const.land_sea)

        if pq_const.oli_tirs:
            mask = Cloud_Shadow(image_stack=nbar_stack[1:,:,:], kelvin_array=kelvin_array,
                                cloud_mask=fmask_cloud_mask, input_dataset=l1t_input_dataset,
                                land_sea_mask=land_sea_mask, contiguity_mask=contiguity_mask,
                                cloud_algorithm='FMASK', growregion=True)
        else: # TM or ETM
            mask = Cloud_Shadow(image_stack=nbar_stack, kelvin_array=kelvin_array,
                                cloud_mask=fmask_cloud_mask, input_dataset=l1t_input_dataset,
                                land_sea_mask=land_sea_mask, contiguity_mask=contiguity_mask,
                                cloud_algorithm='FMASK', growregion=True)

        #bit_index = CONFIG.pqa_test_index['FMASK_SHADOW']
        bit_index = pq_const.fmask_shadow
        result.set_mask(mask, bit_index)
        if CONFIG.debug:
            dump_array(mask,
                       os.path.join(CONFIG.work_path, 'mask_%02d.tif' % bit_index),
                       l1t_input_dataset)
    else: # OLI/TIRS only
        logger.debug('Cloud Shadow Algorithm Not Run! %s sensor not configured for the cloud shadow algorithm.'%l1t_input_dataset.sensor)
        cloud_shadow_logfile = open(os.path.join(pqa_temp_output, 'scene01', 'FMASK_CLOUD_SHADOW_LOGFILE.txt'), 'w')
        cloud_shadow_logfile.write("Test Not Run!\n")
        cloud_shadow_logfile.write("%s sensor not configured for the Cloud Shadow algorithm."%l1t_input_dataset.sensor)
        cloud_shadow_logfile.close()

