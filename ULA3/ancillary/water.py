'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''

import logging
from ULA3.dataset import WaterVapourDataset
from ULA3.meta import print_call

logger = logging.getLogger('root.' + __name__)

@print_call(logger.info)
def get_water_vapour_data(l1t_input_dataset, DIR_WaterVapour, DEBUG):

    water_vapour_dataset = WaterVapourDataset(DIR_WaterVapour, l1t_input_dataset.scene_centre_datetime)
    water_vapour_value = water_vapour_dataset.get_data_value(l1t_input_dataset.lonlats['CENTRE'])

    water_vapour_data = {
        'data_source': 'Water Vapour',
        'data_file': water_vapour_dataset.pathname,
        'value': water_vapour_value
        }

    return water_vapour_data
