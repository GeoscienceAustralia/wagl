'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''

import logging
from ULA3.dataset import OzoneDataset
from ULA3.meta import print_call

logger = logging.getLogger('root.' + __name__)

@print_call(logger.info)
def get_ozone_data(l1t_input_dataset, DIR_Ozone_LUT):

    ozone_dataset = OzoneDataset(DIR_Ozone_LUT, l1t_input_dataset.scene_centre_datetime)
    ozone_value = ozone_dataset.get_data_value(l1t_input_dataset.lonlats['CENTRE'])

    ozone_data = {
        'data_source': 'Ozone',
        'data_file': ozone_dataset.pathname,
        'value': ozone_value
        }

    return ozone_data
