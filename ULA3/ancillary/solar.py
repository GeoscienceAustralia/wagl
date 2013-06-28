'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''

import logging
from ULA3.meta import print_call

logger = logging.getLogger('root.' + __name__)

@print_call(logger.info)
def get_solar_irrad(l1t_input_dataset, solar_irrad_file):
    """
    Extract solar irradiance values from the specified file. One for each band

    """

    refective_bands = l1t_input_dataset.bands('REFLECTIVE')
    fpIn = open(solar_irrad_file)

    line = fpIn.readline()
    if not line:
        logger.error("ERROR: Satellite Solar Irrad file %s empty" % solar_irrad_file )
        return None

    if line.find("band solar irradiance") != 0:
        logger.error('ERROR: not a "band solar irradiance" file')
        return None

    solar_irrad = {}

    band_index = 0
    line = fpIn.readline()
    while line and band_index < len(refective_bands):
        temp = line.split()
        band_number = refective_bands[band_index] # TODO: Need to enhance file for LDCM
        solar_irrad_value = float(temp[1])
        solar_irrad[band_number] = float(temp[1])
        # Skip 0.0 value for thermal band 6 in file
        if solar_irrad_value:
            logger.debug('solar_irrad[%s] = %s', band_number, solar_irrad[band_number])
            band_index += 1
        line = fpIn.readline()

    fpIn.close()

    return solar_irrad





@print_call(logger.info)
def get_solar_distance(solar_dist_file, DOY):
    """
    Extract Earth-Sun distance for this day of the year (varies during orbit).

    """

    fpIn = open(solar_dist_file)

    line = fpIn.readline()
    if not line:
        logger.error("ERROR: Solar Distance file %s empty" % solar_dist_file)
        return 0

    solar_dist = 0
    while line:
        temp = line.split()
        doy = int(temp[0])

        if doy == DOY:
            solar_dist = float(temp[1])
        line = fpIn.readline()

    fpIn.close()

    assert solar_dist, 'Solar distance not found for DOY %s in file %s' % (DOY, solar_dist_file)

    logger.debug("Solar Distance: for DOY:%d %f" % (DOY, solar_dist))

    return solar_dist

