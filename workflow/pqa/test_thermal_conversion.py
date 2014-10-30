#!/bin/env python

import unittest
import numpy
import logging
from constants import PQAConstants
from thermal_conversion import get_landsat_temperature

from EOtools.DatasetDrivers import SceneDataset

class MockSatellite(object):
    def __init__(self):
        self.k = (607.76, 1260.56)
        self.TAG = 'LS5'
        
class MockDataset(object):
    def __init__(self):
        self.satellite = MockSatellite() 
        self.gain = { \
            1: 0.7658267716535433, 2: 1.4481889763779527, 3: 1.043976377952756, \
            4: 0.876023622047244, 5: 0.12035433070866142, 6: 0.0553740157480315, \
            7: 0.0655511811023622} 
        self.bias = { \
            1: -2.2858267716535465, 2: -4.288188976377967, 3: -2.213976377952804, \
            4: -2.386023622047219, 5: -0.4903543307086622, 6: 1.1826259842519669, \
            7: -0.21555118110236293}

class TestThermalConversion(unittest.TestCase):

    def test_dataset_creation(self):
        path = '/g/data/v10/L1/2009-01/LS5_TM_OTH_P51_GALPGS01-002_092_086_20090115'
        l1t_sd = SceneDataset(path)
        print  "l1t_sd.satellite.k=", l1t_sd.satellite.k
        print  "l1t_sd.gain=", l1t_sd.gain
        print  "l1t_sd.bias=", l1t_sd.bias
    
    def test_conversion(self):
        sensor = "TM"
        pq_const = PQAConstants(sensor)

        b1 = [
            [10,11,12],
            [13,14,15],
            [16,17,18]
            ]
        l1t_data = numpy.array([b1, b1, b1, b1, b1, b1, b1])

        kelvin = get_landsat_temperature(l1t_data, MockDataset(), pq_const)
        print l1t_data.shape
        print kelvin.shape
        print kelvin
        self.assertAlmostEqual(215.08184814, kelvin[0][0])
        


if __name__ == '__main__':
    
    logging.basicConfig(level=logging.INFO,
       format='%(asctime)s %(levelname)s %(message)s')
    unittest.main()
