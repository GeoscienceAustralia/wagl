import unittest
import datetime
import logging
import os
import gaip.ancillary as ancillary

AOT_LOADER_PATH = os.path.abspath('../../bin')
AEROSOL_PATH =  '/g/data1/v10/eoancillarydata/aerosol/AATSR/2.0'
DEM_PATH = '/g/data1/v10/eoancillarydata/elevation/world_1deg'
OZONE_PATH = '/g/data1/v10/eoancillarydata/lookup_tables/ozone'
SOLAR_IRRAD_PATH = '/g/data1/v10/eoancillarydata/lookup_tables/solar_irradiance'
SOLAR_DIST_PATH = '/g/data1/v10/eoancillarydata/lookup_tables/earthsun_distance'

class AerosolAncillaryTest(unittest.TestCase):

    def test_run_fail(self):
        dt = datetime.datetime.now()
        self.assertRaises(OSError,
                          ancillary.run_aot_loader('foo', dt, 1., 1.,
                                                   10., 10., 'bar'))

class ElevationAncillaryTest(unittest.TestCase):
    pass

class OzoneAncillaryTest(unittest.TestCase):
    pass

class SolarIrradianceAncillaryTest(unittest.TestCase):
    pass

class SolarDistanceAncillaryTest(unittest.TestCase):
    pass


if __name__ == '__main__':
    log = logging.getLogger()
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass
    log.addHandler(NullHandler())
    unittest.main()
