"""
Manages access methods to test data files
"""

from os.path import join as pjoin, abspath, dirname

DATA_DIR = pjoin(dirname(abspath(__file__)), 'data')

LS5_SCENE1 = pjoin(DATA_DIR, 'LANDSAT5', 'LS5_TM_OTH_P51_GALPGS01-002_090_081_20090407')
LS7_SCENE1 = pjoin(DATA_DIR, 'LANDSAT7', 'LS7_ETM_OTH_P51_GALPGS01-002_090_081_20090415')
LS8_SCENE1 = pjoin(DATA_DIR, 'LANDSAT8', 'LS8_OLITIRS_OTH_P51_GALPGS01-032_090_084_20131011')

LAND_SEA_RASTERS = pjoin(DATA_DIR, 'ancillary', 'Land_Sea_Rasters')
