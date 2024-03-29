"""
Manages access methods to test data files
"""

from os.path import join as pjoin, abspath, dirname

DATA_DIR = pjoin(dirname(abspath(__file__)), "data")
TLE_DIR = pjoin(dirname(abspath(__file__)), "data")

LS5_SCENE1 = pjoin(DATA_DIR, "LANDSAT5", "LS5_TM_OTH_P51_GALPGS01-002_090_081_20090407")
LS7_SCENE1 = pjoin(DATA_DIR, "LANDSAT7", "LS7_ETM_OTH_P51_GALPGS01-002_090_081_20090415")
LS7_SCENERTC2 = pjoin(DATA_DIR, "LANDSAT7", "LE71140812021051EDC00__C2_RT")
LS8_SCENE1 = pjoin(
    DATA_DIR, "LANDSAT8", "LS8_OLITIRS_OTH_P51_GALPGS01-032_090_084_20131011"
)
LS8_SCENE1C2 = pjoin(
    DATA_DIR, "LANDSAT8", "LC08_L1TP_092084_20201029_20201106_02_T1_MTL.txt"
)
LS8_SCENE1C2 = pjoin(DATA_DIR, "LANDSAT8", "unpacked-C2")
LS8_SCENERTC2 = pjoin(DATA_DIR, "LANDSAT8", "LC81060632021051LGN00__C2_RT")

S2A_SCENE1 = pjoin(
    DATA_DIR,
    "SENTINEL2",
    "S2A_MSIL1C_20171207T002051_N0206_R116_T55JEJ_20171207T032513.zip",
)
S2B_SCENE1 = pjoin(
    DATA_DIR,
    "SENTINEL2",
    "S2B_MSIL1C_20170719T000219_N0205_R030_T56JKT_20170719T000218.zip",
)

LAND_SEA_RASTERS = pjoin(DATA_DIR, "ancillary", "Land_Sea_Rasters")

LS7_GAP_MASK = pjoin(
    DATA_DIR, "LANDSAT7", "LE07_L1TP_092084_20110809_20161206_01_T1_RR.tar"
)
LS7_NO_GAP_MASK = pjoin(
    DATA_DIR, "LANDSAT7", "LE07_L1TP_092084_19990925_20170217_01_T1_RR.tar"
)
