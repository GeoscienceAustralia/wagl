#!/usr/bin/env python

"""MTL Metadata module

Author: Alex Ip (alex.ip@ga.gov.au)
"""

import logging
from . import Metadata

logger = logging.getLogger('root.' + __name__)

class MTLMetadata(Metadata):
    """Subclass of Metadata to manage MTL data
    """
    # Class variable holding metadata type string
    _metadata_type_id = 'MTL'
    _filename_pattern = '.*\.mtl' # Default RegEx for finding metadata file.

    def _parse_mtl_string(self, mtl_string, tree_dict=None):
        """Private function to parse a block of MTL text into a nested dict
        Exposed to allow unit testing using a string
        Arguments:
            mtl_string: MTL text block (string)
        Returns:
            MTL parameter nested dictionary.
        """
        def parse_mtl_lines(line_list, tree_dict, level=0):
            """Recursive helper function to parse MTL lines
            """
            while line_list:
                line = line_list.pop(0)
                try:
                    key, value = [s.strip(' "') for s in line.split('=')]
                    key = key.replace(' ', '_')
                    if key == 'GROUP':
                        group = value
                        tree_dict[group] = {}
                        parse_mtl_lines(line_list, tree_dict[group], level + 1)
                    elif key == 'END_GROUP':
                        break
                    else:
                        tree_dict[key] = value

                    logger.debug('%s%s = %s', '  ' * level, key, value)
                except ValueError:
                    pass

        tree_dict = tree_dict or self._metadata_dict

        # Build nested dictionaries by group.
        parse_mtl_lines([line for line in mtl_string.splitlines() if line], tree_dict)


    def read_file(self, filename = None):
        """Function to parse an MTL metadata file and store the results in self._metadata_dict
        Argument:
            filename: MTL Metadata file to be parsed and stored
        Returns:
            nested dict containing metadata
        """
        filename = filename or self._filename
        assert filename, 'Filename must be specified'

        logger.debug('Parsing MTL file %s', filename)

        # Open MTL document
        try:
            infile = open(filename, 'r')
            assert infile is not None, 'Unable to open MTL file ' + filename
            return self._parse_mtl_string(infile.read())
        finally:
            infile.close()

        return self._metadata_dict

    def write_file(self, filename = None):
        """Function write the metadata contained in self._metadata_dict to an MTL file
        Argument:
            filename: Metadata file to be written
        """
        filename = filename or self._filename
        assert filename, 'Filename must be specified'
        # ToDo - needs to be implemented
        pass


def main():
    # Test data from file LS7_ETM_OTH_P51_GALPGS01_092_085_20100315/scene01/L71092085_08520100315_MTL.txt
    TESTMTL = """GROUP = L1_METADATA_FILE
  GROUP = METADATA_FILE_INFO
    ORIGIN = "Image courtesy of the U.S. Geological Survey"
    REQUEST_ID = "0072010074092_00085"
    PRODUCT_CREATION_TIME = 2011-10-25T05:14:51Z
    STATION_ID = "EDC"
    LANDSAT7_XBAND = "2"
    GROUND_STATION = "ASA"
    LPS_PROCESSOR_NUMBER = 3
    DATEHOUR_CONTACT_PERIOD = "1007423"
    SUBINTERVAL_NUMBER = "01"
  END_GROUP = METADATA_FILE_INFO
  GROUP = PRODUCT_METADATA
    PRODUCT_TYPE = "L1T"
    ELEVATION_SOURCE = "SRTM3"
    PROCESSING_SOFTWARE = "LPGS_11.4.0"
    EPHEMERIS_TYPE = "DEFINITIVE"
    SPACECRAFT_ID = "Landsat7"
    SENSOR_ID = "ETM+"
    SENSOR_MODE = "BUMPER"
    ACQUISITION_DATE = 2010-03-15
    SCENE_CENTER_SCAN_TIME = 23:55:11.6634891Z
    WRS_PATH = 92
    STARTING_ROW = 85
    ENDING_ROW = 85
    BAND_COMBINATION = "123456678"
    PRODUCT_UL_CORNER_LAT = -35.0700000
    PRODUCT_UL_CORNER_LON = 144.8700000
    PRODUCT_UR_CORNER_LAT = -35.0700000
    PRODUCT_UR_CORNER_LON = 147.5900000
    PRODUCT_LL_CORNER_LAT = -37.0000000
    PRODUCT_LL_CORNER_LON = 144.8700000
    PRODUCT_LR_CORNER_LAT = -37.0000000
    PRODUCT_LR_CORNER_LON = 147.5900000
    PRODUCT_UL_CORNER_MAPX = 14487000.000
    PRODUCT_UL_CORNER_MAPY = -3507000.000
    PRODUCT_UR_CORNER_MAPX = 14759000.000
    PRODUCT_UR_CORNER_MAPY = -3507000.000
    PRODUCT_LL_CORNER_MAPX = 14487000.000
    PRODUCT_LL_CORNER_MAPY = -3700000.000
    PRODUCT_LR_CORNER_MAPX = 14759000.000
    PRODUCT_LR_CORNER_MAPY = -3700000.000
    PRODUCT_SAMPLES_PAN = 21761
    PRODUCT_LINES_PAN = 15441
    PRODUCT_SAMPLES_REF = 10881
    PRODUCT_LINES_REF = 7721
    PRODUCT_SAMPLES_THM = 5441
    PRODUCT_LINES_THM = 3861
    BAND1_FILE_NAME = "L71092085_08520100315_B10.FST"
    BAND2_FILE_NAME = "L71092085_08520100315_B20.FST"
    BAND3_FILE_NAME = "L71092085_08520100315_B30.FST"
    BAND4_FILE_NAME = "L71092085_08520100315_B40.FST"
    BAND5_FILE_NAME = "L71092085_08520100315_B50.FST"
    BAND61_FILE_NAME = "L71092085_08520100315_B61.FST"
    BAND62_FILE_NAME = "L72092085_08520100315_B62.FST"
    BAND7_FILE_NAME = "L72092085_08520100315_B70.FST"
    BAND8_FILE_NAME = "L72092085_08520100315_B80.FST"
    GCP_FILE_NAME = "L71092085_08520100315_GCP.txt"
    METADATA_L1_FILE_NAME = "L71092085_08520100315_MTL.txt"
    CPF_FILE_NAME = "L7CPF20100101_20100331_08"
  END_GROUP = PRODUCT_METADATA
  GROUP = MIN_MAX_RADIANCE
    LMAX_BAND1 = 191.600
    LMIN_BAND1 = -6.200
    LMAX_BAND2 = 196.500
    LMIN_BAND2 = -6.400
    LMAX_BAND3 = 152.900
    LMIN_BAND3 = -5.000
    LMAX_BAND4 = 241.100
    LMIN_BAND4 = -5.100
    LMAX_BAND5 = 31.060
    LMIN_BAND5 = -1.000
    LMAX_BAND61 = 17.040
    LMIN_BAND61 = 0.000
    LMAX_BAND62 = 12.650
    LMIN_BAND62 = 3.200
    LMAX_BAND7 = 10.800
    LMIN_BAND7 = -0.350
    LMAX_BAND8 = 243.100
    LMIN_BAND8 = -4.700
  END_GROUP = MIN_MAX_RADIANCE
  GROUP = MIN_MAX_PIXEL_VALUE
    QCALMAX_BAND1 = 255.0
    QCALMIN_BAND1 = 1.0
    QCALMAX_BAND2 = 255.0
    QCALMIN_BAND2 = 1.0
    QCALMAX_BAND3 = 255.0
    QCALMIN_BAND3 = 1.0
    QCALMAX_BAND4 = 255.0
    QCALMIN_BAND4 = 1.0
    QCALMAX_BAND5 = 255.0
    QCALMIN_BAND5 = 1.0
    QCALMAX_BAND61 = 255.0
    QCALMIN_BAND61 = 1.0
    QCALMAX_BAND62 = 255.0
    QCALMIN_BAND62 = 1.0
    QCALMAX_BAND7 = 255.0
    QCALMIN_BAND7 = 1.0
    QCALMAX_BAND8 = 255.0
    QCALMIN_BAND8 = 1.0
  END_GROUP = MIN_MAX_PIXEL_VALUE
  GROUP = PRODUCT_PARAMETERS
    CORRECTION_METHOD_GAIN_BAND1 = "CPF"
    CORRECTION_METHOD_GAIN_BAND2 = "CPF"
    CORRECTION_METHOD_GAIN_BAND3 = "CPF"
    CORRECTION_METHOD_GAIN_BAND4 = "CPF"
    CORRECTION_METHOD_GAIN_BAND5 = "CPF"
    CORRECTION_METHOD_GAIN_BAND61 = "CPF"
    CORRECTION_METHOD_GAIN_BAND62 = "CPF"
    CORRECTION_METHOD_GAIN_BAND7 = "CPF"
    CORRECTION_METHOD_GAIN_BAND8 = "CPF"
    CORRECTION_METHOD_BIAS = "IC"
    BAND1_GAIN = "H"
    BAND2_GAIN = "H"
    BAND3_GAIN = "H"
    BAND4_GAIN = "L"
    BAND5_GAIN = "H"
    BAND6_GAIN1 = "L"
    BAND6_GAIN2 = "H"
    BAND7_GAIN = "H"
    BAND8_GAIN = "L"
    BAND1_GAIN_CHANGE = "0"
    BAND2_GAIN_CHANGE = "0"
    BAND3_GAIN_CHANGE = "0"
    BAND4_GAIN_CHANGE = "0"
    BAND5_GAIN_CHANGE = "0"
    BAND6_GAIN_CHANGE1 = "0"
    BAND6_GAIN_CHANGE2 = "0"
    BAND7_GAIN_CHANGE = "0"
    BAND8_GAIN_CHANGE = "0"
    BAND1_SL_GAIN_CHANGE = 0
    BAND2_SL_GAIN_CHANGE = 0
    BAND3_SL_GAIN_CHANGE = 0
    BAND4_SL_GAIN_CHANGE = 0
    BAND5_SL_GAIN_CHANGE = 0
    BAND6_SL_GAIN_CHANGE1 = 0
    BAND6_SL_GAIN_CHANGE2 = 0
    BAND7_SL_GAIN_CHANGE = 0
    BAND8_SL_GAIN_CHANGE = 0
    SUN_AZIMUTH = 53.7730377
    SUN_ELEVATION = 41.5610274
    OUTPUT_FORMAT = "FASTL7A"
  END_GROUP = PRODUCT_PARAMETERS
  GROUP = CORRECTIONS_APPLIED
    STRIPING_BAND1 = "NONE"
    STRIPING_BAND2 = "NONE"
    STRIPING_BAND3 = "NONE"
    STRIPING_BAND4 = "NONE"
    STRIPING_BAND5 = "NONE"
    STRIPING_BAND61 = "NONE"
    STRIPING_BAND62 = "NONE"
    STRIPING_BAND7 = "NONE"
    STRIPING_BAND8 = "NONE"
    BANDING = "N"
    COHERENT_NOISE = "Y"
    MEMORY_EFFECT = "N"
    SCAN_CORRELATED_SHIFT = "N"
    INOPERABLE_DETECTORS = "N"
    DROPPED_LINES = "N"
  END_GROUP = CORRECTIONS_APPLIED
  GROUP = PROJECTION_PARAMETERS
    REFERENCE_DATUM = "WGS84"
    REFERENCE_ELLIPSOID = "WGS84"
    GRID_CELL_SIZE_PAN = 0.0001250
    GRID_CELL_SIZE_THM = 0.0005000
    GRID_CELL_SIZE_REF = 0.0002500
    ORIENTATION = "NUP"
    RESAMPLING_OPTION = "CC"
    SCAN_GAP_INTERPOLATION = 2
    MAP_PROJECTION = "EQR"
  END_GROUP = PROJECTION_PARAMETERS
  GROUP = EQR_PARAMETERS
    LATITUDE_OF_FIRST_STANDARD_PARALLEL = 0.0000000
    LONGITUDE_OF_CENTRAL_MERIDIAN = 0.0000000
    FALSE_EASTING = 0.0000000
    FALSE_NORTHING = 0.0000000
    FALSE_EASTING_NORTHING_UNITS = "meters"
  END_GROUP = EQR_PARAMETERS
END_GROUP = L1_METADATA_FILE
END"""

    # Instantiate empty MTLMetadata object and parse test string
    mtl_object = MTLMetadata()
    mtl_object._parse_mtl_string(TESTMTL)

    assert mtl_object.metadata_dict, 'No metadata_dict created'
    assert mtl_object.tree_to_list(), 'Unable to create list from metadata_dict'
    assert mtl_object.get_metadata('L1_METADATA_FILE,PRODUCT_METADATA,SPACECRAFT_ID'.split(',')), 'Unable to find value L1_METADATA_FILE,PRODUCT_METADATA,SPACECRAFT_ID'
    assert mtl_object.get_metadata('...,SPACECRAFT_ID'.split(',')), 'Unable to find value ...,SPACECRAFT_ID'
    assert not mtl_object.get_metadata('RUBBERCHICKEN'.split(',')), 'Found nonexistent key RUBBERCHICKEN'
    mtl_object.set_metadata_node('L1_METADATA_FILE,PRODUCT_METADATA,SPACECRAFT_ID'.split(','), 'Rubber Chicken')
    assert mtl_object.get_metadata('...,SPACECRAFT_ID'.split(',')), 'Unable to change ...,SPACECRAFT_ID to "Rubber Chicken"'
    mtl_object.merge_metadata_dicts({'RUBBERCHICKEN': 'Rubber Chicken'}, mtl_object.metadata_dict)
    assert mtl_object.get_metadata('RUBBERCHICKEN'.split(',')), 'Unable to find value for key RUBBERCHICKEN'
#    print mtl_object.tree_to_list()

if __name__ == '__main__':
    main()

