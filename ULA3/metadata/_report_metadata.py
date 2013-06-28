'''
Created on Jun 22, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
'''

import re, os, logging
from . import Metadata

logger = logging.getLogger('root.' + __name__)

class ReportMetadata(Metadata):
    """Subclass of Metadata to manage Report data
    """
    # Class variable holding metadata type string
    _metadata_type_id = 'REPORT'
    _filename_pattern = 'report\.txt$' # Default RegEx for finding metadata file.
    _orient_map = {'NUP': 0.0, 'EUP': 90.0, 'SUP': 180.0, 'WUP': -90.0}

    # Define template text for output. N.B: All format specifiers are for strings
    _template_text="""
                  GA NBAR Processing Report
                  -------------------------

Dataset ID:         %(DATASET_ID)-s
Satellite:          %(SATELLITE)-20.20sSensor:             %(SENSOR)-s
Camera Number:      %(CAMERA_NUMBER)-20.20sSensor Mode:        %(SENSOR_MODE)-s

Ground Station:     %(GROUND_STATION)-s

Processing Level:   %(PROCESSING_LEVEL)-20.20sResampling:         %(RESAMPLING)-s
Map Projection:     %(MAP_PROJECTION)-20.20sZone:               %(ZONE)-s
Earth Ellipsoid:    %(EARTH_ELLIPSOID)-20.20sDatum:              %(DATUM)-s

Path/Strip No.:     %(PATHSTRIP_NO)-20.20sRow No:             %(ROW_NO).20s
Image Lines:        %(IMAGE_LINES)-20.20sImage Pixels:       %(IMAGE_PIXELS)-s
Image Orientation:  %(IMAGE_ORIENTATION)-20.20sPixel Size:         %(PIXEL_SIZE)-s
Output Bands:       %(OUTPUT_BANDS)-s

Scene Centre Lat:   %(SCENE_CENTRE_LAT)-20.20sScene Centre Long:  %(SCENE_CENTRE_LONG)-s
Scene Centre Date:  %(SCENE_CENTRE_DATE)-20.20sScene Centre Time:  %(SCENE_CENTRE_TIME)-s

Product Format:     %(PRODUCT_FORMAT)-20.20sInterleaving:       %(INTERLEAVING)-s

Completion Date:    %(COMPLETION_DATE)-20.20sCompletion Time:    %(COMPLETION_TIME)-s

Termination Status: %(TERMINATION_STATUS)-s



                  PRODUCT FORMATTING
                  ------------------

Product Scene Centre Location (lat/long) : %(PRODUCT_SCENE_CENTRE_LOCATION_LATLONG)-s
Product Scene Centre Date/Time (yyyy-mm-dd hh:mm:ss.ss): %(PRODUCT_SCENE_CENTRE_DATETIME_YYYYMMDD_HHMMSSSS)-s

Product Extent:

  Lat:   %(UL_LAT)11.11s ------------ Lat:   %(UR_LAT)11.11s
  Long:  %(UL_LONG)11.11s              Long:  %(UR_LONG)11.11s
  North: %(UL_NORTH)11.11s              North: %(UR_NORTH)11.11s
  East:  %(UL_EAST)11.11s              East:  %(UR_EAST)11.11s
  |                                                |
  |                                                |
  |                                                |
  |                                                |
  |                                                |
  |                                                |
  Lat:   %(LL_LAT)11.11s              Lat:   %(LR_LAT)11.11s
  Long:  %(LL_LONG)11.11s              Long:  %(LR_LONG)11.11s
  North: %(LL_NORTH)11.11s              North: %(LR_NORTH)11.11s
  East:  %(LL_EAST)11.11s ------------ East:  %(LR_EAST)11.11s

"""

    # Need to ensure that any variables referenced in _template_text are keys in _template_dict
    # Note that all values are pre-formatted strings
    _template_dict={
        'CAMERA_NUMBER': 'N/A',
        'COMPLETION_DATE': None,
        'COMPLETION_TIME': None,
        'DATASET_ID': None,
        'DATUM': None,
        'EARTH_ELLIPSOID': None,
        'GROUND_STATION': None,
        'IMAGE_LINES': None,
        'IMAGE_ORIENTATION': None,
        'IMAGE_PIXELS': None,
        'INTERLEAVING': None,
        'LL_EAST': None,
        'LL_LAT': None,
        'LL_LONG': None,
        'LL_NORTH': None,
        'LR_EAST': None,
        'LR_LAT': None,
        'LR_LONG': None,
        'LR_NORTH': None,
        'MAP_PROJECTION': None,
        'OUTPUT_BANDS': None,
        'PATHSTRIP_NO': None,
        'PIXEL_SIZE': None,
        'PROCESSING_LEVEL': None,
        'PRODUCT_FORMAT': None,
        'PRODUCT_SCENE_CENTRE_DATETIME_YYYYMMDD_HHMMSSSS': '0000-00-00 00:00:00.000000',
        'PRODUCT_SCENE_CENTRE_LOCATION_LATLONG': '0 / 0 deg.',
        'RESAMPLING': None,
        'ROW_NO': None,
        'SATELLITE': None,
        'SCENE_CENTRE_DATE': None,
        'SCENE_CENTRE_LAT': None,
        'SCENE_CENTRE_LONG': None,
        'SCENE_CENTRE_TIME': '00:00:00.000000',
        'SCENE_SHIFT': None,
        'SENSOR': None,
        'SENSOR_MODE': None,
        'SUN_AZIMUTH': None,
        'SUN_ELEVATION': None,
        'TERMINATION_STATUS': None,
        'UL_EAST': None,
        'UL_LAT': None,
        'UL_LONG': None,
        'UL_NORTH': None,
        'UR_EAST': None,
        'UR_LAT': None,
        'UR_LONG': None,
        'UR_NORTH': None,
        'ZONE': None
        }

    def __init__(self, source=None):
        """Instantiates ReportMetadata object. Overrides Metadata method
        """
        self._value_dict = ReportMetadata._template_dict.copy() # Ensure that template will always work
        Metadata.__init__(self, source); # Call inherited constructor
        self._metadata_dict = self._value_dict # Start with a fully populated dict

    def parse_report_string(self, report_string, tree_dict=None):

        def get_name_value(input_text):
            """Returns (name, value) tuple or None if no "name : value" found
            """
            m = re.match('^\W*(.*?):\s+(.*?)\s*-*$', input_text) # Look for "name: value" expression
            if m and len(m.groups()) == 2:
                result = m.groups()
                # Remove any non-alphanumeric characters from name and convert to upper case
                result = (re.sub('[^\w_]*', '', re.sub('\s+', '_', result[0].strip().upper())),
                          result[1])
            else:
                result = None
            logger.debug('get_name_value result = %s', repr(result))
            return result

        def add_name_value(name_value, tree_dict):
            """Add (name, value) to tree_dict
            """
            if (name_value):
                logger.debug('tree_dict[%s] = %s', name_value[0], name_value[1])
                tree_dict[name_value[0]] = name_value[1]
                self._value_dict[name_value[0]] = name_value[1]

        logger.debug('parse_report_text(%s, %s) called', report_string, tree_dict)

        tree_dict = tree_dict or self._metadata_dict

        y_prefix = 'UL'
        x_prefix = 'LR'

        section = 0
        for line in [line for line in report_string.splitlines() if line]:
            logger.debug('line = %s', line)
            if section == 0:
                if re.match('(Product Extent)|(PRODUCT FORMATTING)', line):
                    section += 1
                    extent_y = 0
                    continue
                else:
                    lh_name_value = get_name_value(line[0:39])
                    rh_name_value = get_name_value(line[40:])
                    if lh_name_value and rh_name_value: # Two "name: value" expressions found on the same line
                        add_name_value(lh_name_value, tree_dict)
                        add_name_value(rh_name_value, tree_dict)
                    else: # Zero or one "name: value" expression found on the line
                        name_value = get_name_value(line)
                        if name_value: # One "name: value" expression found on the line
                            add_name_value(name_value, tree_dict)
            elif section == 1: # Scene extents
                if extent_y == 0 and re.match('\s*\|\s*\|\s*', line):
                    extent_y += 1
                    continue
                else:
                    for extent_x in range(2): # Read in L & R values
                        name_value = get_name_value(line[extent_x * 27: (extent_x + 1) * 27 - 1])
                        if name_value:
                            # Prefix UL, UR, LL, LR to name
                            name_value = (y_prefix[extent_y] + x_prefix[extent_x] + '_' + name_value[0],
                                          name_value[1])
                            add_name_value(name_value, tree_dict)


    def read_file(self, filename=None):
        """Function to parse a Report metadata file and store the results in self._metadata_dict
        Argument:
            filename: Report Metadata file to be parsed and stored
        Returns:
            nested dict containing metadata
        """
        logger.debug('read_file(%s) called', filename)

        filename = filename or self._filename
        assert filename, 'Filename must be specified'

        # Open Report document
        try:
            infile = open(filename, 'r')
            assert infile is not None, 'Unable to open Report file ' + filename
            return self.parse_report_string(infile.read())
        finally:
            infile.close()

        return self._metadata_dict



    def write_file(self, filename=None, save_backup=False):
        """Function write the metadata contained in self._metadata_dict to an Report file
        Argument:
            filename: Metadata file to be written
        """
        logger.debug('write_file(%s) called', filename)

        filename = filename or self._filename
        assert filename, 'Filename must be specified'

        if save_backup and os.path.exists(filename + '.bck'):
            os.remove(filename + '.bck')

        if os.path.exists(filename):
            if save_backup:
                os.rename(filename, filename + '.bck')
            else:
                os.remove(filename)

        # Open Report document
        try:
            outfile = open(filename, 'w')
            assert outfile is not None, 'Unable to open Report file ' + filename + ' for writing'

            logger.debug('Writing Report file %s', filename)

            for name, value in self.tree_to_tuples():
                self._value_dict[name] = value

            # Convert "NUP" orientation to "0.0 deg"
            if 'ORIENTATION' in self._value_dict.keys() and self._value_dict['ORIENTATION'] in ReportMetadata._orient_map.keys():
                self._value_dict['ORIENTATION'] = '%.1f deg' % ReportMetadata._orient_map[self._value_dict['ORIENTATION']]

            outfile.write(ReportMetadata._template_text % self._value_dict)
        finally:
            outfile.close()


def main():
    r = ReportMetadata()
    r.read_file('/home/alex/nbar/test_data/LS7_ETM_OTH_P51_GALPGS01_092_085_20100315/scene01/report.txt')

if __name__ == '__main__':
    main()

