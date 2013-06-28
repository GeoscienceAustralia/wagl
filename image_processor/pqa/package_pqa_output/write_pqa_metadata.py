'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
Some code adapted from NbarProcessor.py and main.py

Write metadata to PQA output directory
N.B: Requires specified pqa template metadata XML file in same directory
'''
# As at 26/7/12, old code has been imported and partially cleaned up, but not tested

import logging, os, re
from copy import copy
from datetime import datetime

from ULA3 import DataManager
from ULA3.image_processor import ProcessorConfig
from ULA3.dataset import SceneDataset
from ULA3.metadata import Metadata, XMLMetadata
#from ULA3.metadata import ReportMetadata
from ULA3.utils import log_multiline, execute, find_files
from ULA3.meta import print_call
from pqa_log_extractor import PQALogExtractor

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    """Write metadata files to output scene directory.

    N.B: Order of precedence for metadata values is as follows:
        1. SceneDataset instance values including command-line specified values (e.g. "purpose")
        2. Non-empty values in metadata template file (e.g. "INDIVIDUALNAME")
        3. Values read from input scene metadata
    """
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    write_metadata(DATA, CONFIG, os.path.join(os.path.dirname(__file__), CONFIG.PQA_XML_METADATA_TEMPLATE), CONFIG.PQA_PROCESSOR_VERSION)

# TODO: should check of any other things accessed through CONFIG need to be specific to
# the dataset.
@print_call(logger.info)
def write_metadata(DATA, CONFIG, xml_metadata_template, PROCESSOR_VERSION):
    """Write metadata files to output scene directory.

    N.B: Order of precedence for metadata values is as follows:
        1. SceneDataset instance values including command-line specified values (e.g. "purpose")
        2. Non-empty values in metadata template file (e.g. "INDIVIDUALNAME")
        3. Values read from input scene metadata
    """
    logger.info('write_metadata called')

    def getFileSizeMB(path):
        """Gets the size of a file (megabytes).

        Arguments:
            path: file path

        Returns:
            File size (MB)

        Raises:
            OSError [Errno=2] if file does not exist
        """
        return os.path.getsize(path) / (1024*1024)

    def compose_lineage(l1t_input_dataset, tests_run):
        """
        Composes lineage string from values in conf file
        Arguments:
            l1t_input_dataset: SceneDataset object for L1T dataset
            tests_run: List or string consisting of 16 '0' or '1' characters for 16 different tests
        """
        run_string = {'0': 'Not Run',
                      '1': 'Run'}

        test_const_list = list(CONFIG.pqa_test_consts) # Ordered list of all test constants as defined in processor.conf

        # The following code block is a bit ugly because we only have filenames and log files to work with
        # Remove un-used test constant for bit index 5 depending on whether LS5 or LS7
        # TODO: Make this more dynamic and less arbitrary
        if l1t_input_dataset.satellite.TAG == 'LS5':
            test_const_list.remove('SATURATION_61')
        if l1t_input_dataset.satellite.TAG == 'LS7':
            test_const_list.remove('SATURATION_60')
        test_const_dict = {}
        # Invert mapping to look up constants from indices instead of other way around
        for test_const in test_const_list:
            test_const_dict[CONFIG.pqa_test_index[test_const]] = test_const

        lineage_string = """
Pixel Quality SVN Version %s
The pixel quality algorithm assesses quality aspects such as saturation, band/spectral contiguity, land/sea, cloud and cloud shadow.
""" % CONFIG.svn_revision
        for test_index in range(len(tests_run)):
            lineage_string += '%s (Bit %d): ' % (CONFIG.pqa_test_description[test_const_dict[test_index]], test_index)
            if test_index == 6 and l1t_input_dataset.satellite.TAG == 'LS5': #
                lineage_string += 'Duplicated Band60'
            else: # Normal test
                lineage_string += run_string[tests_run[test_index]] # Run / Not Run
            lineage_string += '\n'

        log_multiline(logger.debug, lineage_string, 'lineage_string', '\t')
        return lineage_string

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve input scene dataset'
    logger.debug('SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    pqa_dataset_id = DATA.get_item('pqa_dataset_id.dat', str)
    assert pqa_dataset_id, 'Unable to retrieve pqa_dataset_id string'
    logger.debug('string for pqa_dataset_id retrieved')

    pqa_temp_output = DATA.get_item('pqa_temp_output.dat', str)
    assert pqa_temp_output, 'Unable to retrieve pqa_temp_output string'
    logger.debug('string for pqa_temp_output retrieved')

    pqa_tif_path = find_files(os.path.join(pqa_temp_output, 'scene01'), 'L.*\.tif')[-1]

    pqa_log_info = PQALogExtractor(pqa_tif_path) # Get various values from log file

    new_metadata = Metadata() # Create new master metadata object

    # Open template metadata.xml
    xml_metadata = XMLMetadata(xml_metadata_template)
    # Set reference in master Metadata object
    new_metadata.set_root_metadata(xml_metadata.get_metadata(), 'XML')
    # N.B: ALL template XML metadata values will be defined either as fixed or empty strings

#===============================================================================
#     report_metadata = ReportMetadata()
#
#     # Convert all None values to empty strings for report
#     report_metadata_dict = report_metadata.get_metadata()
#     for key in report_metadata_dict:
#         if report_metadata_dict[key] is None:
#             report_metadata_dict[key] = ''
#
#     # Set reference in master Metadata object
#     new_metadata.set_root_metadata(report_metadata_dict, 'REPORT')
#===============================================================================

    if CONFIG.debug:
        log_multiline(logger.info, new_metadata._metadata_dict, 'Template metadata' '\t')

    # Copy any existing metadata from source dataset. Keep any non-empty data from template.
    new_metadata.merge_metadata_dicts(source_tree=l1t_input_dataset.GetMetadata('XML,EODS_DATASET'),
                                  destination_tree=new_metadata.get_metadata('XML,EODS_DATASET'),
                                  overwrite=True, add_new_nodes=False, keep_existing_data=True)
#===============================================================================
#     new_metadata.merge_metadata_dicts(source_tree=l1t_input_dataset.GetMetadata('REPORT'),
#                                   destination_tree=new_metadata.get_metadata('REPORT'),
#                                   overwrite=True, add_new_nodes=False, keep_existing_data=True)
#===============================================================================
    if CONFIG.debug:
        log_multiline(logger.info, new_metadata._metadata_dict, 'Post-merge metadata' '\t')

    pqa_dataset = copy(l1t_input_dataset) # Create a shallow copy of input dataset so we can safely change instance values
    pqa_dataset.set_metadata_object(new_metadata) # Change scene dataset metadata to newly created object

    # Read in all non-empty values from template, overwriting any old values
    pqa_dataset.read_metadata(ignore_empty_values=True)

    if CONFIG.debug:
        log_multiline(logger.info, pqa_dataset.__dict__, 'Scene Dataset after reading template metadata' '\t')

    # Update creation time if files have been re-written. Overwrite any values in template
    completion_datetime = DATA.get_item('create_datetime', datetime)
    if not completion_datetime:
        completion_datetime = datetime.utcnow()
    if pqa_dataset.completion_datetime != completion_datetime:
        pqa_dataset.completion_datetime = completion_datetime
        pqa_dataset.completion_date = pqa_dataset.completion_datetime.date()
        pqa_dataset.completion_time = pqa_dataset.completion_datetime.time()

    #===========================================================================
    # pqa_dataset.lineage_statement = """Pixel Quality SVN Version 1228
    #      The pixel quality algorithm assesses quality aspects such as saturation, band/spectral contiguity, land/sea, cloud and cloud shadow.
    #                           Tests Run:
    #                           Saturation Band1: {Run}
    #                           Saturation Band2: {Run}
    #                           Saturation Band3: {Run}
    #                           Saturation Band4: {Run}
    #                           Saturation Band5: {Run}
    #                           Saturation Band61: {Run}
    #                           Saturation Band62: {Duplicated Band61}
    #                           Saturation Band7: {Run}
    #                           Band Contiguity: {Run}
    #                           Land/Sea: {Run}
    #                           ACCA: {Run}
    #                           Fmask (Cloud): {Run}
    #                           Cloud Shadow (ACCA): {Run}
    #                           Cloud Shadow (Fmask): {Run}
    #                           Empty Test: {Not Run}
    #                           Empty Test: {Not Run}
    #
    #                           Bits Set:
    #                           Saturation Band1: 0
    #                           Saturation Band2: 1
    #                           Saturation Band3: 2
    #                           Saturation Band4: 3
    #                           Saturation Band5: 4
    #                           Saturation Band61: 5
    #                           Saturation Band62: 6
    #                           Saturation Band7: 7
    #                           Band Contiguity: 8
    #                           Land/Sea: 9
    #                           ACCA: 10
    #                           Fmask (Cloud): 11
    #                           Cloud Shadow (ACCA): 12
    #                           Cloud Shadow (Fmask): 13
    #                           Empty Test: Not Set
    #                           Empty Test: Not Set"""
    #===========================================================================
    pqa_dataset.lineage_statement = compose_lineage(l1t_input_dataset, pqa_log_info.tests_run)

    # Update instance values with new values
    pqa_dataset.dataset_id = pqa_dataset_id
    pqa_dataset.processor_level = 'Pixel Quality'
    pqa_dataset.product_format = 'GEOTIFF'

    pqa_dataset.mission_name = ''; # This should be blanked out

    pqa_dataset.abstract = '%s PQ' % pqa_dataset.satellite.NAME
    pqa_dataset.available_bands = 'PQ'

    # Remove values pulled through from L1
    pqa_dataset.band_gain = None
    pqa_dataset.browsegraphic_filename = None
    pqa_dataset.browsegraphic_description = None
    pqa_dataset.browsegraphic_type = None
    pqa_dataset.browsegraphic_resolution = None
    pqa_dataset.browsegraphic_redband = None
    pqa_dataset.browsegraphic_greenband = None
    pqa_dataset.browsegraphic_blueband = None

    pqa_dataset.alternate_title = None
    pqa_dataset.dq_measure_name = None
    pqa_dataset.dq_quantitative_value = None
    pqa_dataset.dq_quantitative_value_unit = None

    pqa_dataset.environment_description = CONFIG.ENVIRONMENT
    if pqa_dataset.environment_description:
        pqa_dataset.environment_description += ' | '
    pqa_dataset.environment_description += re.sub('\n', '', execute(command_string='uname -a')['stdout'])

    pqa_dataset.algorithm_version = 'SVN version ' + CONFIG.svn_revision

    pqa_dataset.supplementary_information = """
The pixel quality algorithm uses data from both the L1T (Systematic Terrain Correction) and ARG25 (Australian Reflectance Grid 25m) products.

ACCA cloud cover is reported as a percentage of the entire data grid, while Fmask is reported as a percentage of the valid image data only.
CLOUD COVER PERCENTAGE ACCA {ACCA_PERCENT}
CLOUD COVER PERCENTAGE Fmask {Fmask_PERCENT}

Cloud shadow is reported as a percentage of the entire data grid for both analyses.
CLOUD SHADOW PERCENTAGE ACCA {ACCA_CLOUD_SHADOW_PERCENT}
CLOUD SHADOW PERCENTAGE Fmask {Fmask_CLOUD_SHADOW_PERCENT}
""".format(ACCA_PERCENT=round(pqa_log_info.acca_percent, 2),
           Fmask_PERCENT=round(pqa_log_info.fmask_percent, 2),
           ACCA_CLOUD_SHADOW_PERCENT=round(pqa_log_info.acca_cloud_shadow_percent, 2),
           Fmask_CLOUD_SHADOW_PERCENT=round(pqa_log_info.fmask_cloud_shadow_percent, 2)
           )
    log_multiline(logger.debug, pqa_dataset.supplementary_information, 'supplementary_information', '\t')

    # Determine size of root geoTIFF file
    root_tif_filename = find_files(os.path.join(pqa_temp_output, 'scene01'), 'L.*.tif')[-1] #TODO: Make this a bit more robust
    pqa_dataset.file_size = int(round(getFileSizeMB(root_tif_filename)))

    # Use command-line arguments for production parameters if they are defined
    if CONFIG.constraint_id is not None:
        pqa_dataset.constraint_id = int(CONFIG.constraint_id)
    if CONFIG.li_source_description is not None:
        pqa_dataset.li_source_description = CONFIG.li_source_description
    if CONFIG.purpose is not None:
        pqa_dataset.purpose = CONFIG.purpose

    pqa_dataset.title = '%s %s PQ x%03d y%03d %s version 0 status completed' % (
        pqa_dataset.satellite.NAME,
        re.sub('\W+', '', pqa_dataset.satellite.sensor),
        pqa_dataset.path_number,
        pqa_dataset.row_number,
        pqa_dataset.scene_centre_datetime.strftime('%Y%m%d')
        )

    if CONFIG.debug:
        log_multiline(logger.info, pqa_dataset.__dict__, 'Scene Dataset before updating metadata' '\t')

    # Update output metadata values from instance values. Do NOT extend template metadata
    pqa_dataset.update_metadata(create_new_nodes=False)

    #===========================================================================
    # # Change any specific values not defined in scene_dataset.xml here
    # metadata.set_metadata_node('XML,EODS_DATASET,QUICKLOOK,QL_RED', str(pqa_dataset.satellite.rgb_bands[0] // 10))
    # metadata.set_metadata_node('XML,EODS_DATASET,QUICKLOOK,QL_GREEN', str(pqa_dataset.satellite.rgb_bands[1] // 10))
    # metadata.set_metadata_node('XML,EODS_DATASET,QUICKLOOK,QL_BLUE', str(pqa_dataset.satellite.rgb_bands[2] // 10))
    #
    # metadata.set_metadata_node('REPORT,OUTPUT_BANDS', available_bands + ' (' +
    #                           metadata.get_metadata('REPORT,PIXEL_SIZE') +
    #                           ')') #TODO: Fill in pixel size in degrees
    #===========================================================================

    log_multiline(logger.info, new_metadata._metadata_dict, 'Metadata to be output for PQA' '\t')

    xml_metadata.write_file(os.path.join(pqa_temp_output, 'metadata.xml'))
#    report_metadata.write_file(os.path.join(pqa_temp_output, 'scene01', 'report.txt'))






