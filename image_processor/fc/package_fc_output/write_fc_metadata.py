'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
Some code adapted from NbarProcessor.py and main.py

Write metadata to FC output directory
N.B: Requires specified fc template metadata XML file in same directory
'''
# As at 26/7/12, old code has been imported and partially cleaned up, but not tested

import logging, os, re
from copy import copy
from datetime import datetime

from ULA3 import DataManager
from ULA3.image_processor import ProcessorConfig
from ULA3.dataset import SceneDataset
from ULA3.metadata import Metadata
from ULA3.metadata import XMLMetadata
#from ULA3.metadata import ReportMetadata
from ULA3.utils import log_multiline, execute, find_files
from ULA3.meta import print_call

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

    write_metadata(DATA, CONFIG, os.path.join(os.path.dirname(__file__), CONFIG.FC_XML_METADATA_TEMPLATE), CONFIG.FC_PROCESSOR_VERSION)

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

    nbar_input_dataset = DATA.get_item(CONFIG.input['nbar']['path'], SceneDataset)
    assert nbar_input_dataset, 'Unable to retrieve input scene dataset'
    logger.debug('SceneDataset object for %s retrieved', nbar_input_dataset.pathname)

    fc_dataset_id = DATA.get_item('fc_dataset_id.dat', str)
    assert fc_dataset_id, 'Unable to retrieve fc_dataset_id string'
    logger.debug('string for fc_dataset_id retrieved')

    fc_temp_output = DATA.get_item('fc_temp_output.dat', str)
    assert fc_temp_output, 'Unable to retrieve fc_temp_output string'
    logger.debug('string for fc_temp_output retrieved')

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
    new_metadata.merge_metadata_dicts(source_tree=nbar_input_dataset.GetMetadata('XML,EODS_DATASET'),
                                  destination_tree=new_metadata.get_metadata('XML,EODS_DATASET'),
                                  overwrite=True, add_new_nodes=False, keep_existing_data=True)
#===============================================================================
#     new_metadata.merge_metadata_dicts(source_tree=nbar_input_dataset.GetMetadata('REPORT'),
#                                   destination_tree=new_metadata.get_metadata('REPORT'),
#                                   overwrite=True, add_new_nodes=False, keep_existing_data=True)
#===============================================================================
    if CONFIG.debug:
        log_multiline(logger.info, new_metadata._metadata_dict, 'Post-merge metadata' '\t')

    fc_dataset = copy(nbar_input_dataset) # Create a shallow copy of input dataset so we can safely change instance values
    fc_dataset.set_metadata_object(new_metadata) # Change scene dataset metadata to newly created object

    # Read in all non-empty values from template, overwriting any old values
    fc_dataset.read_metadata(ignore_empty_values=True)

    if CONFIG.debug:
        log_multiline(logger.info, fc_dataset.__dict__, 'Scene Dataset after reading template metadata' '\t')

    # Update creation time if files have been re-written. Overwrite any values in template
    completion_datetime = DATA.get_item('create_datetime', datetime)
    if not completion_datetime:
        completion_datetime = datetime.utcnow()
    if fc_dataset.completion_datetime != completion_datetime:
        fc_dataset.completion_datetime = completion_datetime
        fc_dataset.completion_date = fc_dataset.completion_datetime.date()
        fc_dataset.completion_time = fc_dataset.completion_datetime.time()

    # Update instance values with new values
    fc_dataset.dataset_id = fc_dataset_id
    fc_dataset.processor_level = 'Fractional Cover'
    fc_dataset.product_format = 'GEOTIFF'

    fc_dataset.mission_name = ''; # This should be blanked out

    fc_dataset.abstract = '%s FC' % fc_dataset.satellite.NAME
    fc_dataset.available_bands = 'FC'

    # Remove values pulled through from L1
    fc_dataset.band_gain = None
    fc_dataset.browsegraphic_filename = fc_dataset_id + '.jpg'

    fc_dataset.alternate_title = None
    fc_dataset.dq_measure_name = None
    fc_dataset.dq_quantitative_value = None
    fc_dataset.dq_quantitative_value_unit = None

    fc_dataset.environment_description = CONFIG.ENVIRONMENT
    if fc_dataset.environment_description:
        fc_dataset.environment_description += ' | '
    fc_dataset.environment_description += re.sub('\n', '', execute(command_string='uname -a')['stdout'])

    fc_dataset.algorithm_version = 'SVN version ' + CONFIG.svn_revision

#===============================================================================
#     fc_dataset.supplementary_information = """
#
# """.format()
#     log_multiline(logger.debug, fc_dataset.supplementary_information, 'supplementary_information', '\t')
#===============================================================================

    # Determine size of root geoTIFF file
    root_tif_filename = find_files(os.path.join(fc_temp_output, 'scene01'), 'L.*.tif')[-1] #TODO: Make this a bit more robust
    fc_dataset.file_size = int(round(getFileSizeMB(root_tif_filename)))

    # Use command-line arguments for production parameters if they are defined
    if CONFIG.constraint_id is not None:
        fc_dataset.constraint_id = int(CONFIG.constraint_id)
    if CONFIG.li_source_description is not None:
        fc_dataset.li_source_description = CONFIG.li_source_description
    if CONFIG.purpose is not None:
        fc_dataset.purpose = CONFIG.purpose

    fc_dataset.title = '%s %s FC x%03d y%03d %s version 0 status completed' % (
        fc_dataset.satellite.NAME,
        re.sub('\W+', '', fc_dataset.satellite.sensor),
        fc_dataset.path_number,
        fc_dataset.row_number,
        fc_dataset.scene_centre_datetime.strftime('%Y%m%d')
        )

    if CONFIG.debug:
        log_multiline(logger.info, fc_dataset.__dict__, 'Scene Dataset before updating metadata' '\t')

    # Update output metadata values from instance values. Do NOT extend template metadata
    fc_dataset.update_metadata(create_new_nodes=False)

    log_multiline(logger.info, new_metadata._metadata_dict, 'Metadata to be output for FC' '\t')

    xml_metadata.write_file(os.path.join(fc_temp_output, 'metadata.xml'))
#    report_metadata.write_file(os.path.join(fc_temp_output, 'scene01', 'report.txt'))






