'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)
Some code adapted from NbarProcessor.py and main.py

Write metadata to NBAR output directory
N.B: Requires specified nbar template metadata XML file in same directory
'''
# As at 26/7/12, old code has been imported and partially cleaned up, but not tested

import logging, os, re
from copy import copy
from datetime import datetime

from ULA3 import DataManager
from ULA3.image_processor import ProcessorConfig
from ULA3.dataset import SceneDataset
from ULA3.metadata import Metadata, XMLMetadata, ReportMetadata
from ULA3.utils import log_multiline, execute
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

    write_metadata(DATA, CONFIG, os.path.join(os.path.dirname(__file__), CONFIG.NBAR_XML_METADATA_TEMPLATE), CONFIG.NBAR_PROCESSOR_VERSION)

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

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve input scene dataset'
    logger.debug('SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    nbar_dataset_id = DATA.get_item('nbar_dataset_id.dat', str)
    assert nbar_dataset_id, 'Unable to retrieve nbar_dataset_id string'
    logger.debug('string for nbar_dataset_id retrieved')

    nbar_temp_output = DATA.get_item('nbar_temp_output.dat', str)
    assert nbar_temp_output, 'Unable to retrieve nbar_temp_output string'
    logger.debug('string for nbar_temp_output retrieved')

    # Retrieve all required ancillary data values for lineage string
    ozone_data = DATA.get_item('ozone.dat', dict)
    assert ozone_data, 'Unable to retrieve dict for ozone.dat'
    water_vapour_data = DATA.get_item('water_vapour.dat', dict)
    assert water_vapour_data, 'Unable to retrieve dict for water_vapour.dat'
    aerosol_data = DATA.get_item('aerosol.dat', dict)
    assert aerosol_data, 'Unable to retrieve dict for aerosol.dat'
    elevation_data = DATA.get_item('elevation.dat', dict)
    assert elevation_data, 'Unable to retrieve dict for elevation.dat'
    brdf_data = DATA.get_item('brdf.dat', dict)
    assert brdf_data, 'Unable to retrieve dict for brdf.dat'
    solar_dist_data = DATA.get_item('solar_dist.dat', dict)
    assert solar_dist_data, 'Unable to retrieve dict for solar_dist.dat'
    solar_irrad_data = DATA.get_item('solar_irrad.dat', dict)
    assert solar_irrad_data, 'Unable to retrieve dict for solar_irrad.dat'
    thumbnail_data = DATA.get_item('thumbnail.dat', dict)
    assert thumbnail_data, 'Unable to retrieve dict for thumbnail.dat'
    new_metadata = Metadata() # Create new master metadata object

    # Open template metadata.xml
    xml_metadata = XMLMetadata(xml_metadata_template)
    # Set reference in master Metadata object
    new_metadata.set_root_metadata(xml_metadata.get_metadata(), 'XML')
    # N.B: ALL template XML metadata values will be defined either as fixed or empty strings

    report_metadata = ReportMetadata()

    # Convert all None values to empty strings for report
    report_metadata_dict = report_metadata.get_metadata()
    for key in report_metadata_dict:
        if report_metadata_dict[key] is None:
            report_metadata_dict[key] = ''

    # Set reference in master Metadata object
    new_metadata.set_root_metadata(report_metadata_dict, 'REPORT')

    if CONFIG.debug:
        log_multiline(logger.info, new_metadata._metadata_dict, 'Template metadata' '\t')

    # Copy any existing metadata from source dataset. Keep any non-empty data from template.
    new_metadata.merge_metadata_dicts(source_tree=l1t_input_dataset.GetMetadata('XML,EODS_DATASET'),
                                  destination_tree=new_metadata.get_metadata('XML,EODS_DATASET'),
                                  overwrite=True, add_new_nodes=False, keep_existing_data=True)
    new_metadata.merge_metadata_dicts(source_tree=l1t_input_dataset.GetMetadata('REPORT'),
                                  destination_tree=new_metadata.get_metadata('REPORT'),
                                  overwrite=True, add_new_nodes=False, keep_existing_data=True)
    if CONFIG.debug:
        log_multiline(logger.info, new_metadata._metadata_dict, 'Post-merge metadata' '\t')

    nbar_dataset = copy(l1t_input_dataset) # Create a shallow copy of input dataset so we can safely change instance values
    nbar_dataset.set_metadata_object(new_metadata) # Change scene dataset metadata to newly created object

    # Read in all non-empty values from template, overwriting any old values
    nbar_dataset.read_metadata(ignore_empty_values=True)

    if CONFIG.debug:
        log_multiline(logger.info, nbar_dataset.__dict__, 'Scene Dataset after reading template metadata' '\t')

    # Update creation time if files have been re-written. Overwrite any values in template
    completion_datetime = DATA.get_item('create_datetime', datetime)
    if not completion_datetime:
        completion_datetime = datetime.utcnow()
    if nbar_dataset.completion_datetime != completion_datetime:
        nbar_dataset.completion_datetime = completion_datetime
        nbar_dataset.completion_date = nbar_dataset.completion_datetime.date()
        nbar_dataset.completion_time = nbar_dataset.completion_datetime.time()

    # Compose lineage statement here
    # Collate all BRDF info
    brdf_descs = []
    for band in sorted(l1t_input_dataset.bands('REFLECTIVE')):
        brdf_descs.append('BAND %d [%.6f (%s), %.6f (%s), %.6f (%s)]' % (
            band,
            brdf_data[(band, 'geo')]['value'], brdf_data[(band, 'geo')]['data_file'],
            brdf_data[(band, 'iso')]['value'], brdf_data[(band, 'iso')]['data_file'],
            brdf_data[(band, 'vol')]['value'], brdf_data[(band, 'vol')]['data_file']
            )
        )

    nbar_dataset.lineage_statement = '\n'.join([
        "NBAR version %s" % PROCESSOR_VERSION,
        'BRDF (%s) = %s' % (CONFIG.BRDF_TYPE, '\n\t'.join(brdf_descs)),
        'OZONE (%s) = %s' % ('GA 2 deg pixel size', str(ozone_data['value'])),
        'WATER VAPOUR (%s) = %s' % ('NOAA', str(water_vapour_data['value'])),
        'DEM (%s) = %s' % ('GA 1 deg pixel size', str(elevation_data['value'])),
        'AOD (%s) = %s' % (aerosol_data['data_file'], str(aerosol_data['value']))
        ])

    # Update instance values with new values
    nbar_dataset.dataset_id = nbar_dataset_id
    nbar_dataset.processor_level = 'NBAR'
    nbar_dataset.product_format = 'GEOTIFF'

    nbar_dataset.mission_name = ''; # This should be blanked out

    nbar_dataset.abstract = 'Processed Image: Landsat Nadir BRDF-Adjusted Reflectance products ' \
                            'generated by Geoscience Australia. Data Source: GA'

    nbar_dataset.algorithm_version = 'Git version ' + CONFIG.git_version
    nbar_dataset.algorithm_title = 'NBAR'

    # Generate comma-separated list of reflective bands
    available_bands = [nbar_dataset.sensor_band_info(band_number)['NUMBER'] // 10 for band_number in nbar_dataset.bands('REFLECTIVE')]
    if any(value == 0 for value in available_bands): # To escape the fact that a potential band 13 might be used
        available_bands = [nbar_dataset.sensor_band_info(band_number)['NUMBER'] for band_number in nbar_dataset.bands('REFLECTIVE')]
    nbar_dataset.available_bands = ', '.join([str(band_number) for band_number in available_bands])
    #nbar_dataset.available_bands = ', '.join([str(nbar_dataset.sensor_band_info(band_number)['NUMBER'] // 10)
    #    for band_number in nbar_dataset.bands('REFLECTIVE')])

    nbar_dataset.browsegraphic_filename = os.path.basename(thumbnail_data['filename'])
#    nbar_dataset.browsegraphic_description = "Low-res scene preview image"
    nbar_dataset.browsegraphic_type = "JPG"
    nbar_dataset.browsegraphic_resolution = nbar_dataset.pixel_x_size * nbar_dataset.image_pixels / thumbnail_data['size'][0]

    def format_band_num(b):
        """
        Format as string, where zero signifies no value
        @type b: int
        @rtype: str
        """
        return '' if b == 0 else '%d' % b

    nbar_dataset.browsegraphic_redband = format_band_num(thumbnail_data['rgb_bands'][0])
    nbar_dataset.browsegraphic_greenband = format_band_num(thumbnail_data['rgb_bands'][1])
    nbar_dataset.browsegraphic_blueband = format_band_num(thumbnail_data['rgb_bands'][2])

    nbar_dataset.alternate_title = None
    nbar_dataset.dq_measure_name = None
    nbar_dataset.dq_quantitative_value = None
    nbar_dataset.dq_quantitative_value_unit = None

    nbar_dataset.environment_description = CONFIG.ENVIRONMENT
    if nbar_dataset.environment_description:
        nbar_dataset.environment_description += ' | '
    nbar_dataset.environment_description += re.sub('\n', '', execute(command_string='uname -a')['stdout'])

    # Determine size of root geoTIFF file
    #root_tif_filename = os.path.join(nbar_temp_output, 'scene01', '%s_B%2d%s' % (nbar_dataset_id, nbar_dataset.satellite.root_band, '.tif'))
    root_tif_filename = os.path.join(nbar_temp_output, 'scene01', '%s_B%d%s' % (nbar_dataset_id, nbar_dataset.satellite.root_band, '.tif'))
    nbar_dataset.file_size = int(round(getFileSizeMB(root_tif_filename)))

    # Use command-line arguments for production parameters if they are defined
    if CONFIG.constraint_id is not None:
        nbar_dataset.constraint_id = int(CONFIG.constraint_id)
    if CONFIG.li_source_description is not None:
        nbar_dataset.li_source_description = CONFIG.li_source_description
    if CONFIG.purpose is not None:
        nbar_dataset.purpose = CONFIG.purpose

    nbar_dataset.title = '%s %s NBAR x%03d y%03d %s version 1 status completed' % (
        nbar_dataset.satellite.NAME,
        re.sub('\W+', '', nbar_dataset.satellite.sensor),
        nbar_dataset.path_number,
        nbar_dataset.row_number,
        nbar_dataset.scene_centre_datetime.strftime('%Y%m%d')
        )

    if CONFIG.debug:
        log_multiline(logger.info, nbar_dataset.__dict__, 'Scene Dataset before updating metadata' '\t')

    # Update output metadata values from instance values. Do NOT extend template metadata
    nbar_dataset.update_metadata(create_new_nodes=False)

    #===========================================================================
    # # Change any specific values not defined in scene_dataset.xml here
    # metadata.set_metadata_node('XML,EODS_DATASET,QUICKLOOK,QL_RED', str(nbar_dataset.satellite.rgb_bands[0] // 10))
    # metadata.set_metadata_node('XML,EODS_DATASET,QUICKLOOK,QL_GREEN', str(nbar_dataset.satellite.rgb_bands[1] // 10))
    # metadata.set_metadata_node('XML,EODS_DATASET,QUICKLOOK,QL_BLUE', str(nbar_dataset.satellite.rgb_bands[2] // 10))
    #
    # metadata.set_metadata_node('REPORT,OUTPUT_BANDS', available_bands + ' (' +
    #                           metadata.get_metadata('REPORT,PIXEL_SIZE') +
    #                           ')') #TODO: Fill in pixel size in degrees
    #===========================================================================

    log_multiline(logger.info, new_metadata._metadata_dict, 'Metadata to be output for NBAR' '\t')

    xml_metadata.write_file(os.path.join(nbar_temp_output, 'metadata.xml'))
    report_metadata.write_file(os.path.join(nbar_temp_output, 'scene01', 'report.txt'))






