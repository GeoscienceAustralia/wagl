'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)

Determines then sets nbar_dataset_id in data manager
Checks whether repackaging and/or reprocessing required.
Some code adapted from old python/NbarProcessor.py

This module also implements a special case where ancillary data
failure will trigger the removal of an existing output dataset. This is to

There is a bit of ugly stuff in here to cater for datasets processed by the
previous ULA NBAR. These do not have the ground station ID in their names

N.B: This module is a special case. The only sys.exit() call in any of the
image_processor sub-modules occurs here
'''
import logging, os, re, errno, sys
from datetime import datetime
from ULA3 import DataManager
from ULA3.dataset import SceneDataset
from ULA3.metadata import XMLMetadata
from ULA3.utils import find_files, log_multiline
from ULA3.image_processor import ProcessorConfig
from ULA3.utils import execute

# This import is a special case required to allow repackaging
import get_ancillary_data
import package_nbar_output

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    def ga_product_spec(use_station_id=True):
        """Get product specification string.

        Parameter:
            use_station_id: Boolean flag indicating whether station code should be included

        Returns:
            Product specification string
        """

    #===========================================================================
    #    if use_station_id:
    #        m = re.search(CONFIG.INPUT_RE_PATTERN, os.path.basename(CONFIG.input['l1t']['path']))
    #
    #        if m:
    #            try:
    #                _input_label, station_code = m.group(5).split('-')
    #                return '%s-%s' % (CONFIG.NBAR_PRODUCT_SPEC, station_code)
    #            except ValueError:
    #                pass
    #===========================================================================
        if use_station_id:
            station_code = CONFIG.station_code.get(l1t_input_dataset.ground_station)
            if station_code:
                return '%s-%s' % (CONFIG.NBAR_PRODUCT_SPEC, station_code)
            else:
                logger.warning('Unrecognised ground station ID %s', l1t_input_dataset.ground_station)

        return CONFIG.NBAR_PRODUCT_SPEC

    def get_datetime(datetime_string):
        """Returns datetime parsed from string
        Needs to cater for different variations in datetime format
        """
        output_datetime = None
        if datetime_string:
            s = re.search('(\d{4})-?(\d{2})-?(\d{2})(T|\s)(\d{2}:\d{2}:\d{2})', datetime_string)
            if s:
                try:
                    output_datetime = datetime.strptime(
                        s.group(1) + '-' + s.group(2) + '-' + s.group(3) + ' ' + s.group(5),
                        '%Y-%m-%d %H:%M:%S')
                except:
                    pass
        return output_datetime

    def backup_output_path(nbar_output_path, work_path):
        command_string = 'mv -f %s %s/backup' % (nbar_output_path, work_path)

        logger.info('Invoking: %s', command_string)

        result = execute(command_string=command_string,
            cwd=CONFIG.work_path
            )

        if result['stdout']:
            log_multiline(logger.info, result['stdout'], command_string + ' in ' + nbar_temp_output, '\t')

        if result['returncode']:
            log_multiline(logger.error, result['stderr'], 'stderr from ' + command_string, '\t')
            raise Exception('%s failed', command_string)


    def check_dataset_id(use_station_id=True):
        logger.info('  check_dataset_id(%s) called', use_station_id)

        def get_nbar_dataset_id(use_station_id=True):
            # Use ID provided or derive new dataset ID for NBAR output
            return CONFIG.dataset_id or '%s_%s_NBAR_%s_%s_%03d_%03d_%s' % (
                l1t_input_dataset.satellite.TAG,
                re.sub('\W', '', l1t_input_dataset.sensor), # Strip out any non-alphanumeric chars (ETM+ -> ETM)
                CONFIG.NBAR_PRODUCT_CODE,
                ga_product_spec(use_station_id),
                l1t_input_dataset.path_number,
                l1t_input_dataset.row_number,
                l1t_input_dataset.scene_centre_date.strftime('%Y%m%d'))

        def rename_dataset(nbar_output_path, old_dataset_id, new_dataset_id):
            new_nbar_output_path = nbar_output_path.replace(old_dataset_id, new_dataset_id)
            os.rename(nbar_output_path, new_nbar_output_path)
            logger.info('Directory %s renamed to %s', nbar_output_path, new_nbar_output_path)
            nbar_output_path = new_nbar_output_path

            for filename in find_files(root_dir=nbar_output_path,
                       filename_pattern=old_dataset_id + '.*',
                       ):
                new_filename = filename.replace(old_dataset_id, new_dataset_id)
                os.rename(filename, new_filename)
                logger.info('File %s renamed to %s', filename, new_filename)


        nbar_dataset_id = get_nbar_dataset_id(use_station_id)
        logger.info('nbar_dataset_id: %s', nbar_dataset_id)
        DATA.set_item('nbar_dataset_id.dat', nbar_dataset_id)

        # Determine path of temporary output image and create directory
        nbar_temp_output = os.path.join(CONFIG.work_path, nbar_dataset_id)
        logger.info('nbar_temp_output: %s', nbar_temp_output)
        DATA.set_item('nbar_temp_output.dat', nbar_temp_output)

        # Default output_path to standard directory under NBAR_DATA_ROOT if not specified
        nbar_output_path = os.path.abspath(CONFIG.output_path or os.path.join(CONFIG.output_root, nbar_dataset_id))
        logger.info('nbar_output_path: %s', nbar_output_path)
        DATA.set_item('nbar_output_path.dat', nbar_output_path)

        band_file_number = l1t_input_dataset.sensor_band_info(l1t_input_dataset.root_band_number)['NUMBER']
        final_root_tif_filename = os.path.join(nbar_output_path, 'scene01', '%s_B%2d%s' % (nbar_dataset_id, band_file_number, '.tif'))

        if os.path.exists(final_root_tif_filename): # Final output already exists
            logger.debug('Output file %s already exists' % final_root_tif_filename)

            # Use either MTL metadata value or file modification date for input/output datasets
            logger.debug('Input modification datetime (from metadata) = %s', l1t_input_dataset.completion_datetime)
            input_datetime = l1t_input_dataset.completion_datetime or datetime.fromtimestamp(
                os.path.getmtime(l1t_input_dataset.root_dataset_pathname))
            logger.debug('Input modification datetime = %s', input_datetime)

            # Attempt to determine output completion time if "final" XML metadata file exists
            output_datetime = None
            if os.path.exists(os.path.join(nbar_output_path, 'metadata.xml')):
                xml = XMLMetadata(os.path.join(nbar_output_path, 'metadata.xml'))
                output_datetime = (get_datetime(xml.get_metadata('METADATA,CITATION,CITATION_DATE'))
                                   or
                                   get_datetime(xml.get_metadata('EODS_DATASET,MDRESOURCE,CITATION,DATE'))
                                   )
                logger.debug('Output completion datetime (from existing XML file) = %s', output_datetime)
            output_datetime = output_datetime or datetime.fromtimestamp(os.path.getmtime(
                final_root_tif_filename))
            logger.debug('Output completion datetime = %s', output_datetime)

            # Output file is newer than input - no re-processing required
            if output_datetime > input_datetime:
                logger.debug('%s is newer than %s' % (final_root_tif_filename, l1t_input_dataset.root_dataset_pathname))

                if use_station_id: # We are good to go with the current dataset ID
                    # This is slightly dodgy but it will allow us to call this routine to re-package existing output
                    nbar_temp_output = nbar_output_path
                    DATA.set_item('nbar_temp_output.dat', nbar_temp_output)
                    logger.debug('nbar_temp_output set to %s', nbar_temp_output)

                    DATA.set_item('create_datetime', output_datetime) # Keep previous creation datetime

                    if CONFIG.repackage: # Repackaging required
                        logger.debug('%s is being repackaged', l1t_input_dataset.root_dataset_pathname)
                        try:
                            get_ancillary_data.process() # Re-fetch ancillary data required for lineage metadata
                        except (Exception), e:
                            # Remove invalid  dataset from product area - place in work directory
                            # This is required to overcome a bug in ULA2 NBAR where it should have failed on BRDF
                            logger.error('get_ancillary_data failed: %s', e.message)

                            backup_output_path(nbar_output_path, CONFIG.work_path)

                            logger.warning('Invalid output dataset %s moved to %s/backup', nbar_output_path, CONFIG.work_path)
                            raise e

                        package_nbar_output.process() # Skip straight to packaging

                        logger.warning('Repackaging completed successfully')

                    if not CONFIG.subprocess and not CONFIG.subprocess_file:
                        logger.warning('Skipped re-processing of unchanged output')
                        logger.debug('Abandoning further processing')
                        sys.exit(0) # Abort processing
                else: # Found old dataset - Need to rename directory
                    rename_dataset(nbar_output_path, nbar_dataset_id, get_nbar_dataset_id(True))
            else:
                # Dataset must be re-processed completely. Move existing output to work directory.
                logger.warning('Moving out-of-date result dataset %s to %s/backup', nbar_output_path, CONFIG.work_path)
                backup_output_path(nbar_output_path, CONFIG.work_path)


        if use_station_id: # We are good to go with the current dataset ID
            try:
                os.makedirs(os.path.join(nbar_temp_output, 'scene01'))
                logger.info('Created temporary NBAR output directory %s with sub-directory scene01', nbar_temp_output)
            except OSError, e:
                if e.errno != errno.EEXIST:
                    raise
        else: # Now check dataset ID WITH station ID
            check_dataset_id(use_station_id=True)


    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve input scene dataset'
    logger.debug('SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    # Define these here so that they remain in scope
    nbar_dataset_id = None
    nbar_temp_output = None
    nbar_output_path = None

    # Check for existing output dataset first without and then with station ID
    # Don't bother trying twice if dataset_id has been provided as a command line argument
    check_dataset_id(use_station_id=bool(CONFIG.dataset_id))



