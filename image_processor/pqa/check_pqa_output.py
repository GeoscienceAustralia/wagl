'''
Created on Jun 14, 2012

@author: Alex Ip (alex.ip@ga.gov.au)

Determines then sets pqa_dataset_id in data manager
Checks whether repackaging and/or reprocessing required.
Some code adapted from old python/NbarProcessor.py

N.B: This module is a special case. The only sys.exit() call in any of the
image_processor sub-modules occurs here
'''
import logging, os, re, errno, sys
from datetime import datetime
from glob import glob
from ULA3 import DataManager
from ULA3.dataset import SceneDataset
from ULA3.metadata import XMLMetadata
from ULA3.utils import find_files
from ULA3.image_processor import ProcessorConfig
from ULA3.utils import log_multiline, execute

# This import is a special case required to allow repackaging
import package_pqa_output

logger = logging.getLogger('root.' + __name__)

def process(subprocess_list=[], resume=False):
    logger.info('%s.process(%s, %s) called', __name__, subprocess_list, resume)

    CONFIG = ProcessorConfig()
    DATA = DataManager()

    l1t_input_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert l1t_input_dataset, 'Unable to retrieve input scene dataset'
    logger.debug('SceneDataset object for %s retrieved', l1t_input_dataset.pathname)

    # Define these here so that they remain in scope
    pqa_temp_output = None

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

    def backup_output_path(pqa_output_path, work_path):
        command_string = 'mv -f %s %s/backup' % (pqa_output_path, work_path)

        logger.info('Invoking: %s', command_string)

        result = execute(command_string=command_string,
            cwd=CONFIG.work_path
            )

        if result['stdout']:
            log_multiline(logger.info, result['stdout'], command_string + ' in ' + pqa_temp_output, '\t')

        if result['returncode']:
            log_multiline(logger.error, result['stderr'], 'stderr from ' + command_string, '\t')
            raise Exception('%s failed', command_string)


    def check_dataset_id():
        logger.debug('  check_dataset_id(%s) called')

        def get_pqa_dataset_id():
            station_code = CONFIG.station_code.get(l1t_input_dataset.ground_station)
            assert station_code, 'Unrecognised ground station ID %s' % l1t_input_dataset.ground_station

            # Use ID provided or derive new dataset ID for PQA output
            return CONFIG.dataset_id or '%s_%s_PQ_%s_%s-%s_%03d_%03d_%s' % (
                l1t_input_dataset.satellite.TAG,
                re.sub('\W', '', l1t_input_dataset.sensor), # Strip out any non-alphanumeric chars (ETM+ -> ETM)
                CONFIG.PQA_PRODUCT_CODE,
                CONFIG.PQA_PRODUCT_SPEC,
                station_code,
                l1t_input_dataset.path_number,
                l1t_input_dataset.row_number,
                l1t_input_dataset.scene_centre_date.strftime('%Y%m%d'))

        def rename_dataset(pqa_output_path, old_dataset_id, new_dataset_id):
            new_pqa_output_path = pqa_output_path.replace(old_dataset_id, new_dataset_id)
            os.rename(pqa_output_path, new_pqa_output_path)
            logger.info('Directory %s renamed to %s', pqa_output_path, new_pqa_output_path)
            pqa_output_path = new_pqa_output_path

            for filename in find_files(root_dir=pqa_output_path,
                       filename_pattern=old_dataset_id + '.*',
                       ):
                new_filename = filename.replace(old_dataset_id, new_dataset_id)
                os.rename(filename, new_filename)
                logger.info('File %s renamed to %s', filename, new_filename)


        pqa_dataset_id = get_pqa_dataset_id()
        logger.info('pqa_dataset_id: %s', pqa_dataset_id)
        DATA.set_item('pqa_dataset_id.dat', pqa_dataset_id)

        # Determine path of temporary output image and create directory
        pqa_temp_output = os.path.join(CONFIG.work_path, pqa_dataset_id)
        logger.info('pqa_temp_output: %s', pqa_temp_output)
        DATA.set_item('pqa_temp_output.dat', pqa_temp_output)

        # Default output_path to standard directory under PQA_DATA_ROOT if not specified
        pqa_output_path = os.path.abspath(CONFIG.output_path or os.path.join(CONFIG.output_root, pqa_dataset_id))
        logger.info('pqa_output_path: %s', pqa_output_path)
        DATA.set_item('pqa_output_path.dat', pqa_output_path)

        old_dataset_id = re.match('.*(?=_)', os.path.basename(l1t_input_dataset.root_dataset_pathname)).group()
        old_pqa_tif_template = os.path.join(pqa_output_path, 'scene01', '%s_*.tif' % old_dataset_id)
        final_pqa_tif_template = os.path.join(pqa_output_path, 'scene01', '%s_*.tif' % pqa_dataset_id)
        logger.debug('old_dataset_id = %s, old_pqa_tif_template = %s, final_pqa_tif_template = %s', old_dataset_id, old_pqa_tif_template, final_pqa_tif_template)

        final_pqa_tif_filename = None
        filelist = glob(old_pqa_tif_template) or glob(re.sub('.tif$', '.TIF$', old_pqa_tif_template))
        if filelist: # Old format PQA output filename found
            old_pqa_tif_filename = filelist[0] # There should be only one
            m = re.match('(.*/)(.*)(_\d+\.\w+)', old_pqa_tif_filename) # Path, dataset_id and suffix+extension
            final_pqa_tif_filename = m.group(1) + pqa_dataset_id + m.group(3)
            logger.warning('Renaming existing output file from %s to %s', old_pqa_tif_filename, final_pqa_tif_filename)
            os.rename(old_pqa_tif_filename, final_pqa_tif_filename)
        else:
            filelist = glob(final_pqa_tif_template) or glob(re.sub('.tif$', '.TIF$', final_pqa_tif_template))
            if filelist:
                final_pqa_tif_filename = filelist[0] # There should be only one
            else:
                logger.debug('No existing output file found')
        del filelist

        if final_pqa_tif_filename: # Final output already exists
            logger.debug('Output file %s already exists' % final_pqa_tif_filename)

            # Use either MTL metadata value or file modification date for input/output datasets
            logger.debug('Input modification datetime (from metadata) = %s', l1t_input_dataset.completion_datetime)
            input_datetime = l1t_input_dataset.completion_datetime or datetime.fromtimestamp(
                os.path.getmtime(l1t_input_dataset.root_dataset_pathname))
            logger.debug('Input modification datetime = %s', input_datetime)

            # Attempt to determine output completion time if "final" XML metadata file exists
            output_datetime = None
            if os.path.exists(os.path.join(pqa_output_path, 'metadata.xml')):
                xml = XMLMetadata(os.path.join(pqa_output_path, 'metadata.xml'))
                output_datetime = (get_datetime(xml.get_metadata('METADATA,CITATION,CITATION_DATE'))
                                   or
                                   get_datetime(xml.get_metadata('EODS_DATASET,MDRESOURCE,CITATION,DATE'))
                                   )
                logger.debug('Output completion datetime (from existing XML file) = %s', output_datetime)
            output_datetime = output_datetime or datetime.fromtimestamp(os.path.getmtime(
                final_pqa_tif_filename))
            logger.debug('Output completion datetime = %s', output_datetime)

            # Output file is newer than input - no re-processing required
            if output_datetime > input_datetime:
                logger.debug('%s is newer than %s' % (final_pqa_tif_filename, l1t_input_dataset.root_dataset_pathname))

                # This is slightly dodgy but it will allow us to call this routine to re-package existing output
                pqa_temp_output = pqa_output_path
                DATA.set_item('pqa_temp_output.dat', pqa_temp_output)
                logger.debug('pqa_temp_output set to %s', pqa_temp_output)

                DATA.set_item('create_datetime', output_datetime) # Keep previous creation datetime

                if CONFIG.repackage: # Repackaging required
                    logger.debug('%s is being repackaged', l1t_input_dataset.root_dataset_pathname)

                    package_pqa_output.process() # Skip straight to packaging

                    logger.warning('Repackaging completed successfully')

                if not CONFIG.subprocess and not CONFIG.subprocess_file:
                    logger.warning('Skipped re-processing of unchanged output')
                    logger.debug('Abandoning further processing')
                    sys.exit(0) # Abort processing
            else:
                # Dataset must be re-processed completely. Move existing output to work directory.
                logger.warning('Moving out-of-date result dataset %s to %s/backup', pqa_output_path, CONFIG.work_path)
                backup_output_path(pqa_output_path, CONFIG.work_path)


        try:
            os.makedirs(os.path.join(pqa_temp_output, 'scene01'))
            logger.info('Created temporary PQA output directory %s with sub-directory scene01', pqa_temp_output)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise

    check_dataset_id()
