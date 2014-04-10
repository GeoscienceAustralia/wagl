"""
Definition of the class SceneDataset.

"""

import os, re, math, logging
import numpy
from glob import glob
from datetime import datetime, date, time
import xml.dom.minidom
from pprint import pprint
from osgeo import gdal, gdalconst, osr

import ULA3.geodesic as geocentric
from ULA3.metadata import Metadata, XMLMetadata, MTLMetadata, ReportMetadata
from ULA3.utils import find_files
from ULA3.geodesic import Satellite
from ULA3.utils import unicode_to_ascii, ImageShape

from . import Dataset, DSException

logger = logging.getLogger('root.' + __name__)
gdal.UseExceptions()

class SceneDataset(Dataset):
    """
    This class is a model for a for multi-band Landsat scenes that mimics a :py:class:`gdal.Dataset`
    as closely as possible but provides information contained in various meta-data sources also, and
    deals with some peculiarities of some of the data used within GA (particularly pecularities
    caused by the use of the LPGS processing system.

    SceneDataset manages the sub-datasets for the individual bands in a scene and abstracts out the
    differences between FST and TIFF files. SceneDataset also uses an associated XML class
    configuration file (scene_dataset.xml) to determine which values to extract from the metadata
    and the search order for finding the values in all available metadata. These values are assigned
    to appropriately named (and typed) instance variables when SceneDataset.Open(<scene_directory_path>)
    is called. The :py:meth:`update_metadata` will push the current value of the instance variables back
    into ALL metadata, so this can be used to sync values across all metadata sources prior to writing
    out particular metadata files.

    N.B: Dataset-wide methods will currently return values correct for reflective data only. Thermal
    and pan-chromatic data is managed but ignored for dataset-wide methods and properties
    (except for RasterCount) due to the differing resolution of these bands. Band-specific values
    can be accessed for all bands including thermal and pan-chromatic.

    :todo:
        Enable resampling of thermal and pan-chromatic bands with different native resolutions to the
        reflective bands. Note that this is NOT required for the current NBAR processor, but it may
        be required later. This could be implemented using a VRT dataset, but multi-band FST datasets
        would have to be referenced via one additonal VRT dataset per band, resulting in a nested VRT
        dataset.

    :todo:
        One of the major peculiarities noted above is that the output of the LPGS system is 'fast format'.
        This format (as produced by LPGS) cannot handle lats and lons properly and produces some ugly
        files that GDAL does not handle properly, and workarounds are required. The main workarounds
        are with the projection (search for "SPATIAL REFERENCE SPECIAL CASE" within this file). This
        sets some attributes on the instance that (must be) used somewhere, but I can't see where (see
        Roger about this if interested).

        Anyway, the upshot is that many of the methods called on this object will produce invalid results.
        (Roger also volunteered to fix this). Further even when this is fixed **ONE SHOULD NEVER TRY AND
        ACCESS THE UNDERLYING DATA USING GDAL AS IT WILL BE MISREPRESENTED**.
    """

    _TYPE_MAP = {'string': str, 'str': str, 'integer': int, 'int': int, 'float': float,
                 'boolean': bool, 'bool': bool, 'datetime': datetime, 'date': date, 'time': time}

    _CACHED_DATASETS = {}

    def __new__(cls, pathname=None, eAccess=gdalconst.GA_ReadOnly, default_metadata_required=True, utm_fix=False):
        """
        Implements a cache of SceneDatasets.

        This does not guarantee uniqueness. If a user calls :py:meth:`clear_cache`, then the next instance created
        with the same pathname after a delete, will point to a different instance.

        """
        if pathname and eAccess==gdalconst.GA_ReadOnly:
            return cls._CACHED_DATASETS.get(pathname) or super(SceneDataset, cls).__new__(cls, pathname, eAccess, default_metadata_required, utm_fix)
        return super(SceneDataset, cls).__new__(cls, pathname, eAccess, default_metadata_required, utm_fix)

    def __init__(self, pathname=None, eAccess=gdalconst.GA_ReadOnly, default_metadata_required=True, utm_fix=False):
        """Initialise an SceneDataset instance.

        :param pathname:
            Optional pathname of dataset to open.
        :type pathnname:
            str

        :param eAccess:
            Optional gdalconst value specifying file open mode.

        :param default_metadata_required:
            Boolean flag indicating whether all metadata should be required.

        :param utm_fix:
            Flag indicating whether UTM zone 60 should be fudged.
        :type utm_fix:
            bool
        """

        if eAccess==gdalconst.GA_ReadOnly and pathname and SceneDataset._CACHED_DATASETS.get(pathname):
            return

        self.default_metadata_required = default_metadata_required
        self.utm_fix = utm_fix

        def _initialise_from_xml():
            """
            Parse  XML file associated with this module.

            :todo:
                Only load the configuration the first time the class is instantiated.
            """
            logger.debug('  _initialise_from_xml() called')

            class_config_xml_file = re.sub('\.py(c){0,1}$', '.xml', __file__)
            self._class_config_dom_tree = xml.dom.minidom.parse(class_config_xml_file)
            class_nodes = self._class_config_dom_tree.documentElement.getElementsByTagName('CLASS')
            assert class_nodes, 'Missing CLASS element in ' + class_config_xml_file
            assert unicode_to_ascii(class_nodes[0].getAttribute('NAME')) == self.__class__.__name__, 'Invalid class name ' + unicode_to_ascii(class_nodes[0].getAttribute('NAME')) + ' in ' + class_config_xml_file
            class_node = class_nodes[0]
            data_nodes = class_node.getElementsByTagName('DATA')
            assert data_nodes, 'No DATA element defined under CLASS in ' + class_config_xml_file
            data_node = data_nodes[0]

            # Read in list of possible data directories (under root directory)
            datadir_nodes = data_node.getElementsByTagName('DATADIR')
            assert datadir_nodes, 'No DATADIR elements defined under DATA in ' + class_config_xml_file
            logger.debug('%s DATADIR nodes found', len(datadir_nodes))
            self._DATADIRS = []
            for datadir_node in datadir_nodes:
                self._DATADIRS.append(unicode_to_ascii(datadir_node.getAttribute('PATH')))
            logger.debug('self._DATADIRS = %s', self._DATADIRS)

            # Read in list of possible data directories (under root directory)
            filetype_nodes = data_node.getElementsByTagName('FILETYPE')
            assert filetype_nodes, 'No FILETYPE element defined under DATA in ' + class_config_xml_file
            self._FILE_TYPE_INFO = {}
            for filetype_node in filetype_nodes:
                extension = unicode_to_ascii(filetype_node.getAttribute('EXTENSION'))
                assert extension, 'No EXTENSION attribute defined under FILETYPE in ' + class_config_xml_file
                root_suffix = unicode_to_ascii(filetype_node.getAttribute('ROOTSUFFIX'))
                assert root_suffix, 'No ROOTSUFFIX attribute defined under FILETYPE in ' + class_config_xml_file
                product_format = unicode_to_ascii(filetype_node.getAttribute('PRODUCTFORMAT'))
                assert product_format, 'No PRODUCTFORMAT attribute defined under FILETYPE in ' + class_config_xml_file
                self._FILE_TYPE_INFO[extension] = {'ROOTSUFFIX': root_suffix, 'PRODUCTFORMAT': product_format}
            logger.debug('self._FILE_TYPE_INFO = %s', self._FILE_TYPE_INFO)

            # Read in metadata file matching patterns
            self._METADATA_SOURCES = {}
            metadata_nodes = class_node.getElementsByTagName('METADATA_SOURCE')
            assert metadata_nodes, 'No METADATA_SOURCE element defined under CLASS in ' + class_config_xml_file
            logger.debug('%s METADATA_SOURCE nodes found', len(metadata_nodes))
            for metadata_node in metadata_nodes:
                metadata_type = unicode_to_ascii(metadata_node.getAttribute('TYPE'))
                metadata_required = self.get_bool(unicode_to_ascii(metadata_node.getAttribute('REQUIRED'))) # Defaults to False
                self._METADATA_SOURCES[metadata_type] = {'REQUIRED': metadata_required}
                metadata_file_nodes = metadata_node.getElementsByTagName('FILE')
                assert metadata_file_nodes, 'No FILE elements defined under ' + metadata_type + ' in ' + class_config_xml_file
                logger.debug('%s FILE nodes found', len(metadata_file_nodes))
                file_patterns = []
                for metadata_file_node in metadata_file_nodes:
                    file_patterns.append(unicode_to_ascii(metadata_file_node.getAttribute('PATTERN')))
                assert file_patterns, 'No file matching pattern(s) defined for ' + metadata_type + ' in ' + class_config_xml_file
                self._METADATA_SOURCES[metadata_type]['FILE_PATTERNS'] = file_patterns

        try:
            logger.debug('SceneDataset() called')
            self._sub_datasets = [] # List of subdatasets (gdal.Dataset) managed
            self._file_list = [] # List of files managed
            self._raster_dict = {} # Dict containing lookup information for bands
            self._root_dataset_pathname = None # Pathname of root subdataset
            self._data_dir = None # Directory containing image files
            self._sub_dataset_type = None # 'TIF' or 'FST'
            self._xml_uses_attributes = False # Flag indicating whether XML values are stored as tag attributes
            self._metadata = Metadata() # Master metadata object managing all metadata
            self.satellite = None # Satellite object will be set from metadata values
            self._bands = {} # Dict containing lists of bands by type
            self._band_number_map = {} # Map file numbers to band numbers
            self.rgb_bands = [] # List of RGB band numbers in order
            self._metadata_vars = [] # List of variables read from metadata
            self.spatial_ref = None
            self.spatial_ref_geo = None
            self.cxform_to_geo = None
            self._root_band_number = None

            Dataset.__init__(self)

            _initialise_from_xml() # Read class-specific XML configuration

        except (AssertionError, DSException), e:
            logger.error( 'SceneDataset.__init__ error: %s', e.message)
            raise DSException(e.message)

        if pathname:
            self.Open(pathname, eAccess)

    def __del__(self):
        # Use the method clear_cache if you really want to remove the dataset from the cache.
        super(SceneDataset, self).__del__()

    @staticmethod
    def clear_cache(pathname_or_instance=None):
        """
        Drop either a specified dataset or all datasets from the cache.

        :param pathname_or_instance:
            If ``None``, then drop all datasets, otherwise dump the dataset specified by name or instance.
        :type pathname_or_instance:
            Either a :py:class:SceneDataset`, or a pathname for one.
        """
        if not pathname_or_instance:
            SceneDataset._CACHED_DATASETS.clear()
        elif type(pathname_or_instance) == str:
            SceneDataset._CACHED_DATASETS.pop(pathname_or_instance, None)
        elif type(pathname_or_instance) == Dataset:
            SceneDataset._CACHED_DATASETS.pop(pathname_or_instance.pathname, None)

    def Open(self, pathname, eAccess=gdalconst.GA_ReadOnly):
        """
        Non-GDAL function to open a scene directory
        Overrides Dataset.Open function and emulates it as closely as possible

        :param pathname:
            The root pathname of scene (e.g. parent of scene01 if present).
        :type pathname:
            str

        :param eAccess:
            Mode in which to open the dataset (:py:const:`gdalconst.GA_ReadOnly` or :py:const:`gdalconst.GA_Update`).

        :return:
            ``self``
        """

        def _gather_metadata():
            """
            Function to read ALL available metadata into a single nested dict tree structure
            """
            logger.debug('  _gather_metadata() called')
            for metadata_type in ['XML', 'MTL', 'REPORT']:
                for file_pattern in self._METADATA_SOURCES[metadata_type]['FILE_PATTERNS']:
                    file_list = find_files(self._pathname, file_pattern)
                    if file_list:
                        logger.debug('Reading %s metadata file %s', metadata_type, file_list[0])
                        if metadata_type == 'XML':
                            xml_metadata = XMLMetadata(file_list[0])
                            self._xml_uses_attributes = xml_metadata.uses_attributes # Remember whether attributes are used in XML
                            self._metadata.merge_root_metadata_from_object(metadata_object=xml_metadata, overwrite=False)
                        elif metadata_type == 'MTL':
                            self._metadata.merge_root_metadata_from_object(metadata_object=MTLMetadata(file_list[0]))
                        elif metadata_type == 'REPORT':
                            self._metadata.set_root_metadata_from_object(metadata_object=ReportMetadata(file_list[0]))
                    else:
                        logger.debug('No %s files matching "%s" found', metadata_type, file_pattern)

                assert self._metadata.get_metadata([metadata_type]) or not self._METADATA_SOURCES[metadata_type]['REQUIRED'], (
                    'No ' + metadata_type + ' metadata file(s) found for ' + self._pathname)

            #===================================================================
            # # Import image metadata (if any) and store under 'FST' or 'TIF as required
            # image_metadata = self._root_dataset.GetMetadata_Dict()
            # if image_metadata:
            #    self._metadata.set_root_metadata(self._sub_dataset_type, image_metadata)
            #===================================================================

            # Set all instance attributes as defined in class XML file
            self.read_metadata()

            self.satellite = Satellite(self.satellite_name, self.sensor)
            assert self.satellite, 'Unable to create Satellite object for %s %s' % (self.satellite_name, self.sensor)

            self._bands = dict(self.satellite.BAND_TYPES) # Copy dict from satellite
            for band_type in self._bands:
                self._bands[band_type] = []

        def _set_instance_values():
            """
            Function to set any class-specific attributes derived from metadata
            This is called after read_metadata(). Any derived vales should be set here.
            """
            def get_mtl_bias_gain():
                '''
                Calculate bias & gain for each band from MTL metadata
                Sets self.bias & self.gain lists with values for each band
                Code adapted from old ula.metadata.py
                N.B: Will exit with warning if required metadata not found
                '''
                bdict = {}
                gdict = {}

                # Get MIN_MAX_RADIANCE dict from MTL data
                mmr = self._metadata.get_metadata('MTL,L1_METADATA_FILE,MIN_MAX_RADIANCE')
                if mmr:
                    # Get MIN_MAX_RADIANCE dict from MTL data
                    mmpv = self._metadata.get_metadata('MTL,L1_METADATA_FILE,MIN_MAX_PIXEL_VALUE')
                    if mmpv:
                        # TODO: This will probably break with Landsat 8 data
                        for key in sorted(mmr.keys()):
                            m = re.match('LMIN_BAND(\d+)', key)
                            if m:
                                mmr_number = int(m.group(1))
                                if 60 <= mmr_number < 70:
                                    band_file_number = mmr_number
                                else:
                                    band_file_number = mmr_number * 10
#                                logger.info('key = %s, band_file_number = %s', key, band_file_number)
                                assert band_file_number in self.satellite.BAND_TYPES['ALL'], 'Invalid band file number %d' % band_file_number
                                if band_file_number in self._band_number_map:
                                    band_number = self._band_number_map[band_file_number]

#                                    bias = gain = 0.0
#                                    if band_number in self._bands['REFLECTIVE']:
                                    lmin = float(mmr['LMIN_BAND%d' % mmr_number])
                                    lmax = float(mmr['LMAX_BAND%d' % mmr_number])
                                    qcalmin = float(mmpv['QCALMIN_BAND%d' % mmr_number])
                                    qcalmax = float(mmpv['QCALMAX_BAND%d' % mmr_number])
                                    gain = (lmax - lmin) / (qcalmax - qcalmin)
                                    bias = lmax - gain * qcalmax

#                                    logger.info('mmr_number = %s, band_number = %s, lmin = %s, lmax = %s, gain = %s, bias = %s',
#                                                mmr_number, band_number, lmin, lmax, gain, bias)

                                    bdict[band_number] = bias
                                    gdict[band_number] = gain
                    else:
                        logger.debug('No MIN_MAX_PIXEL_VALUE data found in MTL file. Unable to compute bias & gain.')
                else:
                    logger.debug('No MIN_MAX_RADIANCE data found in MTL file. Unable to compute bias & gain.')

                self.bias = bdict
                self.gain = gdict

            def get_mtl_bias_gain_landsat8_lookup():

                bdict = {}
                gdict = {}

                rad_rescale_params = self._metadata.get_metadata('MTL,L1_METADATA_FILE,RADIOMETRIC_RESCALING')

                print
                print rad_rescale_params
                pprint(rad_rescale_params)

                # Need to check if we have a valid return from the metadata lookup
                if rad_rescale_params:
                    for key, value in rad_rescale_params.iteritems():

                        # Bias: RADIANCE_ADD_BAND_X = <value>

                        match_add = re.match('RADIANCE_ADD_BAND_(\d+)', key)
                        if match_add:
                            bdict[ int(match_add.group(1)) ] = float(value)

                        # Gain: RADIANCE_MULT_BAND_X = <value>

                        match_mult = re.match('RADIANCE_MULT_BAND_(\d+)', key)
                        if match_mult:
                            gdict[ int(match_mult.group(1)) ] = float(value)
                else:
                    logger.debug('No RADIOMETRIC_RESCALING data found in MTL file. Unable to compute bias & gain.')

                print
                print 'bdict'
                pprint(bdict)
                print
                print 'gdict'
                pprint(gdict)
                print

                self.bias = bdict
                self.gain = gdict



            if self.scene_centre_date and self.scene_centre_time:
                self.scene_centre_datetime = datetime(
                    self.scene_centre_date.year,
                    self.scene_centre_date.month,
                    self.scene_centre_date.day,
                    self.scene_centre_time.hour,
                    self.scene_centre_time.minute,
                    self.scene_centre_time.second,
                    self.scene_centre_time.microsecond)

                self.scene_start_datetime = self.scene_centre_datetime - self.satellite.acquistion_seconds // 2
                self.scene_end_datetime = self.scene_start_datetime + self.satellite.acquistion_seconds
            else:
                self.scene_centre_datetime = None
                self.scene_start_datetime = None
                self.scene_end_datetime = None

            if self.completion_date and self.completion_time:
                self.completion_datetime = datetime(
                    self.completion_date.year,
                    self.completion_date.month,
                    self.completion_date.day,
                    self.completion_time.hour,
                    self.completion_time.minute,
                    self.completion_time.second,
                    self.completion_time.microsecond)
            else:
                self.completion_datetime = None

            # SPATIAL REFERENCE GENERAL CASE
            # Use the projection and geotransform that GDAL gives us.
            # WARNING: in special case code below, always create a new spatial reference
            # instance. Modifying the instance created here DOES NOT WORK.

            self.spatial_ref = osr.SpatialReference()
            self.spatial_ref.ImportFromWkt(self.GetProjection())

            # Initialise geotransform to None (otherwise it calls the SceneDataset Class function
            # GetGeoTransform() which is supposed to override the Dataset Class function of the same name.
            self.geotransform = None
            self.geotransform = self.GetGeoTransform()

            # SPATIAL REFERENCE SPECIAL CASE
            # FAST-EQR dataset. HRF.FST files coming from EODS may have PROJECTION=EQR,
            # but also contain UTM-like extents and pixel size, which are in the spec
            # but give GDAL fits -- GetProjection() and GetGeoTransform() return
            # unusable values. To fix this situation, create new spatial reference and
            # geotransform objects using values derived from metadata. WARNING: Client
            # code should always use <dataset>.geotransform rather than
            # <dataset>.GetGeoTransform().

            if self.product_format.startswith('FAST') and self.map_projection.startswith('EQR'):
                self.spatial_ref = osr.SpatialReference()
                self.spatial_ref.SetWellKnownGeogCS('WGS84')

                #self.geotransform = [
                #    self.ul_lon - self.pixel_x_size / 2.,
#               #     (self.ur_lon - self.ul_lon) / self.image_pixels,
                #    self.pixel_x_size,
                #    0.0,
                #    self.ul_lat + self.pixel_y_size / 2.,
                #    0.0,
#               #     (self.ll_lat - self.ul_lat) / self.image_lines,
                #    -self.pixel_y_size
                #]

                gt = self.GetGeoTransform()
                self.geotransform = (numpy.array(gt)/100000.).tolist()
                # Need to reset the rotation params. These are probably zero anyway as our imagery 
                # is created North-Up (for L1T), but it is better to be safe. But then again rotation
                # values might be scaled incorrectly???
                self.geotransform[2] = gt[2]
                self.geotransform[4] = gt[4]
                self.coordinate_reference_system = 'WGS84'

                # Force calculation of the corner points - metadata values may be scaled incorrectly
                self.ul_x = None
                self.ul_y = None
                self.ur_x = None
                self.ur_y = None
                self.ll_x = None
                self.ll_y = None
                self.lr_x = None
                self.lr_y = None

                logger.warning('FAST-EQR format, created substitute spatial reference and geotransform')

            # SPATIAL REFERENCE SPECIAL CASE
            # UTM zone outside of GDA94 datum specification, e.g. UTM zone 60 for GDA94.
            # Substitute a WGS84 spatial reference. Differences between GDA94 and WGS84
            # should be sub-pixel, see http://www.icsm.gov.au/gda/wgs84fact.pdf.
            # *************************************************************************
            # *** This workaround is likely wrong, based on further info from Matthew.
            # *** Need ITRF2008 support.
            # *************************************************************************

            _local = self.spatial_ref.GetAttrValue('LOCAL_CS')
            if _local:
                assert self.utm_fix, 'LOCAL_CS [%s] (UTM zone could be outside the datum spec) -- ABORTED' % _local

                m = re.search('UTM Zone (\d+) ', _local)  # hardwired
                if m:
                    zone = int(m.group(1))
                else:
                    assert False, 'Could not parse LOCAL_CS\n%s' % _local
                self.spatial_ref = osr.SpatialReference()
                self.spatial_ref.SetWellKnownGeogCS('WGS84')
                self.spatial_ref.SetUTM(zone, False)
                self.geotransform = [
                    self.ul_lon,
                    (self.ur_lon - self.ul_lon) / self.image_pixels,
                    0.0,
                    self.ul_lat,
                    0.0,
                    (self.ll_lat - self.ul_lat) / self.image_lines,
                ]
                logger.warning('Overrode LOCAL_CS spatial reference with WGS84 (zone=%d)' % zone)

            #
            # TODO ...catch other spatialreference/projection error cases...
            #

            self.spatial_ref_geo = self.spatial_ref.CloneGeogCS()

            self.cxform_to_geo = osr.CoordinateTransformation(self.spatial_ref, self.spatial_ref_geo)
            self.cxform_from_geo = osr.CoordinateTransformation(self.spatial_ref_geo, self.spatial_ref)

            extents        = self.GetExtent()
            array_extents  = numpy.array(extents)
            centre_x       = float(numpy.mean(array_extents[:,0]))
            centre_y       = float(numpy.mean(array_extents[:,1]))
            extents.append([centre_x,centre_y])

            #self.lonlats = {
            #    'CENTRE': (self.scene_centre_long, self.scene_centre_lat),
            #    'UL': (self.ul_lon, self.ul_lat),
            #    'UR': (self.ur_lon, self.ur_lat),
            #    'LL': (self.ll_lon, self.ll_lat),
            #    'LR': (self.lr_lon, self.lr_lat)}

            if self.IsGeographic():
                print 'IsGeographic'
                self.lonlats = {
                    'CENTRE' : (extents[4][0], extents[4][1]),
                    'UL'     : (extents[0][0], extents[0][1]),
                    'UR'     : (extents[2][0], extents[2][1]),
                    'LL'     : (extents[1][0], extents[1][1]),
                    'LR'     : (extents[3][0], extents[3][1])
                               }
                print self.lonlats
                # If the scene is natively in geographics, we shouldn't need to 
                # project the co-ordinates to UTM.

                # Set the georeferenced coordinates of the corner points if we don't already have them.
                # These generally only get set when the product is FAST-EQR when they're forced to None
                if not (self.ul_x and self.ul_y):
                    self.ul_x, self.ul_y = self.lonlats['UL']
                if not (self.ur_x and self.ur_y):
                    self.ur_x, self.ur_y = self.lonlats['UR']
                if not (self.ll_x and self.ll_y):
                    self.ll_x, self.ll_y = self.lonlats['LL']
                if not (self.lr_x and self.lr_y):
                    self.lr_x, self.lr_y = self.lonlats['LR']


                self.scene_centre_x, self.scene_centre_y = self.lonlats['CENTRE']
            else:
                self.coords = {
                    'CENTRE' : (extents[4][0], extents[4][1]),
                    'UL'     : (extents[0][0], extents[0][1]),
                    'UR'     : (extents[2][0], extents[2][1]),
                    'LL'     : (extents[1][0], extents[1][1]),
                    'LR'     : (extents[3][0], extents[3][1])
                              }

                re_prj_extents=[]
                for x,y in extents:
                    new_x, new_y, new_z = self.cxform_to_geo.TransformPoint(x,y)
                    re_prj_extents.append([new_x,new_y])
                    print new_x, new_y

                self.lonlats = {
                    'CENTRE' : (re_prj_extents[4][0], re_prj_extents[4][1]),
                    'UL'     : (re_prj_extents[0][0], re_prj_extents[0][1]),
                    'UR'     : (re_prj_extents[2][0], re_prj_extents[2][1]),
                    'LL'     : (re_prj_extents[1][0], re_prj_extents[1][1]),
                    'LR'     : (re_prj_extents[3][0], re_prj_extents[3][1])
                               }

                # Set the georeferenced coordinates of the corner points if we don't already have them.
                # These generally only get set when the product is FAST-EQR when they're forced to None
                if not (self.ul_x and self.ul_y):
                    self.ul_x, self.ul_y = self.coords['UL']
                if not (self.ur_x and self.ur_y):
                    self.ur_x, self.ur_y = self.coords['UR']
                if not (self.ll_x and self.ll_y):
                    self.ll_x, self.ll_y = self.coords['LL']
                if not (self.lr_x and self.lr_y):
                    self.lr_x, self.lr_y = self.coords['LR']


                self.scene_centre_x, self.scene_centre_y = self.coords['CENTRE']


            #self.scene_centre_x, self.scene_centre_y = self.coord_from_geo(self.lonlats['CENTRE'])
            #self.scene_centre_x, self.scene_centre_y = self.coords['CENTRE']

            # Set the georeferenced coordinates of the corner points if we don't already have them
            #if not (self.ul_x and self.ul_y):
            #    #self.ul_x, self.ul_y = self.coord_from_geo(self.lonlats['UL'])
            #    self.ul_x, self.ul_y = self.coords['UL']
            #if not (self.ur_x and self.ur_y):
            #    #self.ur_x, self.ur_y = self.coord_from_geo(self.lonlats['UR'])
            #    self.ur_x, self.ur_y = self.coords['UR']
            #if not (self.ll_x and self.ll_y):
            #    #self.ll_x, self.ll_y = self.coord_from_geo(self.lonlats['LL'])
            #    self.ll_x, self.ll_y = self.coords['LL']
            #if not (self.lr_x and self.lr_y):
            #    #self.lr_x, self.lr_y = self.coord_from_geo(self.lonlats['LR'])
            #    self.lr_x, self.lr_y = self.coords['LR']

            #self.coords = {
            #    'CENTRE': (self.scene_centre_x, self.scene_centre_y),
            #    'UL': (self.ul_x, self.ul_y),
            #    'UR': (self.ur_x, self.ur_y),
            #    'LL': (self.ll_x, self.ll_y),
            #    'LR': (self.lr_x, self.lr_y)}

            # Pre-compute some useful derived time quantities
            if self.scene_centre_date:
                self.DOY = self.scene_centre_date.timetuple().tm_yday

            if self.scene_centre_time:
                self.decimal_hour = (self.scene_centre_time.hour +
                    (self.scene_centre_time.minute +
                     (self.scene_centre_time.second +
                      self.scene_centre_time.microsecond / 1000000.0) / 60.0) / 60.0)

            if self.scene_centre_date and self.scene_centre_time:
                self.decimal_day = self.DOY + self.decimal_hour / 24.0

            MTL_bias_gain_dict = {
                          'Landsat-5' : get_mtl_bias_gain,
                          'Landsat-7' : get_mtl_bias_gain,
                          'Landsat-8' : get_mtl_bias_gain_landsat8_lookup,
                                 }
            MTL_bias_gain_dict[self.satellite.NAME]()
            #get_mtl_bias_gain()
            #get_mtl_bias_gain_landsat8_lookup()

            # Deal with special case where datum/ellipsoid requires +ve UTM zone
            if self.zone and self.zone < 0 and self.datum == 'GDA94' and self.earth_ellipsoid == 'GRS80':
                self.zone = abs(self.zone)
                self.update_metadata('zone')


        logger.debug('SceneDataset.Open(%s,%s) called', repr(pathname), repr(eAccess))
        self._pathname = os.path.abspath(pathname)
        self._eAccess = eAccess

        try:
            # Apply most basic validity checks here
            assert os.path.isdir(pathname), pathname + ' is not a directory'

            # Read all metadata into nested dict and look up satellite info
            _gather_metadata()

            # Find root dataset
            thermal_dataset = None
            panchromatic_dataset = None
            self._data_dir = None
            for datadir in self._DATADIRS:
                datadir = os.path.abspath(os.path.join(self._pathname, datadir))
                if os.path.isdir(datadir):
                    for f in sorted(os.listdir(datadir)):
                        m = re.search('(\w+)_(\w+).(\w+)$', f)
                        if m:
                            suffix = m.group(2).upper()
                            extension = m.group(3).upper()

                            print 'suffix         ', suffix
                            print 'extension      ', extension
                            print '_FILE_TYPE_INFO', self._FILE_TYPE_INFO
                            print 'extension in _FILE_TYPE_INFO', (extension in self._FILE_TYPE_INFO)
                            print 'satellite.root_band', self.satellite.root_band

                            # Determine the root suffix for the found file type
                            if extension in self._FILE_TYPE_INFO:
                                try:
                                    root_suffix = self._FILE_TYPE_INFO[extension]['ROOTSUFFIX'] % self.satellite.root_band
                                except (TypeError):
                                    root_suffix = self._FILE_TYPE_INFO[extension]['ROOTSUFFIX']

                                print 'root_suffix', root_suffix, (suffix == root_suffix)

                                if suffix == root_suffix:
                                    self._root_dataset_pathname = os.path.abspath(os.path.join(datadir, f))
                                    self._sub_dataset_type = extension
                                    self._data_dir = datadir

                                    print '_root_dataset_pathname', self._root_dataset_pathname

                                    try: # ToDo: Need to deal with open failure here
                                        self._root_dataset = gdal.Open(self._root_dataset_pathname, eAccess)
                                    except (RuntimeError), e: # ToDo: Need to handle specific exception type
                                        self._root_dataset = None
                                        raise DSException(e.message)

                                    self._sub_datasets.append(self._root_dataset)

                                elif extension == 'FST' and suffix == 'HTM': # Found Thermal FST dataset
                                    thermal_dataset = gdal.Open(os.path.join(datadir, f), eAccess)

                                elif extension == 'FST' and suffix == 'HPN': # Found Pan-Chromatic FST dataset
                                    panchromatic_dataset = gdal.Open(os.path.join(datadir, f), eAccess)

                if self._data_dir and self._root_dataset_pathname: # Data directory has been determined
                    break # Stop searching

            print
            print '_data_dir             ', self._data_dir
            print '_root_dataset_pathname', self._root_dataset_pathname

            assert self._data_dir and self._root_dataset_pathname, 'Unable to find root dataset under ' + self._pathname

            # Perform basic checks on root dataset
            assert self._sub_dataset_type, 'Unable to determine dataset type'
            assert self._root_dataset, 'Unable to open root dataset'

            # check all files and set up references to subdatasets. Bands are loaded in numerical order
            for band_index in range(len(self.satellite.BAND_LIST)):
                band_number = band_index + 1
                sensor_band_info = self.satellite.BAND_LIST[band_index]
                band_file_number = sensor_band_info['NUMBER']
                file_list = glob(os.path.join(self._data_dir, '*_B' + str(band_file_number) + '.' + self._sub_dataset_type.lower()))
                file_list += glob(os.path.join(self._data_dir, '*_B' + str(band_file_number) + '.' + self._sub_dataset_type.upper()))
                if len(file_list) == 1:# There can be only one...
                    # Use absolute pathnames to ensure string matches
                    band_file = os.path.abspath(file_list[0])

                    band_dict = {
                        'file': band_file,
                        'band_file_number': band_file_number,
                        'type': sensor_band_info['TYPE']}

                    # TIF only ever one raster per subDataset
                    if self._sub_dataset_type == 'TIF':
                        # Do not duplicate _root_dataset - already opened
                        if band_dict['file'] == self._root_dataset_pathname:
                            band_dict['dataSet'] = self._root_dataset
                        else:
                            band_dict['dataSet'] = gdal.Open(band_dict['file'], eAccess)
                            # ToDo: Perform sanity checks on sub-datasets here
                            # (e.g. same extents)
                            self._sub_datasets.append(band_dict['dataSet'])

                        band_dict['rasterIndex'] = 1

                    # FST has multiple rasters per subDataset
                    elif self._sub_dataset_type == 'FST':

                        # Multiple reflective bands in one dataset
                        if band_file_number in self.satellite.BAND_TYPES['REFLECTIVE'] and self._root_dataset:

                            # Filter *.aux.xml files from the root dataset file list.
                            # These XML files mess up band_dict construction.
                            __flist = [x for x in self._root_dataset.GetFileList() if not x.endswith('.aux.xml')]
                            band_dict['rasterIndex'] = __flist.index(file_list[0])
                            band_dict['dataSet'] = self._root_dataset

                        # Two thermal bands in one dataset
                        elif band_file_number in self.satellite.BAND_TYPES['THERMAL'] and thermal_dataset:
                            band_dict['dataSet'] = thermal_dataset
                            band_dict['rasterIndex'] = thermal_dataset.GetFileList().index(file_list[0])
                            if thermal_dataset not in self._sub_datasets:
                                self._sub_datasets.append(thermal_dataset)

                        # Pan-chromatic band is in its own dataset
                        elif band_file_number in self.satellite.BAND_TYPES['PANCHROMATIC'] and panchromatic_dataset:
                            band_dict['dataSet'] = panchromatic_dataset
                            # Only ever one band in pan-chromatic dataset
                            band_dict['rasterIndex'] = 1
                            if panchromatic_dataset not in self._sub_datasets:
                                self._sub_datasets.append(panchromatic_dataset)

                        else:
                            raise Exception('Invalid band file number ' + str(band_file_number))

                    # Register lookup information for band
                    self._raster_dict[band_number] = band_dict
                    self._band_number_map[band_file_number] = band_number
                    self._bands['ALL'].append(band_number)
                    self._bands[sensor_band_info['TYPE']].append(band_number)
                    if band_file_number == self.satellite.root_band:
                        self._root_band_number = band_number
                else:
                    logger.debug('%s files found for band file number %s' % (len(file_list), band_file_number))

            # Validity checks for required bands
            assert len(self._raster_dict), 'No band information found'
            if len(self._raster_dict) < len(self.satellite.BAND_TYPES['ALL']):
                logger.debug('Only %s/%s bands found' % (len(self._raster_dict), len(self.satellite.BAND_TYPES['ALL'])))

            # Remove duplicates and sort file list
            fileSet = set()
            for subDataset in self._sub_datasets:
                fileSet |= set(subDataset.GetFileList())
            self._file_list = sorted(list(fileSet))

            if self.satellite.rgb_bands:
                self.rgb_bands = [self._band_number_map[band_file_number] for band_file_number in self.satellite.rgb_bands]

            # Set any class-specific attributes derived from metadata
            _set_instance_values()

            # Dataset opened successfully
            return self

        except (AssertionError, DSException), e:
            logger.error( 'SceneDataset.Open(' + pathname + ') error: %s', e.message)
            raise DSException(e.message)

        if eAccess==gdalconst.GA_ReadOnly:
            SceneDataset._CACHED_DATASETS[pathname] = self


    def read_metadata(self, metadata_type_id = None, class_config_dom_tree = None, ignore_empty_values=False):
        """
        Parse an XML specification file and set instance values from metadata performing type conversions from strings.

        :param metadata_type_id:
            Metadata type (e.g. 'XML', 'MTL', 'REPORT', 'TIF', 'FST') or None for all types.
        :type metadata_type_id:
            str

        :param class_config_dom_tree:
            DOM tree containing XML specification file defining the variables to be read from the metadata
            This object contains the name, type, metadata path(s) of a range of variables to be read and returned.
        :type class_config_dom_tree:
            ???

        :param ignore_empty_values:
            Flag indicating whether empty strings should not overwrite existing values.
        :type ignore_empty_values:
            bool
        """
        logger.debug('read_metadata(%s, %s, %s) called',
                     repr(class_config_dom_tree), repr(metadata_type_id), repr(ignore_empty_values))

        class_config_dom_tree = class_config_dom_tree or self._class_config_dom_tree
        var_list = []

        try:
            variable_nodes = class_config_dom_tree.getElementsByTagName('VARIABLE')
            for variable_node in variable_nodes:
                # Variable name should always be defined
                var_name = unicode_to_ascii(variable_node.getAttribute('NAME'))
                assert var_name, 'Variable has no name defined'
                logger.debug('Variable name = %s', var_name)

                # Optional TYPE attribute - defaults to str
                var_type = self.__class__._TYPE_MAP.get(variable_node.getAttribute('TYPE'), str)
                logger.debug('  Variable type = %s', repr(var_type))

                # Optional REQUIRED attribute - defaults to True
                var_required = self.get_bool(unicode_to_ascii(variable_node.getAttribute('REQUIRED')),
                                             self.default_metadata_required)

                logger.debug('  Value required = %s', var_required)

                var_metatada_nodes = variable_node.getElementsByTagName('METADATA')
                assert var_metatada_nodes, 'No metadata search paths defined for ' + var_name
                for var_metatada_node in var_metatada_nodes:
                    var_key = unicode_to_ascii(var_metatada_node.getAttribute('KEY'))
                    assert var_key, 'Metadata key not defined'
                    logger.debug('  Searching metadata under = %s', repr(var_key))

                    # Skip any non-matching searches if metadata_type_id is specified
                    if metadata_type_id and not re.match('^' + metadata_type_id +',', var_key):
                        continue

                    var_string = self._metadata.get_metadata(var_key.split(',')) # Search for metadata

                    if type(var_string) == dict: # Incomplete metadata path - not a leaf node
                        var_string = None

                    if var_string is not None: # Value not found in metadata (empty string allowed)
                        # Check for optional regex and even more optional group
                        var_regex = unicode_to_ascii(var_metatada_node.getAttribute('REGEX')) or None
                        var_regex_group = 0
                        if var_regex:
                            try:
                                var_regex_group = int(unicode_to_ascii(var_metatada_node.getAttribute('GROUP')))
                            except ValueError:
                                var_regex_group = 0

                            logger.debug('  Regular expression = %s, group = %s', var_regex, var_regex_group)
                            s = re.search(var_regex, var_string)
                            if s:
                                var_string = s.group(var_regex_group)
                            else:
                                var_string = None

                            logger.debug('  RegEx processed metadata string = %s', repr(var_string))

                    if var_string is not None: # Value not found in metadata (empty string allowed)
                        break # Stop searching metadata

                assert var_string or ignore_empty_values or not var_required, var_name + ' not found in metadata'
                var_string = var_string or ''

                logger.debug('  Raw metadata string = %s', repr(var_string))

                var_value = None
                if not var_string:
                    # Return empty string for missing string value, None for every other type
                    if var_type == str:
                        var_value = ''
                    else:
                        var_value = None
                elif var_type == int:
                    var_value = int(var_string)
                elif var_type == float:
                    var_value = float(var_string)
                elif var_type == bool:
                    var_value = self.get_bool(var_string)
                elif var_type == datetime or var_type == date or var_type == time:
                    # Try parsing datetime using format string
                    var_format = unicode_to_ascii(var_metatada_node.getAttribute('FORMAT'))
                    assert var_format, 'Format string must be provided for datetime/date/time metadata'
                    try:
                        logger.debug('  Trying datetime format string = %s', repr(var_format))
                        if var_type == datetime:
                            var_value = datetime.strptime(var_string, var_format)
                        elif var_type == date:
                            var_value = datetime.strptime(var_string, var_format).date()
                        elif var_type == time:
                            var_value = datetime.strptime(var_string, var_format).time()
                    except (ValueError), e:
                        logger.debug('  Datetime parsing failed: %s', e.message)

                    assert var_value or ignore_empty_values or not var_required, 'Invalid datetime format'
                else: #str type
                    var_value = var_string

                logger.debug('  Metadata value %s = %s', var_name, repr(var_value))

                if var_value or not ignore_empty_values:
                    var_list.append((var_name, var_value))
                else:
                    logger.debug('  Empty metadata for %s ignored', var_name)

            # All or nothing - only update instance values if all metadata is OK
            logger.info('%s values read from metadata', len(var_list))
            for var_name, var_value in var_list:
                self._metadata_vars.append(var_name)
                self.__setattr__(var_name, var_value)

        except (AssertionError, DSException), e:
            logger.error( 'SceneDataset.read_metadata error: %s', e.message)
            raise DSException(e.message)


    def update_metadata(self,
        var_name_or_var_name_list = None,
        first_match_only = False,
        create_new_nodes = False,
        metadata_type_id = None,
        class_config_dom_tree = None,
        keep_existing_values = False):
        """Function to parse an XML specification file and set values in the metadata_dict from the list of variable
        names specified.
        Argument:
            metadata_type_id: metadata type string (e.g. 'XML', 'MTL', 'REPORT', 'TIF', 'FST') or None for all types
            class_config_dom_tree: DOM tree containing XML specification file defining the variables to be read from the metadata
            This object contains the name, type, metadata path(s) of a range of variables to be read and returned.
            keep_existing_values: Boolean flag indicating whether existing non-empty values in metadata dict should be overwritten

        :todo:
            The documentation for this method is wrong and needs to be fixed.
        """
        logger.debug('update_metadata(%s, %s, %s, %s, %s, %s) called',
                     repr(var_name_or_var_name_list), repr(first_match_only),
                     repr(create_new_nodes), repr(metadata_type_id), repr(class_config_dom_tree),
                     repr(keep_existing_values))

        class_config_dom_tree = class_config_dom_tree or self._class_config_dom_tree

        # Allow single variable to be set
        if type(var_name_or_var_name_list) == str:
            var_list = [var_name_or_var_name_list]
        else:
            var_list = var_name_or_var_name_list

        logger.debug('var_list = %s', var_list)

        vars_found = set()
        vars_changed = {} # Dict keyed by var_name containing lists of keys changed

        try:
            variable_nodes = class_config_dom_tree.getElementsByTagName('VARIABLE')
            for variable_node in variable_nodes:
                # Variable name should always be defined
                var_name = unicode_to_ascii(variable_node.getAttribute('NAME'))
                assert var_name, 'Variable has no name defined'

                # Skip variable_node if not included in list to write
                if var_list and not var_name in var_list:
                    continue

                logger.debug('Variable name = %s', var_name)

                # Optional TYPE attribute - defaults to str
                var_type = self.__class__._TYPE_MAP.get(variable_node.getAttribute('TYPE'), str)
                logger.debug('  Variable type = %s', repr(var_type))

                var_metatada_nodes = variable_node.getElementsByTagName('METADATA')
                assert var_metatada_nodes, 'No metadata search paths defined for ' + var_name

                # Only look at first (preferred) metadata definition if first_match_only is True
                if first_match_only:
                    var_metatada_nodes = var_metatada_nodes[0:1]

                # Set all required metadata values in tree
                for var_metatada_node in var_metatada_nodes:
                    var_key = unicode_to_ascii(var_metatada_node.getAttribute('KEY'))
                    assert var_key, 'Metadata key not defined'

                    # Skip any non-matching searches if metadata_type_id is specified
                    if metadata_type_id and not re.match('^' + metadata_type_id +',', var_key):
                        logger.debug('Skipping non-%s %s', metadata_type_id, var_key)
                        continue

                    if vars_changed.get(var_name) and var_key in vars_changed[var_name]:
                        logger.debug('%s already changed for %s', var_key, var_name)
                        continue

                    logger.debug('  Searching metadata under = %s', repr(var_key))

                    old_var_string = self._metadata.get_metadata(var_key.split(','))

                    assert type(old_var_string) != dict, 'Incomplete metadata path ' + var_key + ' (not a leaf node)' # Could be None

                    if not create_new_nodes and old_var_string is None:
                        logger.debug('%s not found in metadata', var_key)
                        continue

                    logger.debug('  Old metadata string = %s', repr(old_var_string))

                    # Check for optional regex and even more optional group
                    var_regex = unicode_to_ascii(var_metatada_node.getAttribute('REGEX')) or None
                    var_regex_group = 0
                    if var_regex and old_var_string:
#                        assert old_var_string, 'Unable to apply regular expression to non-existent metadata text'
                        try:
                            var_regex_group = int(unicode_to_ascii(var_metatada_node.getAttribute('GROUP')))
                        except ValueError:
                            var_regex_group = 0

                    var_value = self.__getattribute__(var_name)
                    assert var_type == type(var_value) or var_value is None, 'Variable type mismatch for ' + var_name + '. (Expected ' + var_type.__name__ + ', found ' + type(var_value).__name__ + ')'
                    logger.debug('  Metadata value %s = %s', var_name, repr(var_value))

                    # Get format string (Empty string if no format string specified)
                    var_format = unicode_to_ascii(var_metatada_node.getAttribute('FORMAT'))
                    logger.debug('  var_format = %s', repr(var_format))

                    if var_value is None:
                        var_string = '' # Convert None to empty string
                    elif var_type == datetime or var_type == date or var_type == time:
                        assert var_format, 'Format string must be provided for datetime/date/time metadata'

                        # Format datetime using format string
                        logger.debug('  Formatting date/time %s using format string %s', repr(var_value), repr(var_format))
                        var_string = var_value.strftime(var_format)
                    elif var_format: # Format string provided
                        var_string = var_format % var_value
                    else: # Convert to str type
                        var_string = format(var_value)

                    logger.debug('  Metadata string %s = %s', var_name, repr(var_string))

                    # If RegEx exists - do the substitution
                    if var_regex and old_var_string:
                        logger.debug('  Regular expression = %s, group = %s', var_regex, var_regex_group)

                        s = re.search(var_regex, old_var_string)
                        if not s: # Don't change string if we can't find a match in the old one
                            logger.debug('Regex "%s" not found in old string "%s"', var_regex, old_var_string)
                            continue

                        try:
                            var_string = old_var_string[0:s.start(var_regex_group)] + var_string + old_var_string[s.end(var_regex_group): ]
                            logger.debug('  Substituted metadata string %s = %s', var_name, repr(var_string))
                        except (IndexError):
                            logger.debug('Regex group "%s" not found in search result "%s"', var_regex_group, s.groups())
                            continue

                    vars_found.add(var_name)
                    if var_string != old_var_string and not (keep_existing_values and old_var_string):
                        self._metadata.set_metadata_node(var_key.split(','), var_string, overwrite = True)
                        logger.debug('  Metadata %s set to %s', var_key, var_string)

                        if not vars_changed.get(var_name): # No var_key list defined for var_name
                            vars_changed[var_name] = [var_key] # Create new var_key list
                        else:
                            vars_changed[var_name].append(var_key) # Append var_key to existing list

            logger.debug('Specified variables found in metadata : %s', sorted(list(vars_found)))
            logger.debug('Specified variables changed in metadata : %s', sorted(vars_changed.keys()))

            if var_list:
                vars_not_found = set(var_list) - vars_found
            else:
                vars_not_found = set()

            assert not vars_not_found, str(len(vars_not_found)) + ' specified variables not found in metadata : ' + str(sorted(list(vars_not_found)))

            logger.info('%s / %s values changed in metadata', len(vars_changed), len(vars_found))
        except (AssertionError, DSException), e:
            logger.error( 'SceneDataset.update_metadata error: %s', e.message)
            raise DSException(e.message)


    def IsGeographic(self):
        """
        Check whether the projection is geographic.

        This method should really be implemented on an a :py:class:`osgeo.osr.SpatialReference`, and hence would be
        accessed like

            >>> sr = osr.SpatialReference()
            >>> sr.ImportFromWkt(dataset.GetProjection())
            >>> sr.IsGeographic()

        However, I don't trust that the ``dataset.GetProjection()`` part will work for :py:class:`SceneDataset`s
        (and it almost certainly won't for EQR) and since we have a :py:class:`osr.SpatialReference` lying around,
        will access it directly.

        :todo:
            When someone fixes the projection issues on this class, they may want to change this (I reference it
            from :py:mod:`image_processor.nbar.radiative_transfer.run_tc`).
        """
        return self.spatial_ref.IsGeographic()


    def get_bool(self, instring, default=False):
        """
        Translate a string to a Boolean value.
        """
        if instring:
            if  type(instring) == str:
                return instring.lower() in ("yes", "true", "t", "1")
            else:
                return bool(instring)
        else:
            return default


    def orbital_time(self, irow, icol):

        gclat = self.geoc_lat(irow, icol)

        # Orbital arc distance from antinode (radians)
        rho = math.acos(math.sin(gclat) / self.satellite.INCL_SIN)

        # Antinode to pixel travel time (seconds)
        t = (rho + math.pi/2) / self.satellite.OMEGA - self.satellite.SWEEP_PERIOD / 2

        return t


    def bands(self, band_type='ALL'):
        """
        Returns a copy of an internal band list by type (read-only)
        """
        band_list = self._bands.get(band_type)
        if band_list:
            return list(band_list)
        else:
            return []


    def sensor_band_info(self, band_no):
        """
        Get a copy of sensor info dict for specified band number (read-only).

        :param band_no:
            The band got get the sensor dict for.
        :type band_no:
            int
        """
        assert self.satellite, 'Satellite not defined'
        return dict(self.satellite.BAND_LIST[band_no - 1])


    @property
    def bounds_getter(self):
        """
        A function that will produce a :py:class:`ULA3.utils.ImageShape` for this :py:class:`SceneDataset`.
        That function takes a SceneDataset as an argument.

        :todo:
            It may be better to have a function that just returns the bounds, rather than returning a function that
            will return bounds. One (perhaps silly) reason I did not do this is that this class already declares
            the function :py:meth:`SceneDataset.get_bounds`, which does something different. Having said that, I guess
            the presences of this function could still be confusing.
        """
        def gtr(dataset):
            assert isinstance(dataset, SceneDataset), "dataset must be a SceneDataset"
            # We need to do this differently than for gdal.Dataset, because the EQR datasets produced in GA
            # are kind of funky and GDAL doesn't handle them properly. Specifically, the degrees are scaled
            # by 1e5 or 1e6 and GDAL (which is apparently because LPGS was made to work in UTM and uses a format
            # like %0.3f for printing floats... which makes 'a mess' of numbers like 0.00025).
            g = self.geotransform
            return ImageShape(g[0], g[3], dataset.RasterXSize, dataset.RasterYSize, g[1], g[5], dataset.spatial_ref.ExportToWkt())
        return gtr


    @property
    def metadata_vars(self):
        """
        A copy of the list of variables read from metadata.
        """
        return list(self._metadata_vars)


    @property
    def data_dir(self):
        """
        Path of data directory.
        """
        return self._data_dir


    @property
    def root_dataset_pathname(self):
        """
        The path to the directory in which this dataset lives.
        """
        return self._root_dataset_pathname


    @property
    def root_band_number(self):
        """
        The root band number.
        """
        return self._root_band_number


    # Overridden Dataset properties
    @property
    def RasterCount(self):
        """
        The number of bands in the dataset.
        """
        return len(self._raster_dict)


    # Overridden Dataset methods
    def GetRasterBand(self, nBand):
        """Overrides Dataset method."""
        assert nBand in self._bands['ALL'], 'Illegal band number %s' % nBand
        band_dict = self._raster_dict[nBand]
        return band_dict['dataSet'].GetRasterBand(band_dict['rasterIndex'])


    def SetProjection(self, *args):
        """
        Set the projection for the dataset by setting the projection on all bands.
        """
        for dataset in self._sub_datasets:
            dataset.SetProjection(*args)


    def SetGeoTransform(self, *args):
        """
        Set the GeoTransform for this dataset by setting the given geogtransform on all bands.

        :todo:
            Ensure this function works with thermal and pan-chromatic bands.
        """
        for dataset in self._sub_datasets:
            dataset.SetGeoTransform(*args)


    def SetGCPs(self, *args):
        """
        Overrides Dataset method.

        :todo:
            Ensure this function works with thermal and pan-chromatic bands.
        """
        for dataset in self._sub_datasets:
            dataset.SetGCPs(*args);


    def FlushCache(self):
        """Overrides Dataset method"""
        for dataset in self._sub_datasets:
            dataset.FlushCache();


    def GetFileList(self):
        """Overrides Dataset method"""
        return list(self._file_list) # Return a copy of the private file list to prevent interference with it


    def GetSubDatasets(self):
        """Overrides Dataset method"""
        # Return a copy of the private subdataset list to prevent interference with it
        return list(self._sub_datasets)


    def ReadAsArray(self, *args, **kwargs):
        """
        Overrides Dataset method.

        :return:
            A 3D array with one level per reflective or thermal band.

        :todo:
            Allow optional selection of all bands including pan-chromatic
        """
        def get_array(dataset, template_ds):
            if dataset.RasterXSize == template_ds.RasterXSize and dataset.RasterYSize == template_ds.RasterYSize:
                return dataset.ReadAsArray()
            else:
                logger.debug('  Reprojecting %s to the same resolution as %s', dataset.GetFileList()[0], template_ds.GetFileList()[0])
                template_band = template_ds.GetRasterBand(1)
                driver = gdal.GetDriverByName('MEM')
                mem_dataset  = driver.Create("", template_ds.RasterXSize, template_ds.RasterYSize, 1,
                    template_band.DataType)
                mem_dataset.SetGeoTransform(template_ds.GetGeoTransform())
                mem_dataset.SetProjection(template_ds.GetProjection())
                _proj = gdal.ReprojectImage(dataset, mem_dataset)

                return mem_dataset.ReadAsArray()

        if self._sub_dataset_type == 'FST':
            # Return the 3D refelective dataset array for FST
            return self._root_dataset.ReadAsArray(*args, **kwargs)
        elif self._sub_dataset_type == 'TIF':
            datasets = [bd['dataSet'] for bd in self._raster_dict.values() if bd['band_file_number'] not in self.satellite.BAND_TYPES['PANCHROMATIC']]
            # Construct 3D array from first non-panchromatic subdataset
            result_array = numpy.resize(datasets[0].ReadAsArray(),
                                 (len(datasets),
                                  datasets[0].RasterYSize,
                                  datasets[0].RasterXSize))
            # Populate 3D array from all other reflective datasets
            for dataset in datasets[1:]:
                result_array[datasets.index(dataset)] = get_array(dataset, datasets[0])
            return result_array


    def GetMetadata(self, domain = ''):
        """
        Overrides Dataset method

        :param domain:
            Comma-separated string beginning with XML,MTL,REPORT,TIF or FST which contains the path
            to the required domain in a tree structure.

        :return:
            Nested dict reflecting tree structure under specified domain.
        """
        return self.GetMetadata_Dict(domain)


    def GetMetadataItem(self, name, domain = ''):
        """
        Overrides Dataset method
        N.B: domain and name are concatenated for the search

        :param name:
            Name of metadata item to find. Could be comma separated compound key

        :param domain:
            Comma-separated string beginning with XML,MTL,REPORT,TIF or FST which contains the key path to the
            specified domain in a tree-structured format
        """
        metadata = self._metadata.get_metadata(domain.split(',') + name.split(','))
        if metadata:
            if type(metadata) == dict:
                return None # domain + name is not a leaf node
            else:
                return metadata # domain + name is a leaf node
        else:
            return None # domain + name not found


    def GetMetadata_Dict(self, domain = ''):
        """
        Overrides Dataset method.

        :param name:
            Name of metadata item to find.

        :param domain:
            Comma-separated string beginning with XML,MTL,REPORT,TIF or FST which contains the path to
            the domain in a tree-structured format

        :return:
            Nested dict reflecting tree structure under specified domain
        """
        return self._metadata.get_metadata(domain.split(','))


    def GetMetadata_List(self, domain = ''):
        """
        Overrides Dataset method

        :param name:
            Name of metadata item to find.

        :param domain:
            Comma-separated string beginning with XML,MTL,REPORT,TIF or FST which contains the path to
            the domain in a tree-structured format

        :return:
            List of <key>, <Value> pairs for metadata under specified domain tree
        """

        return sorted(self._metadata.tree_to_list(domain))


    def SetMetadata(self, metadata, domain):
        """Overrides Dataset method"""
        # TODO: not sure how this one works
#        if self._eAccess != gdalconst.GA_ReadOnly:
#            self._metadata.set_metadata_value(domain.split(','), metadata, True)
        raise NotImplementedError()

    def SetMetadataItem(self, name, value, domain = ""):
        """Overrides Dataset method"""
        if self._eAccess != gdalconst.GA_ReadOnly:
            self._metadata.set_metadata_value(domain.split(',') + name.split(','), value, True)


    def img_coord(self, pixel, line, d_pixel=0.5, d_line=0.5):
        """
        Convert pixel indices to image coordinates.

        :param pixel:
            pixel index (x)

        :param line:
            line index (y)

        :param d_pixel:
            Fractional pixel width offset.

        :param d_line:
            Fractional pixel height offset.

        :return:
            Image coordinates (2-tuple)
        """

        g = self.geotransform

        return ( (g[0] + g[1]*pixel + g[2]*line) + g[1]*d_pixel,
                 (g[3] + g[4]*pixel + g[5]*line) + g[5]*d_line )


    def geo_coord_deg(self, pixel, line, d_pixel=0.5, d_line=0.5):
        """
        Convert pixel indices to geographic coordinates.

        :param pixel:
            Pixel (x) index.

        :param line:
            Line (y) index.

        :paran d_pixel:
            Pixel width offset.

        :param d_line:
            Pixel height offset.

        :return:
            Geographic coordinates (2-tuple, (lon, lat)) in decimal degrees.
        """

        x, y = self.img_coord(pixel, line, d_pixel, d_line)
        lon, lat, __z = self.cxform_to_geo.TransformPoint(x, y, 0)
        return (lon, lat)


    def geo_coord(self, pixel, line, d_pixel=0.5, d_line=0.5):
        """
        Convert pixel indices to geographic coordinates.

        :param pixel:
            Pixel (x) index.

        :param line:
            Line (y) index.

        :param d_pixel:
            Pixel width offset.

        :param d_line:
            Pixel height offset.

        :return:
            Geographic coordinates (2-tuple, (lon, lat)) in radians.
        """

        x, y = self.img_coord(pixel, line, d_pixel, d_line)
        lon, lat, __z = self.cxform_to_geo.TransformPoint(x, y, 0)
        return (math.radians(lon), math.radians(lat))


    def geo_lon(self, irow, icol):
        """
        Calculate geographic longitude at a (row, column) array location.

        :param irow:
            Row (y) index.

        :param icol:
            Column (x) index.

        :note:
            ``icol`` is the pixel/x index, and ``irow`` is the line/y index.
        """

        return self.geo_coord(icol, irow)[0]


    def geo_lat(self, irow, icol):
        """
        Calculate geographic latitude at a (row, column) array location.

        :param irow:
            Row (y) index.

        :param icol:
            Column (x) index.

        :note:
            ``icol`` is the pixel/x index, and ``irow`` is the line/y index.
        """

        return self.geo_coord(icol, irow)[1]


    def geo_elev(self, irow, icol):
        """
        Calculate geographic elevation at a (row, column) array location.

        :param irow:
            Row (y) index.

        :param icol:
            Column (x) index.

        :note:
            ``icol`` is the pixel/x index, and ``irow`` is the line/y index.
        """

        return self.geo_coord(icol, irow)[2]


    def geoc_lat(self, irow, icol):

        return geocentric.geocentric_lat(self.geo_lat(irow, icol))


    def coord_from_geo(self, lonlat):
        x, y, _z = self.cxform_from_geo.TransformPoint(lonlat[0], lonlat[1], 0)
        return (x, y)


    def get_bounds(self):
        """
        Get scene bounds for satellite grid calculations.

        :return:
            Coordinate tuple containing x and y coordinate ranges:
            ((<LLx>, <LRx>), (<LLy>, <ULy>))

        :notes:
            Function re-written so it returns the values given by GDAL, rather
            than rely on the MTL or XML files, which are unsuited for our
            purposes. 2013/08/16.
        """
        logger.debug('get_bounds() called', )
        if self.spatial_ref.IsProjected():
            return (
                     (self.coords['LL'][0], self.coords['LR'][0]),
                     (self.coords['LL'][1], self.coords['UL'][1])
                   )
            #return (
            #    (min(self.ll_x, self.ul_x), max(self.lr_x, self.ur_x)),
            #    (min(self.ll_y, self.lr_y), max(self.ul_y, self.ur_y)),
            #)
        else:
            return (
                     (self.lonlats['LL'][0], self.lonlats['LR'][0]),
                     (self.lonlats['LL'][1], self.lonlats['UL'][1])
                   )
            #return (
            #    (min(self.ll_lon, self.ul_lon), max(self.lr_lon, self.ur_lon)),
            #    (min(self.ll_lat, self.lr_lat), max(self.ul_lat, self.ur_lat)),
            #)


    def set_metadata_object(self, metadata_object):
        """
        Function to change private metadata object reference.
        """
        assert type(metadata_object) == Metadata
        logger.info('Metadata object changed')
        self._metadata = metadata_object

    def GetExtent(self):
        """
        Returns a list of corner coordinates from a GDAL geotransform object.
        List order [UL,LL,UR,LR]
        """
        gt = self.geotransform
        extents = []
        x_array = [0,self.RasterXSize]
        y_array = [0,self.RasterYSize]

        for px in x_array:
            for py in y_array:
                x = gt[0]+(px*gt[1])+(py*gt[2])
                y = gt[3]+(px*gt[4])+(py*gt[5])
                extents.append([x,y])
        return extents

    def GetProjectionRef(self):
        """
        Overrides the Dataset Class method only if the spatial_ref instance
        has been set.
        """
        if self.spatial_ref:
            return self.spatial_ref.ExportToWkt()
        else:
            return self._root_dataset.GetProjectionRef()

    def GetGeoTransform(self, *args, **kwargs):
        """
        Overrides the Dataset Class method only if the geotransform instance
        has been set.
        """
        if self.geotransform:
            return self.geotransform
        else:
            return self._root_dataset.GetGeoTransform(*args, **kwargs)
    
