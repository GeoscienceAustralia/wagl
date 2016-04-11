"""
BRDF data extraction utilities
------------------------------

The :ref:`nbar-algorithm-label` and :ref:`tc-algorithm-label` algorithms
require estimates of various atmospheric parameters, which are produced using
`MODTRAN <http://modtran5.com/>`_. MODTRAN, in turn, requires `BRDF
<http://en.wikipedia.org/wiki/Bidirectional_reflectance_distribution_function>`_
estimates. The estimates used in the ULA, are based on `MODIS
<http://modis.gsfc.nasa.gov/>`_ and are produced by CSIRO. For more
information, on how these are used, see :download:`this
<auxiliary/li_etal_2010_05422912.pdf>`.

`MODIS <http://modis.gsfc.nasa.gov/>`_, pre Feb 2001, MODIS data was not
available and an alternative method of deriving `BRDF
<http://en.wikipedia.org/wiki/Bidirectional_reflectance_distribution_function>`_
estimates is required.

"""

import subprocess
import datetime
import logging
import math
import os
import re
import numpy as np

from osgeo import gdal
from osgeo import gdalconst
from osgeo import osr
from shapely.geometry import Polygon
from gaip import GriddedGeoBox
from gaip import write_img
from gaip import constants
from gaip import read_subset
from gaip import extract_ancillary_metadata

log = logging.getLogger('root.' + __name__)


class BRDFLoaderError(Exception):

    """
    BRDF Loader Error
    """
    pass


class BRDFLookupError(Exception):

    """
    BRDF Lookup Error
    """
    pass


class BRDFLoader(object):

    """
    Data loader for MCD43A1.005 mosaics.

    This class is used internally to load BRDF values from a file.
    """

    # MCD43A1.005 HDF files contain 3 subdatasets.
    SDS_MAP = {0: 'brdf', 1: 'lat', 2: 'lon'}

    # Format for accessing a SDS using GDAL.
    SDS_FORMAT = 'HDF4_SDS:UNKNOWN:"%s":%d'

    # Hardwired settings for first version of pre-MODIS BRDF database.
    DEFAULTS = {
        'fill_value': -32768,
        'scale_factor': 0.001,
        'add_offset': 0.0,
    }

    def __init__(self, filename, ul=None, lr=None):
        """Initialise a BRDFLoader instance.

        Arguments:
            filename: data file name
            ul: (lon, lat) of ROI upper left corner [2-tuple, floats]
            lr: (lon, lat) of ROI lower right corner [2-tuple, floats]

        """

        self.filename = filename
        self.roi = {'UL': ul, 'LR': lr}

        log.info('%s: filename=%s, roi=%s', self.__class__.__name__,
                 self.filename, str(self.roi))

        if ul is None or lr is None:
            raise BRDFLoaderError('%s: UL and/or LR not defined'
                                  % (self.__class__.__name__,))

        # Initialise data array container.

        self.data = {}

        # Initialise metadata values.

        self.fill_value = self.DEFAULTS['fill_value']
        self.scale_factor = self.DEFAULTS['scale_factor']
        self.scale_factor_err = self.DEFAULTS['scale_factor']
        self.add_offset = self.DEFAULTS['add_offset']
        self.add_offset_err = self.DEFAULTS['add_offset']

        # Load data from HDF file.

        self.load()

        # The region-of-interest (scene) should lie within the HDF extents.

        # New feature: Polygon intersection
        coords = [self.ul, (self.lr[0], self.ul[1]),
                  self.lr, (self.ul[0], self.lr[1])]
        brdf_poly = Polygon(coords)
        coords = [self.roi['UL'], (self.roi['LR'][0], self.roi['UL'][1]),
                  self.roi['LR'], (self.roi['UL'][0], self.roi['LR'][1])]
        roi_poly = Polygon(coords)

        intersection = brdf_poly.intersection(roi_poly)

        if len(intersection.bounds) == 0:
            self.intersects = False
            raise BRDFLoaderError(('%s: Region of interest %s extends beyond '
                                   'HDF domain {UL: %s, LR: %s}')
                                  % (self.__class__.__name__, str(self.roi),
                                     str(self.ul), str(self.lr)))
        else:
            self.intersects = True
            i_ul = (intersection.bounds[0], intersection.bounds[-1])
            i_lr = (intersection.bounds[2], intersection.bounds[1])
            self.roi = {'UL': i_ul, 'LR': i_lr}


    def load(self):
        """
        Open file and load data arrays and required metadata.

        The following SDS structure is assumed in the HDF file:
            SDS 0: BRDF data array (2d)
            SDS 1: Latitude array (1d)
            SDS 2: Longitude array (1d)

        Latitude and longitude values are centre-of-pixel.

        Fill value, scale factor and offset are obtained from the HDF
        metadata.

        """
        # Load sub-datasets.
        for k in self.SDS_MAP:
            sds_file_spec = self.SDS_FORMAT % (self.filename, k)
            fd = gdal.Open(sds_file_spec, gdalconst.GA_ReadOnly)
            if fd is None:
                raise BRDFLoaderError('%s: gdal.Open failed [%s]'
                                      % (self.__class__.__name__,
                                         sds_file_spec))

            self.data[k] = fd.GetRasterBand(1).ReadAsArray()
            _type = type(self.data[k][0, 0])

            log.debug('%s: loaded sds=%d, type=%s, shape=%s',
                      self.__class__.__name__, k, str(_type),
                      str(self.data[k].shape))

            # Populate metadata entries after reading the BRDF data
            # array (SDS 0).

            if k == 0:
                m = fd.GetMetadata_Dict()
                if m:
                    # cast string to dtype
                    self.fill_value = _type(m['_FillValue'])
                    # floats...
                    self.scale_factor = float(m['scale_factor'])
                    self.scale_factor_err = float(m['scale_factor_err'])
                    self.add_offset = float(m['add_offset'])
                    self.add_offset_err = float(m['add_offset_err'])
                else:
                    # cast default fill value to dtype
                    self.fill_value = _type(self.fill_value)

            fd = None

        log.debug('%s: fill_value=%s, scale_factor=%s, add_offset=%s',
                  self.__class__.__name__, str(self.fill_value),
                  str(self.scale_factor), str(self.add_offset))

    @property
    def delta_lon(self):
        """
        Get the longitude grid increment (cell size).

        :return:
            Longitude grid increment (decimal degrees).

        """

        return (self.data[2][0, 1] - self.data[2][0, 0])

    @property
    def delta_lat(self):
        """
        Get the latitude grid increment (cell size).

        :return:
            Latitude grid increment (decimal degrees).

        """

        return (self.data[1][0, 1] - self.data[1][0, 0])

    @property
    def ul(self):
        """
        Get the upper-left (NW) corner-of-pixel coordinates of the data.

        :return:
            Upper-left coordinate tuple: ``(lon, lat)`` (decimal degrees).

        """

        return (self.data[2][0, 0] - self.delta_lon / 2,
                self.data[1][0, 0] - self.delta_lat / 2)

    @property
    def lr(self):
        """
        Get the lower-right (SE) corner-of-pixel coordinates of the data.

        :return:
            Lower-right coordinate tuple: ``(lon, lat)`` (decimal degrees).

        """

        return (self.data[2][0, -1] + self.delta_lon / 2,
                self.data[1][0, -1] + self.delta_lat / 2)

    def mean_data_value(self):
        """
        Calculate the mean numeric BRDF value over the region of interest.

        :return:
            Mean data value (float) with scale and offset applied.

        """

        # Index calculation matches what happens in hdf_extractor.c.
        # TODO: verify correctness

        xmin = (self.roi['UL'][0] - self.ul[0]) / self.delta_lon
        xmax = (self.roi['LR'][0] - self.ul[0]) / self.delta_lon

        imin = max([0, int(math.ceil(xmin))])
        imax = min([self.data[0].shape[1], int(math.ceil(xmax))])

        ymin = (self.roi['UL'][1] - self.ul[1]) / self.delta_lat
        ymax = (self.roi['LR'][1] - self.ul[1]) / self.delta_lat

        jmin = max([0, int(math.ceil(ymin))])
        jmax = min([self.data[0].shape[0], int(math.ceil(ymax))])

        data = np.ma.masked_values(self.data[0][jmin:jmax + 1, imin:imax + 1],
                                   self.fill_value).compressed()

        try:
            # Float the data sum (int16) to calculate the mean.
            dmean = float(np.sum(data)) / data.size
        except ZeroDivisionError:
            dmean = 0.0

        result = self.scale_factor * (dmean - self.add_offset)

        log.debug('%s: ROI=%s, imin=%d, imax=%d, xmin=%f, xmax=%f, '
                  'jmin=%d, jmax=%d, ymin=%f, ymax=%f, dmean=%.12f, '
                  'result=%.12f',
                  self.__class__.__name__, str(self.roi), imin, imax, xmin,
                  xmax, jmin, jmax, ymin, ymax, dmean, result)

        return result

    def convert_format(self, filename, fmt='ENVI'):
        """
        Convert the HDF file to a more spatially recognisable data
        format such as ENVI or GTiff.
        The default format is ENVI (flat bianry file and an
        accompanying header (*.hdr) text file.
        """

        # Get the UL corner of the UL pixel co-ordinate
        ul_lon = self.ul[0]
        ul_lat = self.ul[1]

        # pixel size x & y
        pixsz_x = self.delta_lon
        pixsz_y = self.delta_lat

        # Setup the projection; assuming Geographics WGS84
        # (Tests have shown that this appears to be the case)
        # (unfortunately it is not expicitly defined in the HDF file)
        sr = osr.SpatialReference()
        sr.SetWellKnownGeogCS("WGS84")
        prj = sr.ExportToWkt()

        # Setup the geobox
        dims = self.data[0].shape
        res = (pixsz_x, pixsz_y)
        geobox = GriddedGeoBox(shape=dims, origin=(ul_lon, ul_lat),
                               pixelsize=res, crs=prj)

        # Write the file
        write_img(self.data[0], filename, fmt, geobox=geobox)

    def get_mean(self, array):
        """
        This mechanism will be used to calculate the mean in place in
        place of mean_data_value, which will still be kept until
        the results have been successfully validated.
        """

        valid = array != self.fill_value
        xbar = np.mean(array[valid])

        if not np.isfinite(xbar):
            xbar = 0.0
        else:
            xbar = self.scale_factor * (xbar - self.add_offset)

        return xbar


def get_brdf_dirs_modis(brdf_root, scene_date, pattern=r'\d{4}.\d{2}.\d{2}$'):
    """
    Get list of MODIS BRDF directories for the dataset.

    :param brdf_root:
        BRDF root directory.
    :type brdf_root:
        :py:class:`str`

    :param scene_date:
        Scene Date.
    :type scene_date:
        :py:class:`datetime.date`

    :param pattern:
        A string containing the pattern upon which directories should
        be matched to. Default is '.' to match any charachter.
        Regular expression (re module) will be used for string
        matching.
    :type pattern:
        :py:class:`str`

    :return:
       A string containing the closest matching BRDF directory.

    """

    # MCD43A1.005 db interval half-width (days).
    offset = datetime.timedelta(8)

    def parsedate(s, sep='.'):
        """
        Returns interval midpoint date of a MCD43A1.005/YYYY.MM.DD directory.
        """
        return datetime.date(*[int(x) for x in s.split(sep)]) + offset

    # Compile the search pattern
    brdf_dir_pattern = re.compile(pattern)

    # List only directories that match 'YYYY.MM.DD' format.
    dirs = sorted([d for d in os.listdir(brdf_root)
                   if brdf_dir_pattern.match(d)])

    # Find the N (n_dirs) BRDF directories with midpoints closest to the
    # scene date.
    delta_map = {abs(parsedate(x) - scene_date): x for x in dirs}

    if scene_date < (parsedate(dirs[0]) - offset):
        raise BRDFLookupError('scene date precedes first MODIS date (%s)' \
                              % dirs[0])

    # Return the closest match (the zeroth index)
    result = delta_map[sorted(delta_map)[0]]

    return result


def get_brdf_dirs_pre_modis(brdf_root, scene_date):
    """
    Get list of pre-MODIS BRDF directories for the dataset.

    :param brdf_root:
        BRDF root directory.
    :type brdf_root:
        :py:class:`str`

    :param scene_date:
        Scene Date.
    :type scene_date:
        :py:class:`datetime.date`

    :return:
       A string containing the closest matching BRDF directory.

    """
    # Pre-MODIS db interval half-width (days).
    offset = 8

    # Find the N (=n_dirs) BRDF directories with midpoints closest to the
    # scene date.
    # Pre-MODIS BRDF directories are named 'XXX' (day-of-year).
    # Return a list of n_dirs directories to maintain compatibility with
    # the NBAR code, even though we know that the nearest day-of-year
    # database dir will contain usable data.

    dirs = sorted(os.listdir(brdf_root))
    scene_doy = scene_date.strftime('%j')
    i = int(scene_doy)
    db_doy_max = max(range(1, 365, 8))
    delta_map = {
        min(abs(int(x) + offset - i),
            abs(db_doy_max - abs(int(x) + offset - i))): x for x in dirs
    }

    # Return the closest match (the zeroth index)
    result = delta_map[sorted(delta_map)[0]]

    return result


def get_brdf_data(acquisition, brdf_primary_path, brdf_secondary_path,
                  work_path):
    """
    Calculates the mean BRDF value for each band wavelength of your
    sensor, for each BRDF factor ['geo', 'iso', 'vol'] that covers
    your image extents.

    :param acquisition:
        An instance of an acquisitions object.

    :param brdf_primary_path:
        A string containing the full file system path to your directory
        containing the source BRDF files.  The BRDF directories are
        assumed to be yyyy.mm.dd naming convention.

    :param brdf_secondary_path:
        A string containing the full file system path to your directory
        containing the Jupp-Li backup BRDF data.  To be used for
        pre-MODIS and potentially post-MODIS acquisitions.

    :param work_path:
        A string containing the full file system path to your NBAR
        working directory. Intermediate BRDF files will be saved to
        work_path/brdf_intermediates/.

    :return:
        A dictionary with tuple (band, factor) as the keys. Each key
        represents the band of your satllite/sensor and brdf factor.
        Each key contains a dictionary with the following keys:

            * data_source -> BRDF
            * data_file -> File system path to the location of the selected
              BRDF wavelength and factor combination.
            * value -> The mean BRDF value covering your image extents.
    """
    # Retrieve the satellite and sensor for the acquisition
    satellite = acquisition.spacecraft_id
    sensor = acquisition.sensor_id

    # Get the required BRDF LUT & factors list
    nbar_constants = constants.NBARConstants(satellite, sensor)

    brdf_lut = nbar_constants.get_brdf_lut()
    brdf_factors = nbar_constants.get_brdf_factors()

    # Compute the geobox
    geobox = acquisition.gridded_geo_box()

    # Get the date of acquisition
    dt = acquisition.scene_center_datetime.date()

    # Get the boundary extents of the image
    nw = geobox.ul_lonlat
    se = geobox.lr_lonlat

    # Compare the scene date and MODIS BRDF start date to select the
    # BRDF data root directory.
    # Scene dates outside the range of the CSIRO mosaic data
    # (currently 2000-02-18 through 2013-01-09) should use the pre-MODIS,
    # Jupp-Li BRDF.
    brdf_dir_list = sorted(os.listdir(brdf_primary_path))
    brdf_dir_range = [brdf_dir_list[0], brdf_dir_list[-1]]
    brdf_range = [datetime.date(*[int(x) for x in y.split('.')])
                  for y in brdf_dir_range]

    use_JuppLi_brdf = (dt < brdf_range[0] or dt > brdf_range[1])

    if use_JuppLi_brdf:
        brdf_base_dir = brdf_secondary_path
        brdf_dirs = get_brdf_dirs_pre_modis(brdf_base_dir, dt)
    else:
        brdf_base_dir = brdf_primary_path
        brdf_dirs = get_brdf_dirs_modis(brdf_base_dir, dt)

    # The following hdfList code was resurrected from the old SVN repo. JS
    # get all HDF files in the input dir
    dbDir = os.path.join(brdf_base_dir, brdf_dirs)
    three_tup = os.walk(dbDir)
    hdfList = []
    for (hdfHome, _, filelist) in three_tup:
        for f in filelist:
            if f.endswith(".hdf.gz") or f.endswith(".hdf"):
                hdfList.append(f)

    # Initialise the brdf dictionary to store the results
    brdf_dict = {}

    # Create a BRDF directory in the work path to store the intermediate
    # files such as format conversion and subsets.
    brdf_out_path = os.path.join(work_path, 'brdf_intermediates')
    if not os.path.exists(brdf_out_path):
        os.makedirs(brdf_out_path)

    def find_file(files, band_wl, factor):
        """Find file with a specific name."""
        for f in files:
            if f.find(band_wl) != -1 and f.find(factor) != -1:
                return f
        return None

    # Loop over each defined band and each BRDF factor
    for band in brdf_lut.keys():
        bandwl = brdf_lut[band]  # Band wavelength
        for factor in brdf_factors:
            hdfFileName = find_file(hdfList, bandwl, factor)

            hdfFile = os.path.join(hdfHome, hdfFileName) #FIXME: hdfHome

            # Test if the file exists and has correct permissions
            try:
                with open(hdfFile, 'rb') as f:
                    pass
            except IOError:
                print "Unable to open file %s" % hdfFile

            # Unzip if we need to
            if hdfFile.endswith(".hdf.gz"):
                hdf_file = os.path.join(
                    work_path,
                    re.sub(".hdf.gz", ".hdf",
                           os.path.basename(hdfFile)))
                cmd = "gunzip -c %s > %s" % (hdfFile, hdf_file)
                subprocess.check_call(cmd, shell=True)
            else:
                hdf_file = hdfFile

            # the following now converts the file format and outputs a subset.
            # this should proove useful for debugging and testing.

            # Load the file
            brdf_object = BRDFLoader(hdf_file, ul=nw, lr=se)

            # setup the output filename
            out_fname = '_'.join(['Band', str(band), bandwl, factor])
            out_fname = os.path.join(brdf_out_path, out_fname)

            # Convert the file format
            brdf_object.convert_format(out_fname)

            # Get the intersected roi
            # the intersection is used rather than the actual bounds,
            # as we want to include partial overlaps rather than exclude them
            if not brdf_object.intersects:
                msg = "ROI is outside the BRDF extents!"
                log.error(msg)
                raise Exception(msg)

            roi = brdf_object.roi
            ul_lon, ul_lat = roi['UL']
            ur_lon, ur_lat = (roi['LR'][0], roi['UL'][1])
            lr_lon, lr_lat = roi['UL']
            ll_lon, ll_lat = (roi['UL'][0], roi['LR'][1])

            # Read the subset and geotransform that corresponds to the subset
            subset, geobox_subset = read_subset(out_fname,
                                                (ul_lon, ul_lat),
                                                (ur_lon, ur_lat),
                                                (lr_lon, lr_lat),
                                                (ll_lon, ll_lat))

            # The brdf_object has the scale and offsets so calculate the mean
            # through the brdf_object
            brdf_mean_value = brdf_object.get_mean(subset)

            # Output the brdf subset
            out_fname_subset = out_fname + '_subset'
            write_img(subset, out_fname_subset, geobox=geobox_subset)

            # Remove temporary unzipped file
            if hdf_file.find(work_path) == 0:
                os.remove(hdf_file)

            # Add the brdf filename and mean value to brdf_dict
            res = {'data_source': 'BRDF',
                   'data_file': hdfFile,
                   'value': brdf_mean_value}

            # ancillary metadata tracking
            md = extract_ancillary_metadata(hdfFile)
            for key in md:
                res[key] = md[key]

            brdf_dict[(band, factor)] = res

    # check for no brdf (iso, vol, geo) (0, 0, 0) and convert to (1, 0, 0)
    for band in brdf_lut.keys():
            data = {}
        for factor in brdf_factors:
            data[factor] =  brdf_dict[(band, factor)]
        if all([i == 0 for i in data.values()]):
            brdf_dict[(band, 'iso')] = 1.0

    return brdf_dict
