"""
Utilities for the extraction of BRDF data.

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

import datetime
import logging
import math
import os
import re

import numpy
from osgeo import gdal
from osgeo import gdalconst
from osgeo import osr

from gaip import GriddedGeoBox
from gaip import write_img

logger = logging.getLogger('root.' + __name__)


#TODO: Implement resume
class BRDFLoaderError(Exception):
    """
    :todo:
        Someone who knows what this is used for should document it.
    """
    pass





class BRDFLoader(object):
    """
    Data loader for MCD43A1.005 mosaics.

    This class is used internally to load BRDF values from a file.
    """

    # MCD43A1.005 HDF files contain 3 subdatasets.
    SDS_MAP = { 0: 'brdf', 1: 'lat', 2:'lon' }

    # Format for accessing a SDS using GDAL.
    SDS_FORMAT = 'HDF4_SDS:UNKNOWN:"%s":%d'

    # Hardwired settings for first version of pre-MODIS BRDF database.
    # TODO create these in HDF file metadata when building database
    DEFAULTS = {
        'fill_value': -32768,
        'scale_factor': 0.001,
        'add_offset': 0.0,
    }

    def __init__(self, filename, UL=None, LR=None):
        """Initialise a BRDFLoader instance.

        Arguments:
            filename: data file name
            UL: (lon, lat) of ROI upper left corner [2-tuple, floats]
            LR: (lon, lat) of ROI lower right corner [2-tuple, floats]

        """

        self.filename = filename
        self.roi = {'UL': UL, 'LR': LR}

        logger.info('%s: filename=%s, roi=%s'
                    % (self.__class__.__name__, self.filename, str(self.roi)))

        if UL is None or LR is None:
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

        if (self.roi['UL'][0] < self.UL[0] or self.roi['LR'][0] > self.LR[0] or
            self.roi['UL'][1] > self.UL[1] or self.roi['LR'][1] < self.LR[1]):
            raise BRDFLoaderError(('%s: Region of interest %s extends beyond '
                                   'HDF domain {UL: %s, LR: %s}')
                                  % (self.__class__.__name__, str(self.roi),
                                     str(self.UL), str(self.LR)))

    def load(self):
        """
        Open file and load data arrays and required metadata.

        The following SDS structure is assumed in the HDF file:
            SDS 0: BRDF data array (2d)
            SDS 1: Latitude array (1d)
            SDS 2: Longitude array (1d)

        Latitude and longitude values are centre-of-pixel.

        Fill value, scale factor and offset are obtained from the HDF metadata.

        """

        # Load sub-datasets.

        for k in self.SDS_MAP:
            sds_file_spec = self.SDS_FORMAT % (self.filename, k)
            fd = gdal.Open(sds_file_spec, gdalconst.GA_ReadOnly)
            if fd is None:
                raise BRDFLoaderError('%s: gdal.Open failed [%s]'
                                      % (self.__class__.__name__, sds_file_spec))

            self.data[k] = fd.GetRasterBand(1).ReadAsArray()
            _type = type(self.data[k][0,0])

            logger.debug('%s: loaded sds=%d, type=%s, shape=%s'
                         % (self.__class__.__name__, k, str(_type),
                            str(self.data[k].shape)))

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

        logger.debug('%s: fill_value=%s, scale_factor=%s, add_offset=%s'
                     % (self.__class__.__name__, str(self.fill_value),
                        str(self.scale_factor), str(self.add_offset)))

    @property
    def delta_lon(self):
        """
        Get the longitude grid increment (cell size).

        :return:
            Longitude grid increment (decimal degrees).

        """

        return (self.data[2][0,1] - self.data[2][0,0])

    @property
    def delta_lat(self):
        """
        Get the latitude grid increment (cell size).

        :return:
            Latitude grid increment (decimal degrees).

        """

        return (self.data[1][0,1] - self.data[1][0,0])

    @property
    def UL(self):
        """
        Get the upper-left (NW) corner-of-pixel coordinates of the data.

        :return:
            Upper-left coordinate tuple: ``(lon, lat)`` (decimal degrees).

        """

        return (self.data[2][0,0] - self.delta_lon/2,
                self.data[1][0,0] - self.delta_lat/2)

    @property
    def LR(self):
        """
        Get the lower-right (SE) corner-of-pixel coordinates of the data.

        :return:
            Lower-right coordinate tuple: ``(lon, lat)`` (decimal degrees).

        """

        return (self.data[2][0,-1] + self.delta_lon/2,
                self.data[1][0,-1] + self.delta_lat/2)

    def mean_data_value(self):
        """
        Calculate the mean numeric BRDF value over the region of interest.

        :return:
            Mean data value (float) with scale and offset applied.

        """

        # Index calculation matches what happens in hdf_extractor.c.
        # TODO: verify correctness

        xmin = (self.roi['UL'][0] - self.UL[0]) / self.delta_lon
        xmax = (self.roi['LR'][0] - self.UL[0]) / self.delta_lon

        imin = max([0, int(math.ceil(xmin))])
        imax = min([self.data[0].shape[1], int(math.ceil(xmax))])

        ymin = (self.roi['UL'][1] - self.UL[1]) / self.delta_lat
        ymax = (self.roi['LR'][1] - self.UL[1]) / self.delta_lat

        jmin = max([0, int(math.ceil(ymin))])
        jmax = min([self.data[0].shape[0], int(math.ceil(ymax))])

        _data = numpy.ma.masked_values(self.data[0][jmin:jmax+1, imin:imax+1],
                                       self.fill_value).compressed()

        try:
            # Float the data sum (int16) to calculate the mean.
            dmean = float(numpy.sum(_data)) / _data.size
        except ZeroDivisionError:
            dmean = 0.0

        result = self.scale_factor * (dmean - self.add_offset)

        logger.debug(('%s: ROI=%s, imin=%d, imax=%d, xmin=%f, xmax=%f, '
                      'jmin=%d, jmax=%d, ymin=%f, ymax=%f, dmean=%.12f, '
                      'result=%.12f')
                     % (self.__class__.__name__, str(self.roi),
                        imin, imax, xmin, xmax,
                        jmin, jmax, ymin, ymax, dmean, result))

        return result

    def convert_format(self, filename, format='ENVI'):
        """
        Convert the HDF file to a more spatially recognisable data
        format such as ENVI or GTiff.
        The default format is ENVI (flat bianry file and an
        accompanying header (*.hdr) text file.
        """

        # Get the UL corner of the UL pixel co-ordinate
        ULlon = self.UL[0]
        ULlat = self.UL[1]

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
        geobox = GriddedGeoBox(shape=dims, origin=(ULlon, ULlat),
            pixelsize=res, crs=prj)

        # Write the file
        write_img(self.data[0], filename, format, geobox=geobox)

    def get_mean(self, array):
        """
        This mechanism will be used to calculate the mean in place in
        place of mean_data_value, which will still be kept until
        the results have been successfully validated.
        """

        valid = array != self.fill_value
        xbar = numpy.mean(array[valid])

        if not numpy.isfinite(xbar):
            xbar = 0.0
        else:
            xbar = self.scale_factor * (xbar - self.add_offset)

        return xbar



class BRDFLookupError(Exception):
    """
    :todo:
        Someone who knows what this is used for should document it.
    """
    pass





def get_brdf_dirs_modis(brdf_root, scene_date, pattern='\d{4}.\d{2}.\d{2}$'):
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
        # Returns interval midpoint date of a MCD43A1.005/YYYY.MM.DD directory.
        return datetime.date(*[int(x) for x in s.split(sep)]) + offset

    # Compile the search pattern
    BRDF_DIR_PATTERN = re.compile(pattern)

    # List only directories that match 'YYYY.MM.DD' format.
    dirs = sorted([d for d in os.listdir(brdf_root) if BRDF_DIR_PATTERN.match(d)])

    # Find the N (n_dirs) BRDF directories with midpoints closest to the
    # scene date.
    delta_map = { abs(parsedate(x) - scene_date): x for x in dirs }

    if scene_date < (parsedate(dirs[0]) - offset):
        raise BRDFLookupError('scene date precedes first MODIS date (%s)' % dirs[0])

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
        min( abs(int(x) + offset - i),
             abs(db_doy_max - abs(int(x) + offset - i)) ) : x for x in dirs
    }

    # Return the closest match (the zeroth index)
    result = delta_map[sorted(delta_map)[0]]

    return result

