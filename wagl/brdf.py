#!/usr/bin/env python

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

from __future__ import absolute_import, print_function
import subprocess
import datetime
import logging
import math
import os
from os.path import join as pjoin, basename
import tempfile
import re
from urllib.parse import urlparse
import numpy as np

from osgeo import gdal
from osgeo import gdalconst
from osgeo import osr
from shapely.geometry import Polygon
from shapely import wkt
from wagl.constants import BrdfParameters
from wagl.geobox import GriddedGeoBox
from wagl.hdf5 import H5CompressionFilter, write_h5_image
from wagl.metadata import extract_ancillary_metadata

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

        log.debug('%s: filename=%s, roi=%s', self.__class__.__name__,
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

        if not intersection.bounds:
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
            self.roi_polygon = roi_poly

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

    def convert_format(self, dataset_name, group, attrs=None,
                       compression=H5CompressionFilter.LZF, filter_opts=None):
        """
        Convert the HDF file to a HDF5 dataset.
        """
        if attrs is None:
            attrs = {}

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
        res = (abs(pixsz_x), abs(pixsz_y))
        geobox = GriddedGeoBox(shape=dims, origin=(ul_lon, ul_lat),
                               pixelsize=res, crs=prj)

        # Write the dataset
        attrs['description'] = 'Converted BRDF data from H4 to H5.'
        attrs['crs_wkt'] = prj
        attrs['geotransform'] = geobox.transform.to_gdal()
        write_h5_image(self.data[0], dataset_name, group, compression, attrs,
                       filter_opts)

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


def _date_proximity(cmp_date, date_interpreter=lambda x: x):
    """_date_proximity providers a comparator for an interable
    with an interpreter function. Used to find the closest item
    in a list.

    If two dates are equidistant return the most recent.

    :param cmp_date: date to compare list against
    :param date_interprater: function applied to the list to
        transform items into dates
    """
    def _proximity_comparator(date):
        _date = date_interpreter(date)
        return (
            abs(_date - cmp_date), 
            -1 * _date.year, 
            -1 * _date.month, 
            -1 * _date.day
        )

    return _proximity_comparator


def get_brdf_dirs_modis(brdf_root, scene_date, pattern='%Y.%m.%d'):
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
        A string handed to strptime to interpret directory names into
        observation dates for the brdf ancillary.
    :type pattern:
        :py:class:`str`

    :type pattern:
        :py:class:`str`

    :return:
       A string containing the closest matching BRDF directory.

    """

    # MCD43A1.005 db interval half-width (days).
    offset = datetime.timedelta(8)
    _offset_scene_date = scene_date - offset

    dirs = []
    for dname in sorted(os.listdir(brdf_root)):
        try:
            dirs.append(datetime.datetime.strptime(dname, pattern).date())
        except ValueError:
            pass  # Ignore directories that don't match specified pattern

    return min(dirs, key=_date_proximity(_offset_scene_date)).strftime(pattern)


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
    offset = datetime.timedelta(8)

    # Find the N (=n_dirs) BRDF directories with midpoints closest to the
    # scene date.
    # Pre-MODIS BRDF directories are named 'XXX' (day-of-year).
    # Return a list of n_dirs directories to maintain compatibility with
    # the NBAR code, even though we know that the nearest day-of-year
    # database dir will contain usable data.

    _offset_scene_date = scene_date - offset

    # Build list of dates for comparison
    dir_dates = []

    # Standardise names be prepended with leading zeros
    for doy in sorted(os.listdir(brdf_root), key=lambda x: x.zfill(3)):
        dir_dates.append((str(_offset_scene_date.year), doy))

    # Add boundary entry for previous year
    dir_dates.insert(0, (str(_offset_scene_date.year - 1), dir_dates[-1][1]))
    # Add boundary entry for next year accounting for inserted entry
    dir_dates.append((str(_offset_scene_date.year + 1), dir_dates[1][1]))

    # Interpreter function
    doy_intpr = lambda x: datetime.datetime.strptime(' '.join(x), '%Y %j').date()

    # return directory name without year
    return min(dir_dates, key=_date_proximity(_offset_scene_date, doy_intpr))[1]


def get_brdf_data(acquisition, brdf_primary_path, brdf_secondary_path,
                  compression=H5CompressionFilter.LZF, filter_opts=None):
    """
    Calculates the mean BRDF value for the given acquisition,
    for each BRDF parameter ['geo', 'iso', 'vol'] that covers
    the acquisition's extents.

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

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        A `dict` with the keys:

            * BrdfParameters.ISO
            * BrdfParameters.VOL
            * BrdfParameters.GEO

        Values for each BRDF Parameter are accessed via the key named
        `value`.

    :notes:
        The keywords compression and filter_opts aren't used as we no
        longer save the BRDF imagery. However, we may need to store
        tables in future, therefore they can remain until we know
        for sure they'll never be used.
    """
    def find_file(files, brdf_wl, parameter):
        """Find file with a specific name."""
        for f in files:
            if f.find(brdf_wl) != -1 and f.find(parameter) != -1:
                return f
        return None

    # Compute the geobox
    geobox = acquisition.gridded_geo_box()

    # Get the date of acquisition
    dt = acquisition.acquisition_datetime.date()

    # Compare the scene date and MODIS BRDF start date to select the
    # BRDF data root directory.
    # Scene dates outside the range of the CSIRO mosaic data
    # should use the pre-MODIS, Jupp-Li BRDF.
    brdf_dir_list = sorted(os.listdir(brdf_primary_path))

    try:
        brdf_dir_range = [brdf_dir_list[0], brdf_dir_list[-1]]
        brdf_range = [datetime.date(*[int(x) for x in y.split('.')])
                      for y in brdf_dir_range]

        use_JuppLi_brdf = (dt < brdf_range[0] or dt > brdf_range[1])
    except IndexError:
        use_JuppLi_brdf = True  # use JuppLi if no primary data available

    if use_JuppLi_brdf:
        brdf_base_dir = brdf_secondary_path
        brdf_dirs = get_brdf_dirs_pre_modis(brdf_base_dir, dt)
    else:
        brdf_base_dir = brdf_primary_path
        brdf_dirs = get_brdf_dirs_modis(brdf_base_dir, dt)

    # The following hdflist code was resurrected from the old SVN repo. JS
    # get all HDF files in the input dir
    dbDir = pjoin(brdf_base_dir, brdf_dirs)
    three_tup = os.walk(dbDir)
    hdflist = []
    hdfhome = None

    for (hdfhome, _, filelist) in three_tup:
        for f in filelist:
            if f.endswith(".hdf.gz") or f.endswith(".hdf"):
                hdflist.append(f)

    results = {}
    for param in BrdfParameters:
        hdf_fname = find_file(hdflist, acquisition.brdf_wavelength,
                              param.name.lower())

        hdfFile = pjoin(hdfhome, hdf_fname)

        # Test if the file exists and has correct permissions
        try:
            with open(hdfFile, 'rb') as f:
                pass
        except IOError:
            print("Unable to open file %s" % hdfFile)

        with tempfile.TemporaryDirectory() as tmpdir:
            # Unzip if we need to
            if hdfFile.endswith(".hdf.gz"):
                hdf_file = pjoin(tmpdir, re.sub(".hdf.gz", ".hdf",
                                                basename(hdfFile)))
                cmd = "gunzip -c %s > %s" % (hdfFile, hdf_file)
                subprocess.check_call(cmd, shell=True)
            else:
                hdf_file = hdfFile

            # Load the file
            brdf_object = BRDFLoader(hdf_file, ul=geobox.ul_lonlat,
                                     lr=geobox.lr_lonlat)

            # guard against roi's that don't intersect
            if not brdf_object.intersects:
                msg = "ROI is outside the BRDF extents!"
                log.error(msg)
                raise Exception(msg)

            # calculate the mean value
            brdf_mean_value = brdf_object.mean_data_value()

        # Add the brdf filename and mean value to brdf_dict
        url = urlparse(hdfFile, scheme='file').geturl()
        res = {'data_source': 'BRDF',
               'url': url,
               'value': brdf_mean_value}

        # ancillary metadata tracking
        md = extract_ancillary_metadata(hdfFile)
        for key in md:
            res[key] = md[key]

        results[param] = res

    # check for no brdf (iso, vol, geo) (0, 0, 0) and convert to (1, 0, 0)
    # and strip any file level metadata
    if all([v['value'] == 0 for _, v in results.items()]):
        results[BrdfParameters.ISO] = {'value': 1.0}
        results[BrdfParameters.VOL] = {'value': 0.0}
        results[BrdfParameters.GEO] = {'value': 0.0}

    # add very basic brdf description metadata and the roi polygon
    for param in BrdfParameters:
        results[param]['extents'] = wkt.dumps(brdf_object.roi_polygon)

    return results
