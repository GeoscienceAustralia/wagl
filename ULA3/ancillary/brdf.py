"""
Utilities for the extraction of BRDF data.

The :ref:`nbar-algorithm-label` and :ref:`tc-algorithm-label` algorithms require estimates of various atmospheric
parameters, which are produced using `MODTRAN <http://modtran5.com/>`_. MODTRAN, in turn, requires
`BRDF <http://en.wikipedia.org/wiki/Bidirectional_reflectance_distribution_function>`_ estimates. The estimates used
in the ULA, are based on `MODIS <http://modis.gsfc.nasa.gov/>`_ and are produced by CSIRO. For more information, on
how these are used, see :download:`this <auxiliary/li_etal_2010_05422912.pdf>`.

`MODIS <http://modis.gsfc.nasa.gov/>`_, pre Feb 2001, MODIS data was not available and an alternative method of
deriving `BRDF <http://en.wikipedia.org/wiki/Bidirectional_reflectance_distribution_function>`_ estimates is required.

:todo:
    Someone who knows more about this (particularly the pre 2001 stuff) should document it.
"""

import logging, os, re, commands, datetime, math, numpy
from glob import glob
from osgeo import gdal, gdalconst

logger = logging.getLogger('root.' + __name__)

# keywords for BRDF parameters f0, f1, f2
_FACTOR_LIST = ['geo', 'iso', 'vol']

# Expect BRDF directories of the form 'YYYY.MM.DD'
_BRDF_DIR_PATTERN = re.compile('^\d{4}\.\d{2}.\d{2}$')

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





class BRDFLookupError(Exception):
    """
    :todo:
        Someone who knows what this is used for should document it.
    """
    pass





def get_brdf_dirs_modis(brdf_root, scene_date, n_dirs=2):
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

    :param n_dirs:
        The number of directories to be found (primary + secondaries).
    :type n_dirs:
        :py:class:`int`

    :return:
       List of BRDF directories, [<primary>, <secondary>, ...]

    """

    # MCD43A1.005 db interval half-width (days).
    offset = datetime.timedelta(8)

    def __date(s, sep='.'):
        # Returns interval midpoint date of a MCD43A1.005/YYYY.MM.DD directory.
        return datetime.date(*[int(x) for x in s.split(sep)]) + offset

    # List only directories that match 'YYYY.MM.DD' format.
    dirs = sorted([d for d in os.listdir(brdf_root) if _BRDF_DIR_PATTERN.match(d)])

    # Find the N (n_dirs) BRDF directories with midpoints closest to the
    # scene date.
    delta_map = { abs(__date(x) - scene_date): x for x in dirs }

    if scene_date < (__date(dirs[0]) - offset):
        raise BRDFLookupError('scene date precedes first MODIS date (%s)' % dirs[0])

    return [delta_map[k] for k in sorted(delta_map)[:n_dirs]]





def get_brdf_dirs_pre_modis(brdf_root, scene_date, n_dirs=3):
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

    :param n_dirs:
        Number of directories to be found (primary + secondaries).
    :type n_dirs:
        :py:class:`int`

    :return:
       List of BRDF directories, [<primary>, <secondary>, ...]

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

    return [delta_map[k] for k in sorted(delta_map)[:n_dirs]]





def average_brdf_value(
    ll_lat, ll_lon,
    ur_lat, ur_lon,
    band_string,
    wavelength_range,
    factor,
    work_path,
    brdf_root, brdf_dir):
    """
    Get the average value of subset of a HDF file.

    :param ll_lat:
        The latitude of the lower left corner of the region ('ll' for 'Lower Left').
    :type ll_lat:
        :py:class:`float`

    :param ll_lon:
        The longitude of the lower left corner of the region ('ll' for 'Lower Left').
    :type ll_lon:
        :py:class:`float`

    :param ur_lat:
        The latitude of the upper right corner of the region ('ur' for 'Upper Right').
    :type ur_lat:
        :py:class:`float`

    :param ur_lon:
        The longitude of the upper right corner of the region ('ur' for 'Upper Right').
    :type ur_lon:
        :py:class:`float`

    :param band_string:
        1-based band number string (i.e. in ['1','2','3','4','5','7']).
    :type band_string:
        :py:class:`str`

    :param wavelength_range:
        Not sure.
    :type wavelength_range:
        Not sure

    :param factor:
        Not sure.
    :type factor:
        Not sure

    :param work_path:
        The working directory.
    :type work_path:
        :py:class:`str`

    :param brdf_root:
        Not sure.
    :type brdf_root:
        :py:class:`str`

    :param brdf_dir:
        Not sure.
    :type brdf_dir:
        :py:class:`str`

    :note:
        This uses maximal axis-aligned extents for BRDF mean value calculation. Note that latitude
        min-max logic is valid for the Southern hemisphere only.

    :todo:
        :py:func:`ULA3.ancillary.aerosol.get_aerosol_value_for_region` accepts a default value
        which is returned if all else fails. Should we do a similar thing here? This would make
        the interface more consistent and potentially be more useful to users.

    :todo:
        ``wavelength_range`` should be documented by someone who knows what it is.

    """
    ul_lat = ur_lat
    ul_lon = ll_lon
    lr_lat = ll_lat
    lr_lon = ur_lon

    def findFile(brdf_root, brdf_dir, band_number, factor):
        '''
        Find the BRDF file corresponding to the wavelength of the band
        N.B: band_index argument is the 0-based band index. It is NOT the 1-based band number
        '''
        def find_wavelength_file(band_number, file_list, max_bandwidth_ratio = 5.0):
            '''Determine wavelength range for specified band index
            and find the first file in the list whose average wavelength
            falls into this range.
            N.B: band_index argument is the 0-based band index. It is NOT the 1-based band number
            '''

            fnames     = []
            brdf_upper = []
            brdf_lower = []

            for fname in sorted(file_list):
                s = re.search('.+_(\d+)_(\d+)nm(\w+)f.+\.hdf.*', fname)
                data_wavelength_range = [x/1000.0 for x in (int(s.group(1)), int(s.group(2)))]

                fnames.append(fname)
                brdf_upper.append(data_wavelength_range[1])
                brdf_lower.append(data_wavelength_range[0])

            # The following method is similar to residual analysis and is only a
            # temporary fix. The previous automated method didn't work for Landsat 5,7,8.
            # Temporary until Alex I. gets more time or another developer comes on board to
            # define hardcoded lookup tables within the satellite.xml configuration file.
            spectral_resids = []
            for i in range(len(brdf_upper)):
                r1 = abs(wavelength_range[0] - brdf_lower[i])
                r2 = abs(wavelength_range[1] - brdf_upper[i])
                spectral_resids.append(abs(r1 + r2))

            spectral_resids  = numpy.array(spectral_resids)
            brdf_match_loc   = numpy.argmin(spectral_resids)

            brdf_fname = fnames[brdf_match_loc]

            print
            print 'New BRDF Lookup!!!'
            print
            print 'band_number          ', band_number
            print 'wavelength_range     ', wavelength_range
            print 'brdf_lower           ', brdf_lower
            print 'brdf_upper           ', brdf_upper
            print 'fnames               ', fnames
            print 'spectral_resids      ', spectral_resids
            print 'brdf_match_loc       ', brdf_match_loc
            print 'brdf_fname           ', brdf_fname
            print

            return brdf_fname

        found_file = find_wavelength_file(band_number,
                         glob(os.path.join(brdf_root, brdf_dir, '*_*_*nm*f' + factor + '.hdf*')))
        assert found_file, 'No filename found for %s %s, %s, %s' % (brdf_root, brdf_dir, band_number, factor)
        return found_file # Should only have one match per directory


    initial_hdf_file = findFile(brdf_root, brdf_dir, int(band_string), factor)
    assert initial_hdf_file, 'BRDF file not found for %s, %s, %s, %s' % (brdf_root, brdf_dir, band_string, factor)

    if initial_hdf_file.endswith(".hdf.gz"):
        final_hdf_file = os.path.join(
            work_path,
            re.sub(".hdf.gz", ".hdf", os.path.basename(initial_hdf_file)))
        gunzipCmd = "gunzip -c %s > %s" % (initial_hdf_file, final_hdf_file)
        (status, msg) = commands.getstatusoutput(gunzipCmd)
        assert status == 0, "gunzip failed: %s" % msg
    else:
        final_hdf_file = initial_hdf_file

    # Use maximal axis-aligned extents for BRDF mean value calculation.
    # Note that latitude min-max logic is valid for the Southern
    # hemisphere only.

    nw = ( min(ul_lon, ll_lon),
           max(ul_lat, ur_lat), )

    se = ( max(lr_lon, ur_lon),
           min(ll_lat, lr_lat), )

    brdf_mean_value = BRDFLoader(final_hdf_file, UL=nw, LR=se).mean_data_value()

    # Remove temporary unzipped file
    if final_hdf_file.find(work_path) == 0:
        os.remove(final_hdf_file)

    return (final_hdf_file, brdf_mean_value)


# Maybe a candidate for constants.py
def brdf_wavelength_lut(satellite_sensor):
    """
    Retrieves the BRDF wavelengths for a given satellite-sensor.

    :param satellite_sensor:
        A string containing a valid satellite-sensor combination.
        Valid combinations are:
        landsat5tm
        landsat7etm
        landsat8oli
        landsat8olitirs

    :return:
        A dictionary containing the Band numbers of a sensor as the
        keys, and the BRDF wavelengths as the values.
    """

    input_str = str(satellite_sensor)

    BRDF_LUT = {
        'landsat5tm' : { 1 : '0459_0479nm',
                         2 : '0545_0565nm',
                         3 : '0620_0670nm',
                         4 : '0841_0876nm',
                         5 : '1628_1652nm',
                         7 : '2105_2155nm'
                       },
        'landsat7etm' : { 1 : '0459_0479nm',
                          2 : '0545_0565nm',
                          3 : '0620_0670nm',
                          4 : '0841_0876nm',
                          5 : '1628_1652nm',
                          7 : '2105_2155nm'
                        },
        'landsat8oli' : { 1 : '0459_0479nm',
                          2 : '0459_0479nm',
                          3 : '0545_0565nm',
                          4 : '0620_0670nm',
                          5 : '0841_0876nm',
                          6 : '1628_1652nm',
                          7 : '2105_2155nm'
                        },
        'landsat8olitirs' : { 1 : '0459_0479nm',
                               2 : '0459_0479nm',
                               3 : '0545_0565nm',
                               4 : '0620_0670nm',
                               5 : '0841_0876nm',
                               6 : '1628_1652nm',
                               7 : '2105_2155nm'
                            }
               }.get(input_str, 'Error')
