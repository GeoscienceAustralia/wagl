"""
Calculations for satellite positions and various associated angles. This module needs some reqorking. It has
evolved through several refactorings and has become difficult to follow. Some of the function signatures are
out of date (i.e. they receive arguments that are either unused or could be simplified for the sake of clarity).
The docstrings are also out of date and consequently the documentation is currently incomplete or misleading.

:todo:
    Cleanup, rework, document and write tests for this module.
"""
import math, os, re, logging, ephem, numpy, copy
import xml.dom.minidom
import scipy.interpolate
from datetime import timedelta
from . import dec_datetime
from vincenty import vinc_dist, vinc_pt
from ULA3.utils import unicode_to_ascii, log_multiline

logger = logging.getLogger('root.' + __name__)

class earth(object):

    # Mean radius
    RADIUS = 6371009.0  # (metres)

    # WGS-84
    #RADIUS = 6378135.0  # equatorial (metres)
    #RADIUS = 6356752.0  # polar (metres)

    # Length of Earth ellipsoid semi-major axis (metres)
    SEMI_MAJOR_AXIS = 6378137.0

    # WGS-84
    A = 6378137.0           # equatorial radius (metres)
    B = 6356752.3142        # polar radius (metres)
    F = (A - B) / A         # flattening
    ECC2 = 1.0 - B**2/A**2  # squared eccentricity

    MEAN_RADIUS = (A*2 + B) / 3

    # Earth ellipsoid eccentricity (dimensionless)
    #ECCENTRICITY = 0.00669438
    #ECC2 = math.pow(ECCENTRICITY, 2)

    # Earth rotational angular velocity (radians/sec)
    OMEGA = 0.000072722052





def N(lat, sm_axis=earth.SEMI_MAJOR_AXIS, ecc2=earth.ECC2):
    lat_s = math.sin(lat)
    return sm_axis / math.sqrt(1.0 - ecc2 * lat_s * lat_s)





def geocentric_lat(lat_gd, sm_axis=earth.A, ecc2=earth.ECC2):
    """Convert geodetic latitude to geocentric latitude.

    Arguments:
        lat_gd: geodetic latitude (radians)
        sm_axis: Earth semi-major axis (metres)
        ecc2: Earth squared eccentricity (dimensionless)

    Returns:
        Geocentric latitude (radians)
    """

    #return math.atan((1.0 - ecc2) * math.tan(lat_gd))
    return math.atan2(math.tan(lat_gd), 1.0/(1.0 - ecc2))


def geodetic_lat(lat_gc, sm_axis=earth.A, ecc2=earth.ECC2):
    """Calculate geocentric latitude to geodetic latitude.

    Arguments:
        lat_gc: geocentric latitude (radians)
        sm_axis: Earth semi-major axis (metres)
        ecc2: Earth squared eccentricity (dimensionless)

    Returns:
        Geodetic latitude (radians)
    """

    return math.atan2(math.tan(lat_gc), (1.0 - ecc2))




# These were found in 'earth.py', but appeared not be used
# def geocentric_lat(lat, sm_axis=SEMI_MAJOR_AXIS, ecc2=ECC2):
#     """Calculate geocentric latitude on the earth ellipsoid surface.
#
#     Arguments:
#         lat: geodetic latitude (radians)
#         sm_axis: length of earth ellipsoid semi-major axis (metres)
#         ecc2: squared eccentricity (dimensionless)
#
#     Returns:
#         Geocentric latitude value (radians).
#     """
#
#     lat_c = math.cos(lat)
#     lat_s = math.sin(lat)
#     N = sm_axis / math.sqrt(1.0 - ecc2 * lat_s * lat_s)
#     return lat - math.asin(N * ecc2 * lat_s * lat_c / sm_axis)
#
# def lat_geocentric(gdlat, sm_axis=SEMI_MAJOR_AXIS, ecc2=ECC2):
#     return math.atan((1.0 - ecc2) * math.tan(gdlat))
#
#
# def lat_geodetic(gclat, sm_axis=SEMI_MAJOR_AXIS, ecc2=ECC2):
#     return math.atan(math.tan(gclat) / (1.0 - ecc2))






class Satellite(object):
    """
    Manages data such as orbital parameters or sensor configuration for a satellite. Information on
    individual sensor bands and their respective wavelength ranges is read in for the specified
    sensor so that wavelength-dependent operations such as the acquisition of BRDF ancillary data can be
    performed dynamically rather than having to be hard-coded for specific band numbers.
    """
    _dom_tree = None

    def __init__(self, sat_name, sensor):
        """
        Constructor. Reads information on individual sensor bands and their respective wavelength ranges for
        the specified sensor from satellite.xml.

        :param sat_name:
            The namve of the satellite.

        :param sensor:
            The name of the sensor.

        """
        def parse_list(list_string, element_type=str):
            """Parses a string representation of a flat list into a list
            """
            return [element_type(element.strip()) for element in re.search('\s*(?<=\[)(.*)(?=\])\s*', list_string).group(0).split(',')]

        def get_sat_attributes(sat_node):
            """sets satellite attributes from XML file
            """
            self.NAME = unicode_to_ascii(sat_node.getAttribute('NAME'))
            self.TAG = unicode_to_ascii(sat_node.getAttribute('TAG'))
            assert self.TAG in ['LS5', 'LS7', 'LS8'], 'Unhandled satellite tag: ' + repr(self.TAG)

            # Orbital semi-major axis (metres)
            self.SEMI_MAJOR_AXIS = float(unicode_to_ascii(sat_node.getAttribute('SEMI_MAJOR_AXIS')))

            # Orbital radius (metres)
            self.RADIUS = float(unicode_to_ascii(sat_node.getAttribute('RADIUS')))

            # Orbital altitude (metres)
            self.ALTITUDE = float(unicode_to_ascii(sat_node.getAttribute('ALTITUDE')))

            # Orbital inclination (radians)
            self.INCLINATION = float(unicode_to_ascii(sat_node.getAttribute('INCLINATION')))
            self.INCL_SIN = math.sin(self.INCLINATION)
            self.INCL_COS = math.cos(self.INCLINATION)
            self.INCL_TAN = math.tan(self.INCLINATION)

            # Orbital angular velocity (radians/sec)
            self.OMEGA = float(unicode_to_ascii(sat_node.getAttribute('OMEGA')))

            # Sensor sweep period (sec)
            self.SWEEP_PERIOD = float(unicode_to_ascii(sat_node.getAttribute('SWEEP_PERIOD')))

            # TLE Format
            self.TLE_FORMAT = unicode_to_ascii(sat_node.getAttribute('TLE_FORMAT'))

            # Solar irradiance file
            self.SOLAR_IRRAD_FILE = unicode_to_ascii(sat_node.getAttribute('SOLAR_IRRAD_FILE'))

            # Spectral filter file
            self.SPECTRAL_FILTER_FILE = unicode_to_ascii(sat_node.getAttribute('SPECTRAL_FILTER_FILE'))

            # POOMA number - approximation only
            self.NOMINAL_PIXEL_DEGREES = float(unicode_to_ascii(sat_node.getAttribute('NOMINAL_PIXEL_DEGREES')))

            # Satellite number
            self.NUMBER = unicode_to_ascii(sat_node.getAttribute('NUMBER'))

            # Satellite classification
            self.CLASSIFICATION = unicode_to_ascii(sat_node.getAttribute('CLASSIFICATION'))

            # International designator
            # TLE compatible format: <last two digits of launch year><launch number of the year><piece of the launch>
            self.INTL_DESIGNATOR = unicode_to_ascii(sat_node.getAttribute('INTL_DESIGNATOR'))


        def get_sensor(sat_node, sensor):

            def get_bands(sensor_node):

                for band_node in sensor_node.getElementsByTagName('BAND'):
                    band_dict = {}
                    band_number = int(unicode_to_ascii(band_node.getAttribute('NUMBER')))
                    band_dict['NUMBER'] = band_number
                    self.BAND_TYPES['ALL'].append(band_number)

                    band_dict['NAME'] = unicode_to_ascii(band_node.getAttribute('NAME'))

                    # Store wavelength range as a two-element list
                    band_dict['WAVELENGTH'] = [float(wavelength) for wavelength in unicode_to_ascii(band_node.getAttribute('WAVELENGTH')).split('-')]

                    band_type = unicode_to_ascii(band_node.getAttribute('TYPE')).upper()
                    band_dict['TYPE'] = band_type
                    logger.debug('Band %d is of type %s', band_number, band_type)
                    # Store band in appropriate list according to type (Default to reflective)
                    band_list = self.BAND_TYPES.get(band_type)
                    if band_list is None:
                        band_list = self.BAND_TYPES['REFLECTIVE']
                    band_list.append(band_number)

                    band_dict['RESOLUTION'] = float(unicode_to_ascii(band_node.getAttribute('RESOLUTION')))

                    self.BAND_LIST.append(band_dict)

            self.sensor = None
            for sensor_node in sat_node.getElementsByTagName('SENSOR'):
                sensor_name = unicode_to_ascii(sensor_node.getAttribute('NAME'))
                #if sensor_name.lower() == sensor.lower():
                #if re.sub('\W+', '', sensor_name.lower()) == re.sub('\W+', '', sensor.lower()):
                if re.sub('[+,_,-]', '', sensor_name.lower()) == re.sub('[+,_,-]', '', sensor.lower()):
                    self.sensor = sensor_name

                    self.k = (float(unicode_to_ascii(sensor_node.getAttribute('K1'))),
                              float(unicode_to_ascii(sensor_node.getAttribute('K2'))))
                    break

            assert self.sensor, 'Sensor %s not found in configuration' % sensor

            self.sensor_description = unicode_to_ascii(sensor_node.getAttribute('DESCRIPTION'))
            self.root_band = int(unicode_to_ascii(sensor_node.getAttribute('ROOT_BAND')))
            self.rgb_bands = [int(bandstring) for bandstring in unicode_to_ascii(sensor_node.getAttribute('RGB_BANDS')).split(',')]
            self.acquistion_seconds = timedelta(seconds=float(unicode_to_ascii(sensor_node.getAttribute('ACQUISITION_SECONDS'))))

            get_bands(sensor_node)

        self.sensor = None
        self.sensor_description = None

        self.BAND_LIST = []

        # Define dict to contain all band groups
        self.BAND_TYPES = {}
        self.BAND_TYPES['ALL'] = []
        self.BAND_TYPES['REFLECTIVE'] = []
        self.BAND_TYPES['THERMAL'] = []
        self.BAND_TYPES['PANCHROMATIC'] = []
        self.BAND_TYPES['ATMOSPHERE'] = []
        self.rgb_bands = []

        self._tle_path_dict = {} # Dict to hold any TLE paths looked up by date

        if not Satellite._dom_tree: # If the satellite.xml file hasn't already been opened
            self.config_file = os.path.join(os.path.dirname(__file__), 'satellite.xml')

            logger.debug('Parsing XML file %s', self.config_file)

            # Open XML document using minidom parser
            Satellite._dom_tree = xml.dom.minidom.parse(self.config_file)

        self.NAME_PATTERN = None
        for sat_node in Satellite._dom_tree.getElementsByTagName('SATELLITE'):
            name_pattern = unicode_to_ascii(sat_node.getAttribute('NAME_PATTERN'))
            if re.match(re.compile(name_pattern, re.IGNORECASE), sat_name):
                self.NAME_PATTERN = name_pattern
                break

        assert self.NAME_PATTERN, 'Configuration for ' + sat_name + ' does not exist in ' + self.config_file

        get_sat_attributes(sat_node)
        get_sensor(sat_node, sensor)

#        log_multiline(logger.info, self.__dict__, 'Satellite object for ' + sat_name, '\t')


    def load_tle(self, centre_datetime, data_root, date_radius=45):
        """
        Loads satellite TLE (two-line element) for the given date and time.

        Arguments:
            centre_datetime: datetime.datetime instance
            data_root: root directory for TLE archive files
            date_radius: date radius for TLE search (days)

        Returns:
            ephem EarthSatellite instance

        """

        return (self._load_tle_from_archive(centre_datetime, data_root, date_radius) or
                self._load_tle_from_file_sequence(centre_datetime, data_root, date_radius))


    def _load_tle_from_archive(self, centre_datetime, data_root=None, date_radius=45):
        """Loads TLE (two-line element) for the satellite, date and time from
        an archive file.

        Arguments:
            centre_datetime: datetime.datetime instance
            data_root: directory containing TLE archive files
            date_radius: date radius for TLE search (days)

        Returns:
            ephem EarthSatellite instance

        """

        if data_root is None:
            logger.warning('load_tle_from_archive: data_root=None')
            return None

        TLE_ENTRY_PATTERN_FORMAT = (
            r'^([1])(\s+)'
            r'([%(NUMBER)s%(CLASSIFICATION)s]+)(\s+)([%(INTL_DESIGNATOR)s]+)(\s+)%(YYDDD)s'  # dict subst keys
            r'(\.)(\d+)(\s+)([\s\-]+)(\.)(\d+)(\s+)(\d+)([\-\+])(\d+)(\s+)(\d+)([\-\+])(\d+)(\s+)(\d+)(\s+)(\d+)(\s)'
            r'^([2])(.+)$'
        )

        def _cmp(x, y): return abs(x) - abs(y)  # For sorting date offsets.

        date_deltas = [timedelta(days=d) for d in sorted(range(-date_radius, date_radius), cmp=_cmp)]
        yyddd_list = [x.strftime('%02y%03j') for x in [(centre_datetime + d) for d in date_deltas]]

        tle_archive_path = os.path.join(
                               data_root,
                               self.NAME.replace('-', '').upper(),
                               'TLE', '%s_ARCHIVE.txt' % self.TAG
                           )

        tle_archive_text = ''

        with open(tle_archive_path, 'r') as fd:
            tle_archive_text = fd.read()

        logger.info('loaded TLE archive file: %s' % tle_archive_path)

        _smap = {
            'NUMBER': self.NUMBER,
            'CLASSIFICATION': self.CLASSIFICATION,
            'INTL_DESIGNATOR': self.INTL_DESIGNATOR,
        }

        tle_entry = None
        for yyddd in yyddd_list:
            _smap['YYDDD'] = yyddd
            m = re.search(TLE_ENTRY_PATTERN_FORMAT % _smap, tle_archive_text, re.MULTILINE)
            if m:
                logger.info('loaded TLE archive entry: %s' % centre_datetime)
                # Reconstitute TLE entry from regex match groups.
                tle_text = ''.join(m.groups()[0:6]) + yyddd + ''.join(m.groups()[6:])
                logger.info('TLE TEXT:\n%s' % tle_text)
                tle_lines = tle_text.split('\n')
                return ephem.readtle(self.NAME, tle_lines[0], tle_lines[1])

        logger.error('No TLE found for %s in archive file %s' % (centre_datetime, tle_archive_path))
        return None


    def _load_tle_from_file_sequence(self, centre_datetime, data_root=None, tle_search_range=45):
        """Load a TLE file for the specified datetime.

        Arguments:
            centre_datetime: scene centre datetime (datetime instance)
            data_root: ephemeris data root directory

        Returns:
            ephem EarthSatellite instance
        """

        if data_root is None:
            logger.warning('load_tle_from_sequence: data_root=None')
            return None

        def open_tle(tle_path, centre_datetime):
            """Function to open specified TLE file
            """
            try:
                fd = open(tle_path, 'r')
                tle_text = fd.readlines()
                logger.info('TLE file %s opened', tle_path)

                log_multiline(logger.debug, tle_text, 'TLE FILE CONTENTS', '\t')

                if self.TAG == 'LS5':
                    tle1, tle2 = tle_text[7:9]
                elif self.TAG == 'LS7':
                    tle1, tle2 = tle_text[1:3]

                sat_obj = ephem.readtle(self.NAME, tle1, tle2)

                # Cache TLE filename for specified date
                self._tle_path_dict[centre_datetime.date()] = tle_path

                return sat_obj
            finally:
                fd.close()

        # Check whether TLE path has been opened previously for this date
        tle_path = self._tle_path_dict.get(centre_datetime.date())
        if tle_path:
            return open_tle(tle_path, centre_datetime)

        data_subdir = re.sub('\W', '', self.NAME).upper()

        # Primary TLE: matches scene date.
        scene_doy = centre_datetime.strftime('%j')  # Note format: '%03d'

        tle_dir = os.path.join(
            data_root,
            data_subdir,  # 'LANDSAT5' or 'LANDSAT7'
            'TLE',
            '%s_YEAR' % self.TAG,
            '%4d' % centre_datetime.year
            )

        tle_file = self.TLE_FORMAT % (centre_datetime.year, scene_doy)
        tle_path = os.path.join(tle_dir, tle_file)

        if os.path.exists(tle_path):
            try:
                logger.debug('Opening primary TLE file %s', tle_path)
                return open_tle(tle_path, centre_datetime)
            except (Exception), e: #TODO: Tighten up this except clause - too general
                logger.warning('Unable to open primary TLE file %s: %s', tle_path, e.message)

        # Secondary TLE: closest TLE within specified number of days from scene date.
        for d in xrange(1, tle_search_range):
            ddelta = timedelta(days=d)
            for s in (-1, 1):
                dt = centre_datetime + (ddelta * s)
                tle_dir = os.path.join(data_root, data_subdir,
                    'TLE',
                    '%s_YEAR' % self.TAG,
                    '%4d' % dt.year)
                tle_file = self.TLE_FORMAT % (dt.year, dt.strftime('%j'))
                tle_path = os.path.join(tle_dir, tle_file)
                if os.path.exists(tle_path):
                    try:
                        logger.debug('Opening secondary TLE file %s', tle_path)
                        return open_tle(tle_path, centre_datetime)
                    except (Exception), e: #TODO: Tighten up this except clause - too general
                        logger.warning('Unable to open Secondary TLE file %s: %s', tle_path, e.message)

        logger.error('No valid TLE file found for %s in %s', dt.strftime('%Y-%m-%d'), tle_dir)
        return None





DEBUG_DISPLACEMENTS = False

PI = math.pi
TWO_PI = math.pi * 2
HALF_PI = math.pi / 2

# Number of out-of-bounds displacement passes for sampling
OUT_OF_BOUNDS_DISPLACEMENT_PASSES = 3

def eval_centre_data(scene_dataset, sat_ephem):
    """Sets up data structure for scene centre description"""
    doy = scene_dataset.scene_centre_datetime.timetuple().tm_yday
    
    start_datetime = scene_dataset.scene_centre_datetime - timedelta(0, 12)
    end_datetime = scene_dataset.scene_centre_datetime + timedelta(0, 12)

    centre_data = {
        'datetime': scene_dataset.scene_centre_datetime,
        'day_of_year': doy,
        'decimal_hour': dec_datetime.decimal_hour(scene_dataset.scene_centre_datetime),
        'decimal_day': dec_datetime.day_fraction(scene_dataset.scene_centre_datetime) + doy,
        'lon': float(scene_dataset.lonlats['CENTRE'][0]),
        'lat': float(scene_dataset.lonlats['CENTRE'][1])
    }
        #replaced the following in the above dictionary
        #'lon': float(scene_dataset.scene_centre_long),
        #'lat': float(scene_dataset.scene_centre_lat)

    # CONFIG.debug HACK
    # Adjust scene centre to match EQR test input
    # t = 23:45:43

    #=======================================================================
    # if debug_eqr.ENABLE:
    #    centre_data.lon = debug_eqr.CENTRE_LON
    #    centre_data.lat = debug_eqr.CENTRE_LAT
    #    _dt = datetime.datetime(d[0], d[1], d[2], 23, 45, 43)
    #    _tj = int(_dt.strftime('%j'))
    #    centre_data.datetime = _dt
    #    centre_data.day_of_year = _tj
    #    centre_data.decimal_hour = dec_datetime.hour_fraction(_dt) + _dt.hour
    #    centre_data.decimal_day = dec_datetime.day_fraction(_dt) + _tj
    #    print
    #    print '********************************************'
    #    print 'ula.main.eval_centre_data'
    #    print 'HACKED CENTRE DATA (to match EQR test scene)'
    #    print '********************************************'
    #    pprint(centre_data)
    #=======================================================================

    # CONFIG.debug
    # Initialize satellite obj for t = (scene datetime + delta)
    #debug_dt = centre_data.datetime + datetime.timedelta(60)
    #sat_ephem.compute(debug_dt)
    #print
    #print '*** CONFIG.debug', debug_dt
    #print '*** CONFIG.debug', math.degrees(sat_ephem.sublong.real), math.degrees(sat_ephem.sublat.real)

    sat_ephem.compute(centre_data['datetime'])

    sublon = sat_ephem.sublong.real
    sublat = sat_ephem.sublat.real

    centre_data['eph'] = {
        'lon_deg': math.degrees(sublon),
        'lat_deg': math.degrees(sublat),
        'lon_offset': math.degrees(sublon) - centre_data['lon'],
        'lat_offset': math.degrees(sublat) - centre_data['lat']
    }

    return centre_data


def compute_time_samples(
    lon_array,
    lat_array,
    centre_data,
    l1t_input_dataset, sat_ephem,
    nc=256, delta_tc=0.25, delta_rc=2500.0,
    sample_pts_file='sample_points.csv'):
    """
    Calculate samples for time, view angle and azimuth.

    Samples are obtains by propagating points along constant time curves
    away from the centre line. Centred point-to-point azimuth differencing
    along the centre line is used to obtain the local orthogonal (time)
    direction.

    :param lon_array:
        Longitude array (radians).

    :param lat_array:
        Latitude array (radians).

    :param centre_data:
        Scene centre data dictionary.

    :param cxform_from_geo:
        Cooordinate transform object (from geographic).

    :param cxform_to_geo:
        Cooordinate transform object (to geographic)

    :param spatial_ref:
        Scene spatial reference object.

    :param meta_mtl:
        MTL metadata dictionary.

    :param nc:
        Number of centre line point samples in forward and reverse time.

    :param delta_tc:
        Time interval between centre line samples (decimal seconds).

    :param delta_rc:
        Propagation step size (metres).

    :param sample_points_file:
        Name of sample points CSV file

    :return:
        8-tuple:
            (
                x_samples,     # x coordinates
                y_samples,     # y coordinates
                t_samples,     # time values
                v_samples,     # view angle values
                a_samples,     # azimuth values
                dlat_samples,  # delta_latitude values
                gamma_samples, # gamma values
                slat_samples,  # satellite nadir latitude values
                cline_pts,     # centre line point coordinates
            )

    :todo:
        Roger notes that this function needs revisiting. Some of the arguments are pretty much redudant
        and the calculations, while stable, have evolved over time such the code is not easy to understand
        and could be simplifed.

    """
    def in_bounds(x, y, x_range, y_range, x_eps=0.0, y_eps=0.0):
        """
        Check that 2D coordinates are in bounds.

        :param x:
            X coordinate.

        :param y:
            Y coordinate.

        :param x_range:
            X coordinate range.
        :type x_range:
            2-tuple:
                (min, max)

        :param y_range:
            Y coordinate range
        :type y_range:
            2-tuple:
                (min, max)

        :param x_eps:
            X coordinate epsilon.

        :param y_eps:
            Y coordinate epsilon.

        :return:
            ``True`` if (x, y) is within the specified range, otherwise ``False``.
        """

    #        logger.debug('in_bounds(%s, %s, %s, %s, %s, %s) called', x, y, x_range, y_range, x_eps, y_eps)
        return ( (x_range[0] - x_eps) <= x <= (x_range[1] + x_eps) and
                 (y_range[0] - y_eps) <= y <= (y_range[1] + y_eps) )


    logger.debug('compute_time_samples(%s, %s, %s, %s, %s, %s, %s, %s) called',
                 lon_array, lat_array, centre_data, l1t_input_dataset,
                 nc, delta_tc, delta_rc, sample_pts_file)

    assert lon_array.shape == lat_array.shape, 'ERROR: LON, LAT SHAPE MISMATCH'
    _shape = lon_array.shape

    satellite = l1t_input_dataset.satellite

    HR_RATIO = satellite.ALTITUDE / earth.MEAN_RADIUS
    HR_PLUS_ONE = HR_RATIO + 1

    # X, Y extents from MTL metadata.
    # Use l1t_input_dataset instance values
    x_coord_range, y_coord_range = l1t_input_dataset.get_bounds()
    logger.debug('x_coord_range = %s', x_coord_range)
    logger.debug('y_coord_range = %s', y_coord_range)

    # Calculate ephem centre lon, lat offsets.

    centre_time = centre_data['datetime']

    # WARNING The first sat_ephem.compute() call must use a time argument of the
    # form (y, m, decimal_day), as in the centre line generation loop below.
    # If other arg types are used the Earth object cannot resolve sub-second
    # time intervals.

    sat_ephem.compute((centre_time.year, centre_time.month, dec_datetime.decimal_day_of_month(centre_time)))

    lon_offset = math.radians(centre_data['lon']) - sat_ephem.sublong.real
    lon_offset_deg = math.degrees(lon_offset)

    centre_gc_lat = geocentric_lat(math.radians(centre_data['lat']))

    gc_lat_offset = centre_gc_lat - sat_ephem.sublat.real
    gc_lat_offset_deg = math.degrees(gc_lat_offset)

#===========================================================================
#    if CONFIG.debug:
#        print
#        print '***** CONFIG.debug *****'
#        print 'centre_data'
#        ppr(centre_data)
#        print centre_time
#        print (centre_time.year, centre_time.month, dec_datetime.decimal_day_of_month(centre_time))
#        print 'MTL centre_data.lat (geodetic)', math.radians(centre_data.lat)
#        print 'MTL centre_data.lat (geocentric)', centre_gc_lat
#        print 'MTL centre E,N', cxform_from_geo.TransformPoint(centre_data.lon, centre_data.lat, 0)
#        print 'SAT sublong, sublat (geocentric)', sat_ephem.sublong.real, sat_ephem.sublat.real
#        print 'geocentric offset (rad)', lon_offset, gc_lat_offset
#        print 'geocentric offset (deg)', math.degrees(lon_offset), math.degrees(gc_lat_offset)
#        print 'checksum', centre_gc_lat - (sat_ephem.sublat.real + gc_lat_offset)
#        print 'lat return trip', geocentric['geodetic']_lat(sat_ephem.sublat.real + gc_lat_offset)
#        print
#
#    # Generate centre line points.
#
#    print
#    print 'CENTRE LINE'
#===========================================================================

    tdelta = timedelta(seconds=delta_tc)
    centre_pts = {}

    discontinuity = False
    max_timeu = 0 # Latest UTC PM time detected
    for i in xrange(-nc/2, nc/2):
        t = centre_time + tdelta * i
        timeu = dec_datetime.decimal_hour(t)
        
        # Need to keep checking for discontinuity while one hasn't been detected
        discontinuity = discontinuity or max_timeu - timeu > 12.0
        if discontinuity and timeu < 12.0: # If discontinuity detected and time is AM
            # Offset AM times by 24h to avoid discontinuity
            timeu += 24.0
        else:            
            max_timeu=max(max_timeu, timeu) # Only need to update this if no discontinuity detected yet
            
        # WARNING
        # ephem can't seem to resolve subsecond timedeltas unless a
        # time tuple of the form (year, month, decimal_day) is used.

        sat_ephem.compute((t.year, t.month, dec_datetime.decimal_day_of_month(t)))

        # Correct for ephem coordinate discrepancy.
        lon = sat_ephem.sublong.real + lon_offset
        glat = sat_ephem.sublat.real + gc_lat_offset
        lat = geodetic_lat(glat)

        x, y, _h = l1t_input_dataset.cxform_from_geo.TransformPoint(math.degrees(lon), math.degrees(lat), 0)

        rho = math.acos(math.sin(glat) / math.sin(satellite.INCL_SIN))
        beta = math.atan2(-1.0, satellite.INCL_TAN * math.sin(rho))
        cos_glat = math.cos(glat)
        skew = math.atan2( earth.OMEGA * cos_glat * math.cos(beta),
                           satellite.OMEGA + earth.OMEGA * cos_glat * math.sin(beta) )

        # Renormalize heading and skew to positive values (clockwise from North)
        # to match the displacement sign convention.
        hbeta = (beta if beta >= 0 else (beta + TWO_PI))
        skew = (skew if beta >= 0 else (-skew))

        logger.debug('%-27s %-13s (%12.8f, %12.8f) (%e, %e) %f %s %s',
                  str(t), timeu,
                  math.degrees(lon), math.degrees(lat), x, y, math.degrees(hbeta),
                  in_bounds(x, y, x_coord_range, y_coord_range),
                  ('CENTRE' if t == centre_time else '')
              )

        centre_pts[t] = {
            'in_bounds': in_bounds(x, y, x_coord_range, y_coord_range),
            'geocentric': (lon, glat),
            'geodetic': (lon, lat),
            'projected': (x, y),
            'az_forward': None, # NOT USED
            'az_tangent': None, # NOT USED
            'hbeta': hbeta,
            'skew': skew,
        }

    times = sorted(centre_pts)

    assert any([p['in_bounds'] for p in centre_pts.itervalues()]), 'No points in bounds'
        #===================================================================
        # print
        # print 'ERROR: CENTRE POINT GENERATION FAILED (NO POINTS IN-BOUNDS)'
        # print 'x_coord_range', x_coord_range
        # print 'y_coord_range', y_coord_range
        # sys.exit('CALCULATION ABORTED')
        #===================================================================

    ###########################################################################
    # OLD Forward/centred centre line azimuth differencing.
    # Azimuths obtained by this method are skewed relative to the satellite
    # heading.
    ###########################################################################
    # The azimuth of the current-to-next or prev-to-next vectors approximates
    # the azimuth of the tangent vector to the centre line at each point.
    # Normal directions to the centre line are ~ (tangent_azimuth +- pi/2).
    #
    # Centre line direction via forward azimuth differencing (p(i) -> p(i+1))
    # OPT az_forward not used
    #for t in times[0:-1]:
    #    p = centre_pts[t         ]['geodetic']
    #    q = centre_pts[t + tdelta]['geodetic']
    #    cdist, cfaz, craz = vinc_dist(earth.F, earth.A, p[1], p[0], q[1], q[0])
    #    #print 'FDIFF', t, cfaz
    #    centre_pts[t]['az_forward'] = cfaz
    #
    # Centre line direction via centred azimuth differencing (p(i-1) -> p(i+1))
    #
    #for t in times[1:-1]:
    #    p = centre_pts[t - tdelta]['geodetic']
    #    q = centre_pts[t + tdelta]['geodetic']
    #    cdist, cfaz, craz = vinc_dist(earth.F, earth.A, p[1], p[0], q[1], q[0])
    #    #print 'CDIFF', t, cfaz
    #    centre_pts[t]['az_tangent'] = cfaz
    #
    # CONFIG.debug
    #print
    #for t in times:
    #    if centre_pts[t]['in_bounds']:
    #        azfd = centre_pts[t]['az_forward']
    #        azcd = centre_pts[t]['az_tangent']
    #        azhb = centre_pts[t]['hbeta']
    #        skew = centre_pts[t]['skew']
    #        print t, azfd, azcd, skew, azhb, (azhb + skew - azcd)
    ###########################################################################

    # Displace centre line points along curves of constant time.
    # Displacement direction is ~ (heading +- pi/2).

    n_samples_buf = 32000

    x_samples = numpy.zeros((n_samples_buf,))
    y_samples = numpy.zeros_like(x_samples)
    t_samples = numpy.zeros_like(x_samples)
    v_samples = numpy.zeros_like(x_samples)
    a_samples = numpy.zeros_like(x_samples)
    slat_samples = numpy.zeros_like(x_samples)
    dlat_samples = numpy.zeros_like(x_samples)
    gamma_samples = numpy.zeros_like(x_samples)

    sample_pts_list = []
    sample_pt_index = 0

    # Initialise container for L/R azimuth samples.

    x_lr_samples = numpy.zeros((n_samples_buf/2,))

    az_samples = {
                     'L': {
                         'count': 0,
                         'ipass': [],
                         'x': x_lr_samples,
                         'y': numpy.zeros_like(x_lr_samples),
                         'a': numpy.zeros_like(x_lr_samples),
                     },
                     'R': {
                         'count': 0,
                         'ipass': [],
                         'x': numpy.zeros_like(x_lr_samples),
                         'y': numpy.zeros_like(x_lr_samples),
                         'a': numpy.zeros_like(x_lr_samples),
                     },
                 }

    vu_samples = {
                     'L': {
                         'count': 0,
                         'ipass': [],
                         'x': x_lr_samples,
                         'y': numpy.zeros_like(x_lr_samples),
                         'v': numpy.zeros_like(x_lr_samples),
                     },
                     'R': {
                         'count': 0,
                         'ipass': [],
                         'x': numpy.zeros_like(x_lr_samples),
                         'y': numpy.zeros_like(x_lr_samples),
                         'v': numpy.zeros_like(x_lr_samples),
                     },
                 }

    # Main displacement loop.

    disp_map = {'L': HALF_PI, 'R': -HALF_PI}

    for disp_key, disp_az_delta in disp_map.iteritems():

        t_list = copy.deepcopy(times)
        t_used = []
        t_used.append(t_list.pop(0))
        t_used.append(t_list.pop())

        disp_pts = {}
        disp_count = 1

        # CONFIG.debug
        # Re-initialize sample array for each side of the centre line.
        #sample_pts_list = []
        #sample_pt_index = 0

        pts = copy.deepcopy(centre_pts)

        # Propagate sample points away from the centre line in the given
        # direction.

        out_pass_count = 0
        while len(t_used) < len(times):

            if DEBUG_DISPLACEMENTS:
                logger.debug('DISPLACEMENT PASS: %s %s', disp_key, disp_count)

            # Store the starting index for each pass.

            az_samples[disp_key]['ipass'].append( az_samples[disp_key]['count'] )
            vu_samples[disp_key]['ipass'].append( vu_samples[disp_key]['count'] )

            discontinuity = False
            max_timeu = 0
            for t in sorted(t_list):

                disp_direction = pts[t]['hbeta'] + disp_az_delta
                disp_distance = (5.0 if disp_count == 1 else delta_rc)

                lat, lon, _raz = vinc_pt(earth.F, earth.A,
                                        pts[t]['geodetic'][1], pts[t]['geodetic'][0],
                                        disp_direction, disp_distance)

                x, y, _h = l1t_input_dataset.cxform_from_geo.TransformPoint(math.degrees(lon), math.degrees(lat), 0)
                flag = in_bounds(x, y, x_coord_range, y_coord_range)
                clon, clat = centre_pts[t]['geodetic']

                (centre_arc_distance,
                 az_to_centre,
                 _az_from_centre) = vinc_dist(earth.F, earth.A, lat, lon, clat, clon)

                latgc = geocentric_lat(lat)
                clatgc = geocentric_lat(clat)

                dlat_geodetic = lat - clat
                dlat = latgc - clatgc
                dlon = lon - clon

                # Central angle: spherical cosine formula
                # Note: Geocentric latitude here, not spc colatitude.
                q = (math.sin(latgc) * math.sin(clatgc) +
                     math.cos(latgc) * math.cos(clatgc) * math.cos(dlon))
                gamma_sc = math.acos(q)

                # Central angle: haversine formula
                q = (math.sin(dlat/2)**2 +
                     math.cos(latgc) * math.cos(clatgc) * math.sin(dlon/2)**2)
                gamma_ha = math.asin(math.sqrt(q)) * 2

                # Central angle: radians = arc_length / radius
                gamma_arc = centre_arc_distance / earth.MEAN_RADIUS

                # Choose a gamma...
                gamma = gamma_ha

                # Satellite azimuth: The angle between the local meridian and
                # point-to-centre arc (beta, Earth angle) is the complement
                # of the azimuth. Counterclockwise rotation => positive beta.
                # NOT REQUIRED *** USE VINCENTY AZIMUTH INSTEAD OF (PI - BETA)
                #beta = math.copysign(math.acos(math.tan(dlat) / math.tan(gamma)), disp_az_delta)

                # Satellite view angle (signed)
                # Interpolating unsigned sample values with scipy.interpolate.griddata
                # produces a result without any v=0 pixels.
                alpha = math.atan2(math.sin(gamma), HR_PLUS_ONE - math.cos(gamma))
                alpha = math.copysign(alpha, disp_az_delta)

                # Time

                timeu = dec_datetime.decimal_hour(t)
        
                # Need to keep checking for discontinuity while one hasn't been detected
                discontinuity = discontinuity or max_timeu - timeu > 12.0
                if discontinuity and timeu < 12.0: # If discontinuity detected and time is AM
                    # Offset AM times by 24h to avoid discontinuity
                    timeu += 24.0
                else:            
                    max_timeu=max(max_timeu, timeu) # Only need to update this if no discontinuity detected yet

                # Satellite heading
                # (beta, from www.eoc.csiro.au/hswww/oz_pi/util/heading.pdf)

                rho = math.acos(math.sin(latgc) / math.sin(satellite.INCL_SIN))
                beta = math.atan2(-1.0, satellite.INCL_TAN * math.sin(rho))
                cos_latgc = math.cos(latgc)
                skew = math.atan2( earth.OMEGA * cos_latgc * math.cos(beta),
                                   satellite.OMEGA + earth.OMEGA * cos_latgc * math.sin(beta) )

                # Renormalize values to match the displacement sign convention.
                hbeta = (beta if beta >= 0 else (beta + TWO_PI))
                skew = (skew if beta >= 0 else (-skew))

                disp_pts[t] = {
                                  'in_bounds': flag,
                                  'geodetic': (lon, lat),
                                  'geocentric': (lon, latgc),
                                  'projected': (x, y),
                                  'centre_arc_distance': centre_arc_distance,
                                  'gamma': gamma,
                                  'alpha': alpha,
                                  'az_forward': None, # NOT USED
                                  'az_tangent': None, # NOT USED
                                  'hbeta': hbeta
                              }

                x_samples[sample_pt_index] = x
                y_samples[sample_pt_index] = y
                t_samples[sample_pt_index] = timeu
                v_samples[sample_pt_index] = alpha
                a_samples[sample_pt_index] = az_to_centre
                dlat_samples[sample_pt_index] = math.copysign(dlat, disp_az_delta)
                slat_samples[sample_pt_index] = clatgc
                gamma_samples[sample_pt_index] = gamma

                sample_pt_index += 1
                assert sample_pt_index < n_samples_buf, 'SAMPLE BUFFER OVERFLOW'

                # Segregate azimuth samples for L/R interpolation.

                iaz = az_samples[disp_key]['count']
                az_samples[disp_key]['x'][iaz] = x
                az_samples[disp_key]['y'][iaz] = y
                az_samples[disp_key]['a'][iaz] = az_to_centre
                az_samples[disp_key]['count'] += 1

                # Segregate view angle samples for L/R interpolation.

                ivu = vu_samples[disp_key]['count']
                vu_samples[disp_key]['x'][ivu] = x
                vu_samples[disp_key]['y'][ivu] = y
                vu_samples[disp_key]['v'][ivu] = alpha
                vu_samples[disp_key]['count'] += 1

                s_pt = (
                    '%10.6f, %10.6f, %e, %e, %e, %e # '
                    'd=%8.6f, gs=%8.6f gh=%8.6f ga=%8.6f b=%8.6f az=%8.6f azv=%8.6f # %s %s' ) % (
                        math.degrees(lon), math.degrees(lat), x, y, timeu,
                        alpha, centre_arc_distance, gamma_sc, gamma_ha, gamma_arc,
                        beta, (math.pi - beta), az_to_centre, disp_key, flag
                    )

                if DEBUG_DISPLACEMENTS:
                    logger.debug('s_pt - %s', s_pt)

                sample_pts_list.append(s_pt)

            # Centre line tangent direction.

            t_used.append(t_list.pop(0))
            t_used.append(t_list.pop())

            ###########################################################################
            # OLD Forward/centred centre line azimuth differencing.
            # Azimuths obtained by this method are skewed relative to the satellite
            # heading.
            ###########################################################################
            #for t in sorted(t_list):
            #    # direction from previous point (t - tdelta) to next point (t + tdelta)
            #    q = disp_pts[t + tdelta]['geodetic']
            #    p = disp_pts[t - tdelta]['geodetic']
            #    cdist, cfaz, craz = vinc_dist(earth.F, earth.A, p[1], p[0], q[1], q[0])
            #    disp_pts[t].az_tangent = cfaz # NOT USED
            ###########################################################################

            pts.update(disp_pts)
            disp_count += 1

            # Stop propagation when all sample points go out of scene.
            # griddata will occasionally produce a sliver of NaN values at the
            # UL or LR scene corner if displacements stop immediately when all
            # points are out of scene. Going 2-3 passes beyond the boundary
            # ensures there are enough sample points for interpolation.

            out_of_bounds = all([(not p['in_bounds']) for p in disp_pts.itervalues()])

            if out_of_bounds:
                out_pass_count += 1
                if out_pass_count > OUT_OF_BOUNDS_DISPLACEMENT_PASSES:
                    logger.info('%s PROPAGATION STOPPED (count = %d)', disp_key, disp_count)
                    break

    # Add centre points to samples

    centre_line_spline_pts = []

    eps_y = ((y_coord_range[1] - y_coord_range[0]) / (_shape[0] - 1)) * 4

    for t in sorted(centre_pts):

        lon, lat = centre_pts[t]['geodetic']
        x, y = centre_pts[t]['projected']
        flag = in_bounds(x, y, x_coord_range, y_coord_range)
        timeu = dec_datetime.decimal_hour(t)
        alpha = 0.0

        x_samples[sample_pt_index] = x
        y_samples[sample_pt_index] = y
        t_samples[sample_pt_index] = timeu
        v_samples[sample_pt_index] = alpha

        # Azimuth is singular along the centre line.
        # Omit, set to zero, or set to numpy.nan?
        # Setting to numpy.nan leads to griddata/cubic and bispl failure.
        #a_samples[sample_pt_index] = numpy.nan
        #a_samples[sample_pt_index] = 0.0

        sample_pt_index += 1

        # Add CSV entry.

        s_pt = ( '%f, %f, %e, %e, %e, %e # '
                 'd=%f, gs=%f gh=%f ga=%f az=%f azv=%f # %s %s' ) % (
                  math.degrees(lon), math.degrees(lat), x, y, timeu, alpha,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  #centre_arc_distance, gamma_sc, gamma_ha, gamma_arc, (math.pi - beta), az_to_centre
                  'CENTRE', flag
               )

        if DEBUG_DISPLACEMENTS:
            logger.debug('s_pt = %s', s_pt)

        sample_pts_list.append(s_pt)

        if in_bounds(x, y, x_coord_range, y_coord_range, 0.0, eps_y):
            centre_line_spline_pts.append((timeu, x, y))

    # Write points to a file.

    if sample_pts_file:
        header_txt = [
            '# SAMPLE_POINTS',
            '# %s' % l1t_input_dataset.pathname,
            '# %s points' % len(sample_pts_list),
            '# lon, lat, x, y, tu, alpha, # centre_arc_distance gamma_sc gamma_ha gamma_arc beta az az_vincenty # disp_direction in_scene'
            '#\n',
        ]
        try:
            fd = open(sample_pts_file, 'w')
            fd.write('\n'.join(header_txt))
            fd.write('\n'.join(sample_pts_list))
        finally:
            fd.close()
        logger.info('CREATED DISPLACEMENT POINTS FILE: %s', sample_pts_file)

    logger.debug('TIME AND ANGLE SAMPLE POINTS (n=%d, buffer_size=%d)',
              len(sample_pts_list), n_samples_buf
          )
    logger.debug('PARAMETERS: nc = %d, delta_tc = %f, delta_rc = %f',
              nc, delta_tc, delta_rc
          )
    logger.debug('X     %s', x_samples)
    logger.debug('Y     %s', y_samples)
    logger.debug('T     %s', t_samples)
    logger.debug('V     %s', v_samples)
    logger.debug('A     %s', a_samples)
    logger.debug('DLAT  %s', dlat_samples)
    logger.debug('SLAT  %s', slat_samples)
    logger.debug('GAMMA %s', gamma_samples)

    # Slice unused (zero) elements from the sample value buffers.
    # Unused elements will bias the interpolation result.

    n = sample_pt_index
    assert n > 0, 'ERROR: NO SAMPLE POINTS'

    return (
        x_samples[:n],
        y_samples[:n],
        t_samples[:n],
        v_samples[:n],
        a_samples[:n],
        dlat_samples[:n],
        gamma_samples[:n],
        slat_samples[:n],
        centre_line_spline_pts,
        az_samples, vu_samples
    )


def create_centre_line_file(clpts_tuple, shape, x_range, y_range, cxform_to_geo, view_angle_max, cl_file_path='CENTRELINE_DEBUG'):
    """Creates centre line file by spline fit to sampled centre points.

    Arguments:
        clpts_tuple: list containing centre line coordinates,
                     [(<time>, <x>, <y>), ...].
        shape: scene shape (2-tuple: (<rows>, <ncols>))
        x_range: x-coordinate range (2-tuple: (<min>, <max>))
        y_range: y-coordinate range (2-tuple: (<min>, <max>))
        cxform_to_geo: coordinate transform (to geographic)
        cl_file_path: path of centre line file
    """
    logger.debug('create_centre_line_file(%s, %s, %s, %s, %s, %s) called',
                 clpts_tuple, shape, x_range, y_range,
                 cxform_to_geo, cl_file_path)

    # Unpack X and Y coordinates from tuple.
    # Time (element 0 of each tuple) is not currently used.
    xc = [_t[1] for _t in clpts_tuple]
    yc = [_t[2] for _t in clpts_tuple]

    logger.debug('xc = %s', xc)
    logger.debug('yc = %s', yc)

    # Interpolation coordinates for each row.
    ycinterp = numpy.linspace(y_range[0], y_range[1], shape[0])

    logger.debug('CENTRE LINE: X = f(Y) splrep/splev')
    logger.debug('Y range = %s', y_range)
    logger.debug('ycinterp %s %s', ycinterp.size, ycinterp)

    # splrep, splev require data values in ascending order.
    # Note: xc, yc are lists, not numpy arrays.
    if xc[0] > xc[1]:
        xc.reverse()
        yc.reverse()

    #
    # Interpolate the centre line sample points across the scene using
    # scipy.splrep/splev.
    #
    # Smoothing factors <= 1.0 hit the points using 2500m samples, but
    # occasionally produce curves with divergence at one or both ends.
    # For now, use a hardwired smoothing factor in the range 10-100.
    #
    # TODO adaptive smoothing factor?
    #

    logger.info('interpolating centre line samples using scipy.interpolate.splrep/splrev')

    smoothing_factor = 10.0

    logger.info('smoothing factor = %f', smoothing_factor)

    xcrep = scipy.interpolate.splrep(yc, xc, s=smoothing_factor)
    xcspl = scipy.interpolate.splev(ycinterp, xcrep)
    xdeltas = scipy.interpolate.splev(yc, xcrep) - xc
    mean_abs_xdelta = numpy.mean(numpy.abs(xdeltas))

    # Sanity check the terminal fit points against the sample points
    # nearest the top and bottom scene edges.

    delta_beg = xcspl[0] - xc[0]
    delta_end = xcspl[-1] - xc[-1]

    logger.info('endpoint deltas = %s, %s', delta_beg, delta_end)
    logger.debug('xcspl.size = %s', xcspl.size)
    logger.debug('xcspl = %s', xcspl)
    logger.info('xdeltas:\n%s', xdeltas)
    logger.info('mean_abs_xdelta = %s', mean_abs_xdelta)

    #
    # Fit a cubic curve to the centre line samples if the splrep/splev
    # result looks unstable.
    #
    # TODO improve stability metric?
    #

    tpoint_delta_tol = 1000.0

    if abs(delta_beg) > tpoint_delta_tol or abs(delta_end) > tpoint_delta_tol:

        kdeg = 3
        logger.info('interpolating centre line samples using numpy.polyfit (deg=%d)', kdeg)

        pfyx = numpy.polyfit(yc, xc, kdeg)
        pf_func = numpy.poly1d(pfyx)
        xcfit = pf_func(ycinterp)

        xdeltas = pf_func(yc) - xc
        mean_abs_xdelta = numpy.mean(numpy.abs(xdeltas))

        delta_beg = xcfit[0] - xc[0]
        delta_end = xcfit[-1] - xc[-1]

        logger.info('endpoint deltas = %s, %s', delta_beg, delta_end)
        logger.debug('xcfit.size = %s', xcfit.size)
        logger.debug('xcfit = %s', xcfit)
        logger.info('xdeltas:\n%s', xdeltas)
        logger.info('mean_abs_xdelta = %s', mean_abs_xdelta)

        xcspl = xcfit  # reference for downstream code

    # Reverse coordinate arrays so that regular iteration gives top-down
    # order, as needed in the centreline file.

    yyc = list(ycinterp)
    yyc.reverse()
    xxc = list(xcspl)
    xxc.reverse()

    def _index(q, q_range, n):
        xq =  (((q - q_range[0]) / (q_range[1] - q_range[0])) * n)
        return int(round(xq))

    cl_record_fmt = '%12d %12d %20.12f %20.12f'
    cl_text = [
        '%16.8f' % view_angle_max,
        '%12d %12d' % shape,
    ]

    for ic, yci in enumerate(yyc):
        xci = xxc[ic]
        lon, lat, _z = cxform_to_geo.TransformPoint(xci, yci, 0)

        # centreline record format:
        # '<irow (1-based)> <jcentre (1-based)> <lat (deg)> <lon (deg)>'

        rec = cl_record_fmt % (ic + 1, _index(xci, x_range, shape[1]), lat, lon)
        cl_text.append(rec)

    cl_text.append('\n')

    with open(cl_file_path, 'w') as _fd:
        _fd.write('\n'.join(cl_text))

    logger.info('wrote centreline file: %s', cl_file_path)

    # Generate centre line index array.

    cl_indices = numpy.zeros((shape[0],), dtype=numpy.int)

    for ic, yci in enumerate(yyc):
        cl_indices[ic] = _index(xxc[ic], x_range, shape[1])

    return cl_indices


def interpolate_sample_points(xy, f, mg, method='linear', fill_value=numpy.nan):
    """Interpolate sample points in 2D.

    Arguments:
        xy: sample coordinates (2-tuple: (x, y))
        f: function sample values (1D array)
        mg: interpolation meshgrid (result of numpy.meshgrid())
        method: interpolation method: 'cubic', 'linear' (default), or 'nearest'
        fill_value: fill value for points outside convex hull

    Returns:
        Interpolated grid.
    """
    logger.debug('interpolate_sample_points(%s, %s, %s, %s, %s) called', xy, f, mg, method, fill_value)

    assert xy[0].shape == xy[1].shape == f.shape, 'ERROR: NON-CONFORMING SAMPLE ARRAYS'

    # OLD
    # Interpolation grid coordinates
    #    - shape = (nrows[Y], ncols[X])
    #    - Y runs max-to-min, to give northing = up
    #mg = numpy.meshgrid( numpy.linspace(x_range[0], x_range[1], shape[1]),
    #                     numpy.linspace(y_range[1], y_range[0], shape[0]) )

    return scipy.interpolate.griddata(xy, f, mg, method, fill_value)


def eval_sat_grids(
    lon_array,
    lat_array,
    centre_data,
    l1t_input_dataset,
    sat_ephem,
    view_angle_max,
    work_path,
    debug):
    """
    Computes time, satellite view angle and satellite azimuth grids.

    :param lon_array:
        Grid of longitudes for each cell (radians).
    :type lon_array:
        :py:class:`numpy.ndarray` convertible to ???

    :param lat_array:
        Grid of latitudes for each pixel (radians).
    :type lat_array:
        :py:class:`numpy.ndarray` convertible to ???

    :param centre_data:
        Scene centre data dictionary.
    :type centre_data:
        ???

        cxform_from_geo: cooordinate transform object (from geographic)
        cxform_to_geo: cooordinate transform object (to geographic)
        spatial_ref: spatial reference object

    :return:
        3-tuple of :py:class:`ndarray`s (time, view_angle (radians), azimuth (radians))

    :todo:
        Better description of param ``centre_data`` is required.

    """
    logger.debug('eval_sat_grids(%s, %s, %s, %s) called', lon_array, lat_array, centre_data, l1t_input_dataset)

    dlat = slat = gamma = None

    xgr, ygr = l1t_input_dataset.get_bounds()

    # Generate sample values along constant time curves away from
    # the centre line.

    _x, _y, _t, _v, _a, _dlat, _gamma, _slat, cline_points, az_samples, vu_samples = compute_time_samples(
        lon_array, lat_array, centre_data, l1t_input_dataset, sat_ephem,
        sample_pts_file=os.path.join(work_path, 'sample_points.csv')
        )

    # Create centre line file and L/R mask array.

    clfpath = os.path.join(work_path, 'CENTRELINE')
    cl_indices = create_centre_line_file(cline_points, lon_array.shape, xgr, ygr, l1t_input_dataset.cxform_to_geo, view_angle_max, clfpath)

    cmask = numpy.zeros(lon_array.shape).astype(numpy.bool)
    for i in xrange(cmask.shape[0]):
        cmask[i, cl_indices[i]:] = True

    # Define interpolation grid to match projected scene extents.
    # Note: Y max-to-min produces northing up.
    nrows = lon_array.shape[0]
    ncols = lon_array.shape[1]

    xi = numpy.linspace(xgr[0], xgr[1], ncols)
    yi = numpy.linspace(ygr[1], ygr[0], nrows)  # Max-to-min

    xy = (_x, _y)
    mx, my = numpy.meshgrid(xi, yi)
    mg = (mx, my)

    # Cull out-of-scene sample points for bspl/spline interpolation.

    #in_scene = (_x >= xgr[0]) & (_x <= xgr[1]) & (_y >= ygr[0]) & (_y <= ygr[1])
    #
    #_xs = _x.compress(in_scene)
    #_ys = _y.compress(in_scene)
    #_ts = _t.compress(in_scene)
    #_vs = _v.compress(in_scene)
    #_dlats = _dlat.compress(in_scene)
    #_slats = _slat.compress(in_scene)
    #_gammas = _gamma.compress(in_scene)
    #_as = _a.compress(in_scene)
    #
    #xis = numpy.linspace(xgr[0], xgr[1], ncols)
    #yis = numpy.linspace(ygr[0], ygr[1], nrows)  # Min-to-max
    #
    #print 'sat_grids.eval: culled %d of %d out-of-scene sample points' % (_x.size - _xs.size, _x.size)

    # Interpolate L/R satellite azimuth separately and combine results
    # using the centreline mask.

    az_grids = dict.fromkeys(az_samples)

    for key in az_grids:
        az = az_samples[key]
        # Slice off unused elements.
        n = az['count']
        az['x'] = az['x'][:n]
        az['y'] = az['y'][:n]
        az['a'] = az['a'][:n]
        grid = interpolate_sample_points((az['x'], az['y']), az['a'], mg, 'linear')
        grid = numpy.nan_to_num(grid)
        az_grids[key] = grid
        if debug:
            logger.debug('AZIMUTH: %s %s %s', key, az_grids[key].shape, az_grids[key])
        del grid

    az = numpy.where(cmask, az_grids['R'], az_grids['L'])
    del az_grids

    logger.debug('L+R AZIMUTH: %s %s %s', az.shape, az, cmask)

    # Back fill any pixels that have been stepped over.

    i_miss, j_miss = numpy.where(az == 0)

    for n, i0 in enumerate(i_miss):
        j0 = j_miss[n]
        az[i0, j0] = az[i0, j0-1]

    logger.debug('Zero values after back fill: %s', numpy.where(az == 0))

    # Interpolate L/R satellite view angle separately, and combine results
    # using the centreline mask.

    vu_grids = dict.fromkeys(vu_samples)

    for key in vu_grids:
        vu = vu_samples[key]
        # Slice off unused elements.
        n = vu['count']
        vu['x'] = vu['x'][:n]
        vu['y'] = vu['y'][:n]
        vu['v'] = vu['v'][:n]
        grid = interpolate_sample_points((vu['x'], vu['y']), vu['v'], mg, 'linear')
        grid = numpy.nan_to_num(grid)
        vu_grids[key] = grid
        if debug:
            logger.debug('VIEW ANGLE: %s %s %s', key, vu_grids[key].shape, vu_grids[key])
        del grid

    vu = numpy.where(cmask, vu_grids['R'], vu_grids['L'])
    numpy.sqrt(vu*vu, vu)
    del vu_grids

    logger.debug('L+R VIEW ANGLE: %s %s %s', vu.shape, vu, cmask)

    # Time

    #t = bispl_sample_points(_xs, _ys, _ts, xis, yis)
    #t.astype(numpy.float32).tofile('TIME_BS.bin')
    #t = interpolate_sample_points(xy, _t, mg, 'linear')
    #t.astype(numpy.float32).tofile('TIME_GL.bin')
    #t = interpolate_sample_points(xy, _t, mg, 'cubic')
    #t.astype(numpy.float32).tofile('TIME_GC.bin')

    t = interpolate_sample_points(xy, _t, mg, 'linear')
    
    # Restore sharp discontinuity by removing 24h offset for AM times
    t[t > 24.0] -= 24.0 

    # View angle (unsigned)

    #v = bispl_sample_points(_xs, _ys, _vs, xis, yis)
    #v.astype(numpy.float32).tofile('VIEW_BS.bin')
    #v = interpolate_sample_points(xy, _v, mg, 'linear')
    #v.astype(numpy.float32).tofile('VIEW_GL.bin')
    #v = interpolate_sample_points(xy, _v, mg, 'cubic')
    #v.astype(numpy.float32).tofile('VIEW_GC.bin')

    #v = interpolate_sample_points(xy, _v, mg, 'linear')
    #numpy.sqrt(v*v, v)
    v = vu

    # Satellite zenith
    # OPTIMISATION: satellite zenith is not required for NBAR
    #z = numpy.arcsin(numpy.sin(v) * HR_PLUS_ONE)

    # Satellite azimuth: version 1
    # Centre line discontinuity causes interpolation headaches.

    #az_glin = interpolate_sample_points(xy, _a, mg, 'linear')
    #az_gcub = interpolate_sample_points(xy, _a, mg, 'cubic')
    #az_bspl = bispl_sample_points(_x, _y, _a, xis, yis, kx=3, ky=3, s=1.0)
    #
    # TEST
    #azs = bispl_sample_points(_xs, _ys, _as, xis, yis, kx=3, ky=3, s=1.0)
    #azs.astype(numpy.float32).tofile('AZ_BISPL.bin')

    # Satellite azimuth: version 2
    # Interpolate the intermediate grids
    #     dlat   = pixel to satellite nadir latitude difference
    #     gamma  = satellite nadir to pixel great circle arc
    # and use them to calculate the azimuth. Take care to maintain
    # the sign convention as used in the sampling function.

    ##dlat = interpolate_sample_points(xy, _dlat, mg, 'cubic')
    #dlat = bispl_sample_points(_xs, _ys, _dlats, xis, yis)
    ##gamma = interpolate_sample_points(xy, _gamma, mg, 'cubic')
    #gamma = bispl_sample_points(_xs, _ys, _gammas, xis, yis)
    #beta = numpy.sign(dlat) * numpy.arccos(numpy.tan(dlat) / numpy.tan(gamma))
    #az = math.pi - beta

    # Satellite azimuth: version 3
    # Equation 2.7 of Green, "Spherical Astronomy".
    # Unstable: Produces cos(A) values >> 1 close to the centre line as
    # gamma -> 0.

    #slat = interpolate_sample_points(xy, _slat, mg, 'linear')
    ##slat = bispl_sample_points(_xs, _ys, _slats, xis, yis)
    ## geocentric latitude
    #phigc = numpy.arctan((1.0 - earth.ECC2) * numpy.tan(lat))
    #az3 = ( (numpy.sin(slat) - numpy.cos(gamma) * numpy.sin(phigc)) /
    #        (numpy.sin(gamma) * numpy.cos(phigc)) )
    #az3u = numpy.arccos(az3)
    #az3s = numpy.where(numpy.sign(az3) < 0, az3u, (math.pi*2 - az3u))

    if debug:
        logger.debug('DEBUG *** eval_sat_grids()')
        logger.debug('TIME = %s', t)
        logger.debug('VIEW = %s', v)
#            logger.debug('ZENITH = %s', z)
        logger.debug('AZIMUTH (hbeta, cmask) = %s', az)

        #print 'DLAT'
        #print dlat
        #print 'GAMMA'
        #print gamma
        #print 'BETA'
        #print beta
        #print 'AZIMUTH (via beta)'
        #print az
        #print 'AZIMUTH (direct, griddata/linear)'
        #print az_glin
        #print 'AZIMUTH (direct, griddata/cubic)'
        #print az_gcub
        #print 'AZIMUTH (direct, bispl)'
        #print az_bspl
        #print 'AZIMUTH (direct, bispl, k=1)'
        #print azs

    #return (t, v, z, az)
    return (t, v, az)


