"""
Defines the acquisition classes for the sentinel satellite program for wagl
"""

from xml.etree import ElementTree
import zipfile
import os

from dateutil import parser
import numexpr
import pyproj

import numpy
import pandas
from scipy import interpolate

from .base import Acquisition

def find_all_in(path, s):
    """
    Search through `path` and its children for all occurances of
    files with `s` in their name. Returns the (possibly empty) list
    of file paths
    """
    result = []
    for root, _, files in os.walk(path):
        for f in files:
            if s in f:
                result.append(os.path.join(root, f))
    return result


def s2_index_to_band_id(band_index):
    """s2_index_toBand_id returns the band_id from the band index 

    :param band_index: band index (0-12 inclusive) referencing sensors in ESA's metadata
    :return: band_id for the band_index
    """

    return {
        0: '1', 1: '2', 2: '3', 3: '4', 4: '5', 5: '6', 6: '7',
        7: '8', 8: '8A', 9: '9', 10: '10', 11: '11', 12: '12'
    }[int(band_index)]


class Sentinel2Acquisition(Acquisition):

    """ Base class for a Sentinel-2 acquisition. """

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):

        super(Sentinel2Acquisition, self).__init__(pathname, uri,
                                                   acquisition_datetime,
                                                   band_name=band_name,
                                                   band_id=band_id,
                                                   metadata=metadata)

    def _get_gps_xml(self):
        """Returns in memory XML tree for gps coordinates"""
        # open the zip archive and get the xml root
        archive = zipfile.ZipFile(self.pathname)
        xml_files = [s for s in archive.namelist() if
                     ("DATASTRIP" in s) & (".xml" in s)]

        # there could be several matches; loop till we find one with GPS data
        for xml_file in xml_files:
            xml_root = ElementTree.XML(archive.read(xml_file))

            gps_list = xml_root.findall('./*/Ephemeris/GPS_Points_List')
            if gps_list:
                break

        try:
            gps = gps_list[0]
        except IndexError:
            # we can't process without the gps file, as we aren't
            # collecting the TLE data for S2A.
            msg = "No GPS data found."
            raise Exception(msg)

        return gps

    def read_gps_file(self):
        """
        Returns the recorded gps data as a `pandas.DataFrame`.
        """
        # coordinate transform
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

        # output container
        data = {'longitude': [],
                'latitude': [],
                'altitude': [],
                'timestamp': []}

        gps = self._get_gps_xml()
        # there are a few columns of data that could be of use
        # but for now, just get the location and timestamp from the
        # gps points list
        for point in gps.iter('GPS_Point'):
            x, y, z = [float(i) / 1000 for i in
                       point.findtext('POSITION_VALUES').split()]
            gps_time = parser.parse(point.findtext('GPS_TIME'))

            # coordinate transformation
            lon, lat, alt = pyproj.transform(ecef, lla, x, y, z, radians=False)

            data['longitude'].append(lon)
            data['latitude'].append(lat)
            data['altitude'].append(alt)
            data['timestamp'].append(gps_time)

        return pandas.DataFrame(data)

    def _get_solar_zenith_xml(self):
        """Returns an in memory XML tree for the granule to retrieve solar zenith"""
        archive = zipfile.ZipFile(self.pathname)
        xml_root = ElementTree.XML(archive.read(self.granule_xml))
        return xml_root

    def _retrieve_solar_zenith(self):
        """
        We can't use our own zenith angle to invert TOAr back to
        radiance, even though our angle array is more accurate.
        This is because the correct radiance measurement won't be
        guaranteed if a different value is used in the inversion.

        Code adapted from https://github.com/umwilm/SEN2COR.
        """
        def rbspline(y_coords, x_coords, zdata):
            """
            A wrapper providing the call signiture for RectBivariateSpline.
            Code adapted from https://github.com/umwilm/SEN2COR.
            """
            y = numpy.arange(zdata.shape[0], dtype=numpy.float32)
            x = numpy.arange(zdata.shape[1], dtype=numpy.float32)
            func = interpolate.RectBivariateSpline(x, y, zdata)

            return func(y_coords, x_coords)

        xml_root = self._get_solar_zenith_xml()

        # read the low res solar zenith data
        search_term = './*/Tile_Angles/Sun_Angles_Grid/Zenith/Values_List'
        values = xml_root.findall(search_term)[0]

        dims = (len(values), len(values[0].text.split()))
        data = numpy.zeros(dims, dtype='float32')
        for i, val in enumerate(values.iter('VALUES')):
            data[i] = val.text.split()

        # correct solar_zenith dimensions
        if self.lines < self.samples:
            last_row = int(data[0].size * float(self.lines) /
                           float(self.samples) + 0.5)
            solar_zenith = data[0:last_row, :]
        elif self.samples < self.lines:
            last_col = int(data[1].size * float(self.samples) /
                           float(self.lines) + 0.5)
            solar_zenith = data[:, 0:last_col]
        else:
            solar_zenith = data

        numpy.absolute(solar_zenith, out=solar_zenith)
        numpy.clip(solar_zenith, 0, 70.0, out=solar_zenith)

        # interpolate across the full dimensions of the acquisition
        dims = solar_zenith.shape
        y = numpy.arange(self.lines) / (self.lines - 1) * dims[0]
        x = numpy.arange(self.samples) / (self.samples - 1) * dims[1]

        solar_zenith = numpy.float32(rbspline(y, x, solar_zenith))
        self._solar_zenith = numpy.radians(solar_zenith, out=solar_zenith)

    def radiance_data(self, window=None, out_no_data=-999):
        """
        Return the data as radiance in watts/(m^2*micrometre).

        Sentinel-2a's package is a little convoluted with the various
        different scale factors, and the code for radiance inversion
        doesn't follow the general standard.
        Hence the following code may look a little strange.

        Code adapted from https://github.com/umwilm/SEN2COR.
        """
        # retrieve the solar zenith if we haven't already done so
        if self._solar_zenith is None:
            self._retrieve_solar_zenith()

        # Python style index
        if window is None:
            idx = (slice(None, None), slice(None, None))
        else:
            idx = (slice(window[0][0], window[0][1]),
                   slice(window[1][0], window[1][1]))

        # coefficients
        # pylint: disable=unused-argument,unused-variable
        sf = numpy.float32(1 / (self.c1 * self.qv))
        pi_d2 = numpy.float32(numpy.pi * self.d2)
        esun = numpy.float32(self.solar_irradiance / 10)
        solar_zenith = self._solar_zenith[idx]
        rsf = numpy.float32(self.radiance_scale_factor)

        # toa reflectance
        data = self.data(window=window)

        # check for no data
        no_data = self.no_data if self.no_data is not None else 0
        nulls = data == no_data

        # inversion
        expr = "((data * esun * cos(solar_zenith) * sf) / pi_d2) * rsf"
        radiance = numexpr.evaluate(expr)
        radiance[nulls] = out_no_data

        return radiance

    def close(self):
        """
        Set self._solar_zenith back to None to reclaim memory.
        """
        self._solar_zenith = None


class Sentinel2aAcquisition(Sentinel2Acquisition):

    """ Sentinel-2a acquisition. """

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Sentinel2aAcquisition, self).__init__(pathname, uri,
                                                    acquisition_datetime,
                                                    band_name=band_name,
                                                    band_id=band_id,
                                                    metadata=metadata)

        self.platform_id = 'SENTINEL_2A'
        self.sensor_id = 'MSI'
        self.tle_format = 'S2A%4d%sASNNOR.S00'
        self.tag = 'S2A'
        self.altitude = 786000.0
        self.inclination = 1.721243708316808
        self.omega = 0.001039918
        self.semi_major_axis = 7167000.0
        self.maximum_view_angle = 20.0

        self._norad_id = 40697
        self._classification_type = 'U'
        self._international_designator = '15028A'

        self._gps_file = True
        self._solar_zenith = None


class Sentinel2bAcquisition(Sentinel2Acquisition):

    """ Sentinel-2b acquisition. """

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Sentinel2bAcquisition, self).__init__(pathname, uri,
                                                    acquisition_datetime,
                                                    band_name=band_name,
                                                    band_id=band_id,
                                                    metadata=metadata)

        self.platform_id = 'SENTINEL_2B'
        self.sensor_id = 'MSI'
        self.tle_format = 'S2B%4d%sASNNOR.S00'
        self.tag = 'S2B'
        self.altitude = 786000.0
        self.inclination = 1.721243708316808
        self.omega = 0.001039918
        self.semi_major_axis = 7167000.0
        self.maximum_view_angle = 20.0

        self._norad_id = 42063
        self._classification_type = 'U'
        self._international_designator = '17013A'

        self._gps_file = True
        self._solar_zenith = None


class Sentinel2aSinergiseAcquisition(_Sentinel2SinergiseAcquisition, Sentinel2aAcquisition):
    pass

class Sentinel2bSinergiseAcquisition(_Sentinel2SinergiseAcquisition, Sentinel2bAcquisition):
    pass


class _Sentinel2SinergiseAcquisition(Sentinel2Acquisition):

    """
    _Sentinel2SinergiseAcquisition provides a customisation on retrieval
    for the ancillary information that is retrieved from AWS S3 Bucket
    (s3://sentinel-s2-l1c).

    It is assumed that the path a local directory containing the granule
    location is provided as the acquisition location and that the metadata
    for the datastrip is contained in a subdirectory of called 'datastrip'.
    """

    def _get_gps_xml(self):
        """Returns in memory XML tree for gps coordinates"""
        xml_files = [s for s in find_all_in(self.pathname, 'metadata.xml') if
                     'datastrip' in s]

        # there could be several matches; loop till we find one with GPS data
        for xml_file in xml_files:
            xml_root = None
            with open(xml_file, 'rb') as fd:
                xml_root = ElementTree.XML(fd.read())

            gps_list = xml_root.findall('./*/Ephemeris/GPS_Points_List')
            if gps_list:
                break

        try:
            gps = gps_list[0]
        except IndexError:
            # we can't process without the gps file, as we aren't
            # collecting the TLE data for S2A.
            msg = "No GPS data found."
            raise Exception(msg)

        return gps

    def _get_solar_zenith_xml(self):
        """Returns an in memory XML tree for the granule to retrieve solar zenith"""
        xml_root = None
        with open(self.granule_xml, 'rb') as fd:
            xml_root = ElementTree.XML(fd.read())
        return xml_root
