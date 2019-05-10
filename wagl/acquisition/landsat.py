"""
Defines the acquisition classes for the landsat satellite program for wagl
"""
import math

from .base import Acquisition


class LandsatAcquisition(Acquisition):

    """A Landsat acquisition."""

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        self.min_radiance = 0
        self.max_radiance = 1
        self.min_quantize = 0
        self.max_quantize = 1
        self.__data = {}  # Imagery data cache

        super(LandsatAcquisition, self).__init__(pathname, uri,
                                                 acquisition_datetime,
                                                 band_name=band_name,
                                                 band_id=band_id,
                                                 metadata=metadata)

        self._gain = ((self.max_radiance - self.min_radiance) /
                      (self.max_quantize - self.min_quantize))
        self._bias = self.max_radiance - (self.gain * self.max_quantize)

    @property
    def gain(self):
        """
        A multiplier used for scaling the data.
        """
        return self._gain

    @property
    def bias(self):
        """
        An additive used for scaling the data.
        """
        return self._bias

    def data(self, out=None, window=None, masked=False):
        """
        Retrieves data from source imagery or internal cache if tar
        """
        file_suffix = self.uri.split('!')[0].split('.', 1)
        
        # Check if source imagery directly accessible
        if len(file_suffix) == 1 or 'tar' not in file_suffix[1]:
            return super().data(out, window, masked)

        # Check if source imagery is cached
        if self.__data.get((masked, )) is None:
            self.__data[(masked, )] = super().data(masked=masked)

        # Retrieve data from cache
        if window:
           out = self.__data[(masked, )][
              window[0][0]:window[0][1], window[1][0]:window[1][1]].copy()
        else:
           out = self.__data[(masked, )].copy()

        return out

    def radiance_data(self, window=None, out_no_data=-999, esun=None):
        """
        Return the data as radiance in watts/(m^2*micrometre).
        """
        data = self.data(window=window)

        # check for no data
        no_data = self.no_data if self.no_data is not None else 0
        nulls = data == no_data

        # gain & offset; y = mx + b
        radiance = self.gain * data + self.bias

        # set the out_no_data value inplace of the input no data value
        radiance[nulls] = out_no_data

        return radiance

    def close(self):
        """ Clears acquisition level cache """
        self.__data = {}
        super().close()


class Landsat5Acquisition(LandsatAcquisition):

    """ Landsat 5 acquisition. """

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Landsat5Acquisition, self).__init__(pathname, uri,
                                                  acquisition_datetime,
                                                  band_name=band_name,
                                                  band_id=band_id,
                                                  metadata=metadata)

        self.platform_id = 'LANDSAT_5'
        self.sensor_id = 'TM'
        self.tle_format = 'l5_%4d%s_norad.txt'
        self.tag = 'LS5'
        self.altitude = 705000.0
        self.inclination = 1.7139133254584316445390643346558
        self.omega = 0.001059
        self.radius = 7285600.0
        self.semi_major_axis = 7083160.0
        self.maximum_view_angle = 9.0

        self._norad_id = 14780
        self._classification_type = 'U'
        self._international_designator = '84021A'


class Landsat7Acquisition(LandsatAcquisition):

    """ Landsat 7 acquisition. """

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Landsat7Acquisition, self).__init__(pathname, uri,
                                                  acquisition_datetime,
                                                  band_name=band_name,
                                                  band_id=band_id,
                                                  metadata=metadata)

        self.platform_id = 'LANDSAT_7'
        self.sensor_id = 'ETM+'
        self.tle_format = 'L7%4d%sASNNOR.S00'
        self.tag = 'LS7'
        self.altitude = 705000.0
        self.inclination = 1.7139133254584316445390643346558
        self.omega = 0.001059
        self.radius = 7285600.0
        self.semi_major_axis = 7083160.0
        self.maximum_view_angle = 9.0

        self._norad_id = 25682
        self._classification_type = 'U'
        self._international_designator = '99020A'


class Landsat8Acquisition(LandsatAcquisition):

    """ Landsat 8 acquisition. """

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Landsat8Acquisition, self).__init__(pathname, uri,
                                                  acquisition_datetime,
                                                  band_name=band_name,
                                                  band_id=band_id,
                                                  metadata=metadata)

        self.platform_id = 'LANDSAT_8'
        self.sensor_id = 'OLI'
        self.tle_format = 'L8%4d%sASNNOR.S00'
        self.tag = 'LS8'
        self.altitude = 705000.0
        self.inclination = 1.7139133254584316445390643346558
        self.omega = 0.001059
        self.radius = 7285600.0
        self.semi_major_axis = 7083160.0
        self.maximum_view_angle = 9.0

        self._norad_id = 39084
        self._classification_type = 'U'
        self._international_designator = '13008A'

    def radiance_data(self, window=None, out_no_data=-999, esun=None):
        """
        This method overwrites the parent's method 'radiance_data' for Landsat8
        acquistions to compute radiance by first calculating top of the
        atmosphere reflectance and then converting to radiance using 'esun" value
        if it is not None. If esun is None, parents (LandsatAcquisition)'s
        radiance_data computation is used.

        Return the data as radiance in watts/(m^2*micrometre).
        """

        if not esun:
            return super(Landsat8Acquisition, self).radiance_data(window, out_no_data, esun)

        data = self.data(window=window)

        # check for no data
        no_data = self.no_data if self.no_data is not None else 0
        nulls = data == no_data

        print(data)

        gain = self.reflectance_mult
        bias = self.reflectance_add

        print(gain, bias, esun)
        toa_reflectance = gain * data + bias

        radiance = toa_reflectance * esun / math.pi

        # set the out_no_data value inplace of the input no data value
        radiance[nulls] = out_no_data

        return radiance


ACQUISITION_TYPE = {
    'Landsat5_TM': Landsat5Acquisition,
    'Landsat7_ETM+': Landsat7Acquisition,
    'LANDSAT_5_TM': Landsat5Acquisition,
    'LANDSAT_7_ETM+': Landsat7Acquisition,
    'LANDSAT_8_OLI': Landsat8Acquisition,
    'LANDSAT_8_OLI_TIRS': Landsat8Acquisition
}
