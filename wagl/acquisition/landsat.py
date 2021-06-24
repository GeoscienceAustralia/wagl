"""
Defines the acquisition classes for the landsat satellite program for wagl
"""
import gzip
import math
from pathlib import Path
import tarfile
import tempfile
import numpy
import rasterio

from .base import Acquisition


class LandsatAcquisition(Acquisition):

    """A Landsat acquisition."""

    def __init__(
        self,
        pathname,
        uri,
        acquisition_datetime,
        band_name="BAND 1",
        band_id="1",
        metadata=None,
    ):
        self.min_radiance = 0
        self.max_radiance = 1
        self.min_quantize = 0
        self.max_quantize = 1
        self.__data = {}  # Imagery data cache

        super(LandsatAcquisition, self).__init__(
            pathname,
            uri,
            acquisition_datetime,
            band_name=band_name,
            band_id=band_id,
            metadata=metadata,
        )

        self._gain = (self.max_radiance - self.min_radiance) / (
            self.max_quantize - self.min_quantize
        )
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
        file_suffix = self.uri.split("!")[0].split(".", 1)

        # Check if source imagery directly accessible
        if len(file_suffix) == 1 or "tar" not in file_suffix[1]:
            return super().data(out, window, masked)

        # Check if source imagery is cached
        if self.__data.get((masked,)) is None:
            self.__data[(masked,)] = super().data(masked=masked)

        # Retrieve data from cache
        if window:
            out = self.__data[(masked,)][
                window[0][0] : window[0][1], window[1][0] : window[1][1]
            ].copy()
        else:
            out = self.__data[(masked,)].copy()

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
        """Clears acquisition level cache"""
        self.__data = {}
        super().close()


class Landsat5Acquisition(LandsatAcquisition):

    """Landsat 5 acquisition."""

    def __init__(
        self,
        pathname,
        uri,
        acquisition_datetime,
        band_name="BAND 1",
        band_id="1",
        metadata=None,
    ):
        super(Landsat5Acquisition, self).__init__(
            pathname,
            uri,
            acquisition_datetime,
            band_name=band_name,
            band_id=band_id,
            metadata=metadata,
        )

        self.platform_id = "LANDSAT_5"
        self.sensor_id = "TM"
        self.tle_format = "l5_%4d%s_norad.txt"
        self.tag = "LS5"
        self.altitude = 705000.0
        self.inclination = 1.7139133254584316445390643346558
        self.omega = 0.001059
        self.radius = 7285600.0
        self.semi_major_axis = 7083160.0
        self.maximum_view_angle = 9.0

        self._norad_id = 14780
        self._classification_type = "U"
        self._international_designator = "84021A"


class Landsat7Acquisition(LandsatAcquisition):

    """Landsat 7 acquisition."""

    def __init__(
        self,
        pathname,
        uri,
        acquisition_datetime,
        band_name="BAND 1",
        band_id="1",
        metadata=None,
    ):
        super(Landsat7Acquisition, self).__init__(
            pathname,
            uri,
            acquisition_datetime,
            band_name=band_name,
            band_id=band_id,
            metadata=metadata,
        )

        self.platform_id = "LANDSAT_7"
        self.sensor_id = "ETM+"
        self.tle_format = "L7%4d%sASNNOR.S00"
        self.tag = "LS7"
        self.altitude = 705000.0
        self.inclination = 1.7139133254584316445390643346558
        self.omega = 0.001059
        self.radius = 7285600.0
        self.semi_major_axis = 7083160.0
        self.maximum_view_angle = 9.0

        self._norad_id = 25682
        self._classification_type = "U"
        self._international_designator = "99020A"

        self._gap_mask = None

    def _extract_gap_mask(self):
        """
        Extract the gap mask from within the tar.
        The gap mask is a TIF that is then pushed through gzip,
        horrible, and in order to read it, we need to unpack from
        the tar, decompress using gzip, then write to disk so GDAL
        can read it.
        """
        # mask files are contained in a sub-directory named 'gap-mask'
        path = Path(self.uri.replace("!", "")[6:])
        parts = path.name.split("_")

        # non-thermal bands are named
        # LE07_L1TP_092084_20110809_20161206_01_T1_B4.TIF
        # thermal bands are named
        # LE07_L1TP_092084_20110809_20161206_01_T1_B6_VCID_1.TIF
        # so insert before index 7
        # (not ideal, but USGS could change the whole convention anyway
        parts.insert(7, "GM")
        mask_name = Path("gap_mask", "{}.gz".format("_".join(parts)))

        # open tarfile
        with tarfile.open(str(path.parent)) as tf:
            mem = tf.getmember(str(mask_name))
            fobj = tf.extractfile(mem)
            with gzip.open(fobj) as gz:
                with tempfile.TemporaryDirectory(suffix=".gap-mask") as tmpd:
                    out_fname = Path(tmpd, "gap-mask.tif")
                    with open(out_fname, "wb") as src:
                        src.write(gz.read())

                    # read gap mask into memory
                    with rasterio.open(out_fname) as src:
                        data = src.read(1)

                    self._gap_mask = data == 0

    def radiance_data(self, window=None, out_no_data=-999, esun=None):
        """
        This method overwrites the parent's method to handle a special
        case for Landsat-7. The data supplied by USGS contains ~2pixel
        interpolation for the SLC-OFF data. We need to read the
        gap_mask file in order to mask out the interpolation fill.
        Return the data as radiance in watts/(m^2*micrometre) with
        interpolated pixels masked out.
        """
        # retrieve the gap mask if we haven't already done so
        if self._gap_mask is None:
            try:
                self._extract_gap_mask()
            except (FileNotFoundError, KeyError):
                # we might be dealing with an acquisition that is not
                # a tarfile (legacy), or an acquisition with no gap masks
                self._gap_mask = numpy.zeros((self.lines, self.samples), dtype="bool")

        # Python style index
        if window is None:
            idx = (slice(None, None), slice(None, None))
        else:
            idx = (slice(window[0][0], window[0][1]), slice(window[1][0], window[1][1]))

        data = self.data(window=window)

        # check for no data
        no_data = self.no_data if self.no_data is not None else 0
        nulls = data == no_data

        # read same block for the mask
        mask = self._gap_mask[idx]

        # gain & offset; y = mx + b
        radiance = self.gain * data + self.bias

        # set the out_no_data value inplace of the input no data value
        radiance[nulls] = out_no_data
        radiance[mask] = out_no_data

        return radiance

    def close(self):
        """Clears acquisition level cache"""
        self._gap_mask = None
        self.__data = {}
        super().close()


class Landsat8Acquisition(LandsatAcquisition):

    """Landsat 8 acquisition."""

    def __init__(
        self,
        pathname,
        uri,
        acquisition_datetime,
        band_name="BAND 1",
        band_id="1",
        metadata=None,
    ):
        super(Landsat8Acquisition, self).__init__(
            pathname,
            uri,
            acquisition_datetime,
            band_name=band_name,
            band_id=band_id,
            metadata=metadata,
        )

        self.platform_id = "LANDSAT_8"
        self.sensor_id = "OLI"
        self.tle_format = "L8%4d%sASNNOR.S00"
        self.tag = "LS8"
        self.altitude = 705000.0
        self.inclination = 1.7139133254584316445390643346558
        self.omega = 0.001059
        self.radius = 7285600.0
        self.semi_major_axis = 7083160.0
        self.maximum_view_angle = 9.0

        self._norad_id = 39084
        self._classification_type = "U"
        self._international_designator = "13008A"

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
            return super(Landsat8Acquisition, self).radiance_data(
                window, out_no_data, esun
            )

        data = self.data(window=window)

        # check for no data
        no_data = self.no_data if self.no_data is not None else 0
        nulls = data == no_data

        gain = self.reflectance_mult
        bias = self.reflectance_add

        toa_reflectance = gain * data + bias

        radiance = toa_reflectance * esun / math.pi

        # set the out_no_data value inplace of the input no data value
        radiance[nulls] = out_no_data

        return radiance


ACQUISITION_TYPE = {
    "Landsat5_TM": Landsat5Acquisition,
    "Landsat7_ETM+": Landsat7Acquisition,
    "LANDSAT_5_TM": Landsat5Acquisition,
    "LANDSAT_7_ETM+": Landsat7Acquisition,
    "LANDSAT_8_OLI": Landsat8Acquisition,
    "LANDSAT_8_OLI_TIRS": Landsat8Acquisition,
}
