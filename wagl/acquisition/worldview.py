"""
Acquisitions for WorldView-2 satellite.
"""

import rasterio

from .base import Acquisition


class WorldView2Acquisition(Acquisition):
    """ A WorldView-2 acquisition. """

    platform_id = "WORLDVIEW_2"

    def __init__(
        self,
        pathname,
        uri,
        acquisition_datetime,
        band_name="BAND-C",
        band_id="1",
        metadata=None,
    ):
        super().__init__(
            pathname,
            uri,
            acquisition_datetime,
            band_name=band_name,
            band_id=band_id,
            metadata=metadata,
        )

        self.tag = "WV2"

        self._norad_id = 35946
        self.altitude = 7000000.0
        self.semi_major_axis = 7144000.0
        self._international_designator = "09055A"
        self.inclination = 1.7174
        self.omega = 0.0010451
        self._classification_type = "U"

        self.maximum_view_angle = 20.0

    def close(self):
        super().close()


class WorldView2MultiAcquisition(WorldView2Acquisition):
    """ A multi-band WorldView-2 acquisition. """

    band_names = [
        "BAND-{}".format(name) for name in ["C", "B", "G", "Y", "R", "RE", "N", "N2"]
    ]
    sensor_id = "MUL"

    def __init__(
        self,
        pathname,
        uri,
        acquisition_datetime,
        band_name="BAND-C",
        band_id="1",
        metadata=None,
    ):
        super().__init__(
            pathname,
            uri,
            acquisition_datetime,
            band_name=band_name,
            band_id=band_id,
            metadata=metadata,
        )

    def close(self):
        super().close()

    def data(self, out=None, window=None, masked=False):
        """
        Return `numpy.array` of the data for this acquisition.
        If `out` is supplied, it must be a numpy.array into which
        the Acquisition's data will be read.
        """
        with rasterio.open(self.uri) as ds:
            data = ds.read(int(self.band_id), out=out, window=window, masked=masked)

        return data

    def radiance_data(self, window=None, out_no_data=-999, esun=None):
        """
        Return the data as radiance in watts/(m^2*micrometre).
        """
        data = self.data(window=window)

        # check for no data
        no_data = self.no_data if self.no_data is not None else 0
        nulls = data == no_data

        radiance = (
            self.gain * data * (self.abs_cal_factor / self.effective_bandwidth)
            + self.offset
        )

        # set the out_no_data value inplace of the input no data value
        radiance[nulls] = out_no_data

        return radiance
