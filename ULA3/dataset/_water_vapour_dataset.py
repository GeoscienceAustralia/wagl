#!/usr/bin/env python
"""
WaterVapourDataset class for ancillary data
Subclassed from Dataset

Author: Alex Ip (alex.ip@ga.gov.au)
Includes code from Roger Edberg's original anci_geotiff_loader.py
"""
import os, logging
from datetime import datetime, timedelta
from ULA3.dataset import Dataset

logger = logging.getLogger('root.' + __name__)

class WaterVapourDataset(Dataset):
    """Subclass of Dataset to handle Water vapour ancillary data

    Additional functionality to deal with band structure in water vapour
    GeoTIFF files.
    """
    def __init__(self, wv_data_root, ancillary_datetime):
        """Initialise an AncillaryDataset instance.

        Arguments:
            wv_data_root: root directory for water vapour ancillary data
            ancillary_datetime: Python datetime object
        """
        logger.debug('WaterVapourDataset(%s, %s)', wv_data_root, ancillary_datetime)
        self.ancillary_datetime = ancillary_datetime
        pathname = os.path.join(wv_data_root,
            "pr_wtr.eatm.%s.tif" % ancillary_datetime.strftime('%Y'))

        logger.debug('Water vapour ancillary data file = %s', pathname)

        assert os.path.exists(pathname), 'Water vapour ancillary data file ' + pathname + ' does not exist'
        super(WaterVapourDataset, self).__init__(pathname)

    def get_band_index(self, ancillary_datetime=None):
        """Calculate the data band number based on the datetime string.

        Code and data layout are not documented in original extractVapour.c
        program, so we clone the code here.

        N.B: Water vapour datasets each cover a single year over the whole of Australia
        With day-of-year and hour in separate bands.

        Arguments:
            ancillary_datetime: datetime specifier. Year is already determined

        Returns:
            Band index for the specified ancillary_datetime.
        """
        logger.debug('get_band_index(%s) called', ancillary_datetime)
        ancillary_datetime = ancillary_datetime or self.ancillary_datetime

        assert ancillary_datetime.timetuple().tm_year == self.ancillary_datetime.timetuple().tm_year, 'Year mismatch'

        doy = ancillary_datetime.timetuple().tm_yday
        hour = ancillary_datetime.timetuple().tm_hour
        band_number = (int(doy) - 1) * 4 + int((hour + 3) / 6)

        # Hack to fix boundary condition: 1 Jan, 0-3 hours
        if band_number == 0 and doy == 1:
            band_number = 1

        # Note: RasterCount == tiff->samples_per_pixel

        logger.debug('band_number = %s', band_number)

	####### hack to enable NBAR NRT
	if band_number > (self.RasterCount +1):
	    rasterdoy = (((self.RasterCount)-(int((hour + 3) / 6)))/4)+1
    	    if (doy-rasterdoy) < 7:
                band_number = (int(rasterdoy) - 1) * 4 + int((hour + 3) / 6)
	    #band_number = (self.RasterCount)
	####### end of hack

        assert band_number in range(1, self.RasterCount + 1), 'Invalid Water Vapour band number: %s' % band_number
        return band_number

    def get_data_value(self, longlat, ancillary_datetime = None, nBand = None):
        """Extract water vapour value for a given latitude, longitude, and band.
        Overrides Dataset.get_data_value

        Arguments:
            longlat: tuple containing (longitude, latitude) as decimal degrees
            lat: latitude (decimal degrees)
            datetime: datetime specifier
            nBand: Optional band number

        Returns:
            Data value at the specified (lon, lat) coordinates.
        """
        logger.debug('get_data_value(%s, %s, %s) called', longlat, ancillary_datetime, nBand)

        ancillary_datetime = ancillary_datetime or self.ancillary_datetime
        nBand = nBand or self.get_band_index(ancillary_datetime)

        wv_value_raw = super(WaterVapourDataset, self).get_data_value(longlat, nBand)
        logger.debug('Raw (scaled) water vapour value for %s = %s', longlat, wv_value_raw)

        return wv_value_raw * 0.1 # Convert from kg/m^2 to g/cm^2

def main():
    """Main program for command-line utility.

    """
    import argparse

    argp = argparse.ArgumentParser('Extract ancillary data value')
    argp.add_argument('--dataroot', action='store', default=None, help='Water vapour data root')
    argp.add_argument('--lat', action='store', default=None, help='DD.XXXX')
    argp.add_argument('--lon', action='store', default=None, help='DD.XXXX')
    argp.add_argument('--datetime', action='store', default=None,
                                    help='DATETIME string = "DOY:HH.XXXX"')
    argp.add_argument('--debug', action='store_true', default=False)

    args = argp.parse_args()

    if args.debug:
        print '\n\n*** DEBUG ***'
        print 'ARGS', args

    dataroot = args.dataroot or '/short/v10/eoancillary/water_vapour'

    if args.datetime:
        # Construct a Python datetime (this year) from "DOY:HH.XXXX"
        doy, hour = args.datetime.split(':', 1)
        input_datetime = datetime(2010, 1, 1) + timedelta(int(doy) - 1 + float(hour) / 24.0)
    else:
        input_datetime = datetime(2010, 1, 28, 2, 13, 18, 492769)


    ancillaryDataset = WaterVapourDataset(wv_data_root=dataroot, ancillary_datetime=input_datetime)

    if args.lon:
        lon = float(args.lon)
    else:
        lon = 114.1056

    if args.lat:
        lat = float(args.lat)
    else:
        lat = -23.1045

    data_value = ancillaryDataset.get_data_value((lon, lat))

    print
    print 'Water vapour:', data_value


if __name__ == '__main__':
    main()
