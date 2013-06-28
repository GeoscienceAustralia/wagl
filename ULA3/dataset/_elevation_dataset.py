#!/usr/bin/env python
"""
ElevationDataset class for ancillary data
Subclassed from Dataset

Author: Alex Ip (alex.ip@ga.gov.au)
"""
import os, logging
from datetime import datetime, timedelta

from . import Dataset

logger = logging.getLogger('root.' + __name__)

class ElevationDataset(Dataset):
    """
    Convenience class for handling elevation ancillary data

    """
    def __init__(self, dem_data_root):
        """
        Initialise an AncillaryDataset instance.

        :param dem_data_root:
            root directory for Elevation ancillary data

        """
        pathname = os.path.join(dem_data_root, "DEM_one_deg.tif")

        assert os.path.exists(pathname), 'Elevation ancillary data file ' + pathname + ' does not exist'
        super(ElevationDataset, self).__init__(pathname)

    def get_data_value(self, longlat):
        """Extract Elevation value for a given latitude & longitude.
        Overrides Dataset.get_data_value

        Arguments:
            longlat: tuple containing (longitude, latitude) as decimal degrees

        Returns:
            Data value at the specified (lon, lat) coordinates.
        """

        dem_value_raw = super(ElevationDataset, self).get_data_value(longlat)

        return dem_value_raw * 0.001  # Apply scale to return correct units





def main():
    """Main program for command-line utility.

    """
    import argparse

    argp = argparse.ArgumentParser('Extract ancillary data value')
    argp.add_argument('--dataroot', action='store', default=None, help='Elevation data root')
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


    ancillaryDataset = ElevationDataset(wv_data_root=dataroot, ancillary_datetime=input_datetime)

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
    print 'Elevation:', data_value


if __name__ == '__main__':
    main()
