import unittest as ut
from ULA3.dataset import OzoneDataset
from datetime import datetime, timedelta

class OzoneDatasetTestCase(ut.TestCase):
    """
    This was just copied from the bottom of :py:mod:`ULA3.dataset._ozone_dataset`
    (where it was the body of a function ``main``, which was called when the module was run as a script).
    """

    import argparse

    argp = argparse.ArgumentParser('Extract ancillary data value')
    argp.add_argument('--dataroot', action='store', default=None, help='Ozone data root')
    argp.add_argument('--lat', action='store', default=None, help='DD.XXXX')
    argp.add_argument('--lon', action='store', default=None, help='DD.XXXX')
    argp.add_argument('--datetime', action='store', default=None,
                                    help='DATETIME string = "DOY:HH.XXXX"')
    argp.add_argument('--debug', action='store_true', default=False)

    args = argp.parse_args()

    if args.debug:
        print '\n\n*** DEBUG ***'
        print 'ARGS', args

    dataroot = args.dataroot or '/short/v10/eoancillary/ozone'

    if args.datetime:
        # Construct a Python datetime (this year) from "DOY:HH.XXXX"
        doy, hour = args.datetime.split(':', 1)
        input_datetime = datetime(2010, 1, 1) + timedelta(int(doy) - 1 + float(hour) / 24.0)
    else:
        input_datetime = datetime(2010, 1, 28, 2, 13, 18, 492769)


    ancillaryDataset = OzoneDataset(ozone_data_root=dataroot, ancillary_datetime=input_datetime)

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
    print 'Ozone:', data_value