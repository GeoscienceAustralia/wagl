"""
Core code.
"""
import os
import re
import gaip
import copy
import json
import datetime

from functools import total_ordering
from os.path import isdir, join as pjoin, dirname

REF, THM, PAN, ATM, BQA = range(5)

BAND_TYPE = {
    'Reflective': REF,
    'Thermal': THM,
    'Panchromatic': PAN,
    'Atmosphere': ATM,
    'Quality': BQA
}

with open(pjoin(dirname(__file__), 'sensors.json')) as fo:
    SENSORS = json.load(fo)


@total_ordering
class Acquisition(object):

    """Acquisition information.
    """

    def __init__(self, metadata):
        for k, v in metadata.iteritems():
            self.__dict__.update(v)

    def __eq__(self, other):
        return self.band_name == other.band_name

    def __lt__(self, other):
        return self.band_name < other.band_name

    def __repr__(self):
        return 'Acquisition(band_name=' + self.band_name + ')'

    def data(self):
        return gaip.data(self)

    def data_and_box(self):
        return gaip.data_and_box(self)

    def gridded_geo_box(self):
        return gaip.gridded_geo_box(self)





class LandsatAcquisition(Acquisition):

    def __init__(self, metadata):
        super(LandsatAcquisition, self).__init__(metadata)

    @property
    def samples(self):
        if self.band_type == REF:
            return self.product_samples_ref
        if self.band_type == THM:
            return self.product_samples_thm
        if self.band_type == PAN:
            return self.product_samples_pan

    @property
    def lines(self):
        if self.band_type == REF:
            return self.product_lines_ref
        if self.band_type == THM:
            return self.product_lines_thm
        if self.band_type == PAN:
            return self.product_lines_pan

    @property
    def grid_cell_size(self):
        if self.band_type == REF:
            return self.grid_cell_size_ref
        if self.band_type == THM:
            return self.grid_cell_size_thm
        if self.band_type == PAN:
            return self.grid_cell_size_pan

    @property
    def height(self):
        return self.samples

    @property
    def width(self):
        return self.lines

    @property
    def resolution(self):
        return self.grid_cell_size

    @property
    def scene_center_datetime(self):
        return datetime.datetime.combine(self.acquisition_date, self.scene_center_time)

    @property
    def min_radiance(self):
        return self.lmin

    @property
    def max_radiance(self):
        return self.lmax

    @property
    def min_reflectance(self):
        return getattr(self, 'qcalmin')

    @property
    def max_reflectance(self):
        return getattr(self, 'qcalmax')


class Landsat5Acquisition(LandsatAcquisition):

    def __init__(self, metadata):
        super(Landsat5Acquisition, self).__init__(metadata)

    @property
    def scene_center_time(self):
        return self.scene_center_scan_time

    @property
    def date_acquired(self):
        return self.acquisition_date


class Landsat7Acquisition(LandsatAcquisition):

    def __init__(self, metadata):
        super(Landsat7Acquisition, self).__init__(metadata)

    @property
    def scene_center_time(self):
        return self.scene_center_scan_time

    @property
    def date_acquired(self):
        return self.acquisition_date


class Landsat8Acquisition(LandsatAcquisition):

    def __init__(self, metadata):
        super(Landsat8Acquisition, self).__init__(metadata)

    @property
    def samples(self):
        if self.band_type == REF:
            return self.reflective_samples
        if self.band_type == ATM:
            return self.reflective_samples
        if self.band_type == BQA:
            return self.reflective_samples
        if self.band_type == PAN:
            return self.panchromatic_samples

    @property
    def lines(self):
        if self.band_type == REF:
            return self.reflective_lines
        if self.band_type == ATM:
            return self.reflective_lines
        if self.band_type == BQA:
            return self.reflective_lines
        if self.band_type == PAN:
            return self.panchromatic_lines

    @property
    def grid_cell_size(self):
        if self.band_type == REF:
            return self.grid_cell_size_reflective
        if self.band_type == ATM:
            return self.grid_cell_size_reflective
        if self.band_type == BQA:
            return self.grid_cell_size_reflective
        if self.band_type == PAN:
            return self.grid_cell_size_panchromatic

    @property
    def acquisition_date(self):
        return self.date_acquired

    @property
    def min_radiance(self):
        return getattr(self, 'radiance_minimum')

    @property
    def max_radiance(self):
        return getattr(self, 'radiance_maximum')

    @property
    def lmin(self):
        """The spectral radiance that is scaled to QCALMIN in
        watts/(meter squared * ster * micrometers). """
        return getattr(self, 'radiance_minimum')

    @property
    def lmax(self):
        """The spectral radiance that is scaled to QCALMAX in
        watts/(meter squared * ster * micrometers). """
        return getattr(self, 'radiance_maximum')

    @property
    def qcalmin(self):
        """The minimum quantized calibrated pixel value."""
        return getattr(self, 'quantize_cal_min')

    @property
    def qcalmax(self):
        """The maximum quantized calibrated pixel value."""
        return getattr(self, 'quantize_cal_max')

    @property
    def min_reflectance(self):
        return getattr(self, 'reflectance_minimum')

    @property
    def max_reflectance(self):
        return getattr(self, 'reflectance_maximum')

    @property
    def zone_number(self):
        return getattr(self, 'utm_zone')


ACQUISITION_TYPE = {
    'Landsat5_TM': Landsat5Acquisition,
    'Landsat7_ETM+': Landsat7Acquisition,
    'LANDSAT_8_OLI': Landsat8Acquisition
}


def find_in(path, s):
    """Search through `path` and its children for the first occurance of a
    file with `s` in its name. Returns the path of the file or `None`. """
    for root, dirs, files in os.walk(path):
        for f in files:
            if s in f:
                return os.path.join(root, f)
    return None


def acquisitions(path):
    """Obtain a list of Acquisition objects from `path`. The argument `path`
    can be a MTL file or a directory name. If `path` is a directory then the 
    MTL file will be search for in the directory and its children."""

    if os.path.isdir(path):
        filename = find_in(path, 'MTL')
    else:
        filename = path

    if filename is None:
        raise OSError("Cannot find MTL file in %s" % path)

    dirname = os.path.dirname(os.path.abspath(filename))

    data = gaip.load_mtl(filename)
    bandfiles = [k for k in data['PRODUCT_METADATA'].keys() if 'band' in k
                 and 'file_name' in k]
    bands = [b.replace('file_name', '').strip('_') for b in bandfiles]

    # We now create an acquisition object for each band and make the
    # parameters names nice.

    acquisitions = []
    for band in bands:

        # create a new copy
        new = copy.deepcopy(data)

        # remove unnecessary values
        for kv in new.values():
            for k in kv.keys():
                if band in k:
                    rm = [k.replace(band, b) for b in bands if b != band]
                    for r in rm:
                        try:
                            del kv[r]
                        except KeyError:
                            pass
                    # replace name
                    newkey = k.replace(band, '').strip('_')
                    kv[newkey] = kv[k]
                    del kv[k]

        # set path
        new['PRODUCT_METADATA']['dir_name'] = dirname

        # set band name and number
        new['PRODUCT_METADATA']['band_name'] = band
        try:
            band_num = int(band.replace('band', '').strip('_'))
            new['PRODUCT_METADATA']['band_num'] = band_num
        except ValueError:
            new['PRODUCT_METADATA']['band_num'] = None

        product = new['PRODUCT_METADATA']
        spacecraft = product['spacecraft_id']
        sensor = product['sensor_id']

        new['SPACECRAFT'] = {}
        db = SENSORS[spacecraft]
        for k, v in db.iteritems():
            if k is not 'sensors':
                new['SPACECRAFT'][str(k)] = v

        new['SENSOR_INFO'] = {}
        db = db['sensors'][sensor]
        for k, v in db.iteritems():
            if k is not 'bands':
                new['SENSOR_INFO'][str(k)] = v

        bandname = band.replace('band', '').strip('_')
        new['BAND_INFO'] = {}
        db = db['bands'][bandname]
        for k, v in db.iteritems():
            new['BAND_INFO'][str(k)] = v
        band_type = db['type_desc']
        new['BAND_INFO']['band_type'] = BAND_TYPE[band_type]

        try:
            acqtype = ACQUISITION_TYPE[spacecraft + '_' + sensor]
        except KeyError:
            acqtype = LandsatAcquisition

        acquisitions.append(acqtype(new))

    return sorted(acquisitions)

