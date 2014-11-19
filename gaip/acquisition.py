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
        return getattr(self, 'lmin_' + self.band_name)

    @property
    def max_radiance(self):
        return getattr(self, 'lmax_' + self.band_name)

    @property
    def min_reflectance(self):
        return getattr(self, 'qcalmin_' + self.band_name)

    @property
    def max_reflectance(self):
        return getattr(self, 'qcalmax_' + self.band_name)


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
        return getattr(self, 'radiance_minimum_' + self.band_name)

    @property
    def max_radiance(self):
        return getattr(self, 'radiance_maximum_' + self.band_name)

    @property
    def min_reflectance(self):
        return getattr(self, 'reflectance_minimum_' + self.band_name)

    @property
    def max_reflectance(self):
        return getattr(self, 'reflectance_maximum_' + self.band_name)


ACQUISITION_TYPE = {
    'Landsat5_TM': Landsat5Acquisition,
    'Landsat7_ETM+': Landsat7Acquisition,
    'LANDSAT_8_OLI': Landsat8Acquisition
}


def parse_type(s):
    """Parse the string `s` and return a native python object."""

    strptime = datetime.datetime.strptime

    def yesno(s):
        if len(s) == 1:
            if s == 'Y':
                return True
            if s == 'N':
                return False
        raise ValueError

    def none(s):
        if len(s) == 4 and s == 'NONE':
            return None
        raise ValueError

    parsers = [int,
               float,
               lambda x: strptime(x, '%Y-%m-%dT%H:%M:%SZ'),
               lambda x: strptime(x, '%Y-%m-%d').date(),
               lambda x: strptime(x[0:15], '%H:%M:%S.%f').time(),
               lambda x: yesno(x.strip('"')),
               lambda x: none(x.strip('"')),
               lambda x: str(x.strip('"'))]

    for parser in parsers:
        try:
            return parser(s)
        except ValueError:
            pass
    raise ValueError


def load_mtl(filename, root='L1_METADATA_FILE', pairs=r'(\w+)\s=\s(.*)'):
    """Parse an MTL file and return dict-of-dict's containing the metadata."""

    def parse(lines, tree, level=0):
        while lines:
            line = lines.pop(0)
            match = re.findall(pairs, line)
            if match:
                key, value = match[0]
                if key == 'GROUP':
                    tree[value] = {}
                    parse(lines, tree[value], level + 1)
                elif key == 'END_GROUP':
                    break
                else:
                    tree[key.lower()] = parse_type(value)

    tree = {}
    with open(filename, 'r') as fo:
        parse(fo.readlines(), tree)

    return tree[root]


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

    data = load_mtl(filename)
    bandfiles = [k for k in data['PRODUCT_METADATA'].keys() if 'band' in k
                 and 'file_name' in k]

    # We now create an acquisition object for each band and make the
    # parameters names nice.

    acquisitions = []
    for bandfile in bandfiles:
        band = bandfile.replace('file_name', '').strip('_')

        # create a new copy
        new = copy.deepcopy(data)

        # set path
        new['PRODUCT_METADATA']['dir_name'] = dirname

        # set band name and number
        new['PRODUCT_METADATA']['band_name'] = band
        try:
            band_num = int(band.replace('band', '').strip('_'))
            new['PRODUCT_METADATA']['band_num'] = band_num
        except ValueError:
            new['PRODUCT_METADATA']['band_num'] = None

        # replace filename
        product = new['PRODUCT_METADATA']
        filename = product[bandfile]
        for k in product.keys():
            if 'band' in k and 'file_name' in k:
                product.pop(k, None)
            product['file_name'] = filename

        def replace_param(section_name):
            try:
                section = new[section_name]
                for k in section.keys():
                    if k.endswith(band):
                        newkey = k.replace('_' + band, '')
                        section[newkey] = section[k]
                    elif k.startswith(band):
                        newkey = k.replace(band + '_', '')
                        section[newkey] = section[k]
                    section.pop(k, None)
            except KeyError:
                pass

        replace_param('PRODUCT_PARAMETERS')
        replace_param('CORRECTIONS_APPLIED')

        spacecraft = product['spacecraft_id']
        sensor = product['sensor_id']

        new['SPACECRAFT'] = {}
        db = SENSORS[spacecraft]
        for k, v in db.iteritems():
            if k is not 'sensors':
                new['SPACECRAFT'][k] = v

        new['SENSOR_INFO'] = {}
        db = db['sensors'][sensor]
        for k, v in db.iteritems():
            if k is not 'bands':
                new['SENSOR_INFO'][k] = v

        bandname = band.replace('band', '').strip('_')
        new['BAND_INFO'] = {}
        db = db['bands'][bandname]
        for k, v in db.iteritems():
            new['BAND_INFO'][k] = v
        band_type = db['type_desc']
        new['BAND_INFO']['band_type'] = BAND_TYPE[band_type]

        try:
            acqtype = ACQUISITION_TYPE[spacecraft + '_' + sensor]
        except KeyError:
            acqtype = LandsatAcquisition

        acquisitions.append(acqtype(new))

    return sorted(acquisitions)
