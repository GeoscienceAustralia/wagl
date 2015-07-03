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

def fixname(s):
    """Fix satellite name. Performs 'Landsat7' to 'LANDSAT_7' but also
    handles 'LANDSAT_8' to 'LANDSAT_8'"""
    return re.sub(r'([a-zA-Z]+)_?(\d)',
                  lambda m: m.group(1).upper() + '_' + m.group(2), s)

@total_ordering
class Acquisition(object):

    """Acquisition metadata."""

    def __init__(self, metadata):
        for v in metadata.values():
            self.__dict__.update(v)

    def __eq__(self, other):
        return self.band_name == other.band_name

    def __lt__(self, other):
        return self.sortkey() < other.sortkey()

    def __repr__(self):
        return 'Acquisition(band_name=' + self.band_name + ')'

    def sortkey(self):
        """Representation used for sorting objects."""
        if isinstance(self.band_num, int):
            return "%03d" % (self.band_num, )
        else:
            return self.band_name.replace('band', '')


    def data(self, out=None, window=None, masked=False):
        """
        Return `numpy.array` of the data for this acquisition.
        If `out` is supplied, it must be a numpy.array into which
        the Acquisition's data will be read.
        """
        return gaip.data(self, out=out, window=window, masked=masked)

    def data_and_box(self, out=None, window=None, masked=False):
        """
        Return a tuple comprising the `numpy.array` of the data for this
        Acquisition and the `GriddedGeoBox` describing the spatial extent.
        If `out` is supplied, it must be a numpy.array into which
        the Acquisition's data will be read.
        for this acquisition.
        """
        return gaip.data_and_box(self, out=out, window=window, masked=masked)

    def gridded_geo_box(self):
        """Return the `GriddedGeoBox` for this acquisition."""
        return gaip.gridded_geo_box(self)

    @property
    def no_data(self):
        """
        Return the no_data value for this acquisition.
        """
        return gaip.no_data(self)


class LandsatAcquisition(Acquisition):

    """A Landsat acquisition."""

    def __init__(self, metadata):
        super(LandsatAcquisition, self).__init__(metadata)

    @property
    def samples(self):
        """The number of samples (aka. `width`)."""
        if self.band_type == REF:
            return self.product_samples_ref
        if self.band_type == THM:
            return self.product_samples_thm
        if self.band_type == PAN:
            return self.product_samples_pan

    @property
    def lines(self):
        """The number of lines (aka. `height`)."""
        if self.band_type == REF:
            return self.product_lines_ref
        if self.band_type == THM:
            return self.product_lines_thm
        if self.band_type == PAN:
            return self.product_lines_pan

    @property
    def grid_cell_size(self):
        """The resolution of the cell."""
        if self.band_type == REF:
            return self.grid_cell_size_ref
        if self.band_type == THM:
            return self.grid_cell_size_thm
        if self.band_type == PAN:
            return self.grid_cell_size_pan

    @property
    def height(self):
        """The height of the acquisition (aka. `lines`)."""
        return self.lines

    @property
    def width(self):
        """The width of the acquisition (aka. `samples`)."""
        return self.samples

    @property
    def resolution(self):
        """The resolution of the acquisition (aka. `grid_cell_size`)."""
        return self.grid_cell_size

    @property
    def scene_center_date(self):
        """The acquisition date."""
        return self.acquisition_date

    @property
    def scene_centre_date(self):
        """The acquisition date."""
        return self.acquisition_date

    @property
    def scene_center_datetime(self):
        """The acquisition time."""
        return datetime.datetime.combine(self.acquisition_date,
                                         self.scene_center_time)

    @property
    def scene_centre_datetime(self):
        """The acquisition time."""
        return self.scene_center_datetime

    @property
    def scene_centre_time(self):
        """The acquisition time."""
        return self.scene_center_time

    @property
    def min_radiance(self):
        """The minimum radiance (aka. `lmin`)."""
        return self.lmin

    @property
    def max_radiance(self):
        """The maximum radiance (aka. `lmax`)."""
        return self.lmax

    @property
    def decimal_hour(self):
        """The time in decimal."""
        time = self.scene_centre_time
        return (time.hour + (time.minute + (time.second
                                            + time.microsecond / 1000000.0)
                             / 60.0) / 60.0)

    @property
    def gain(self):
        """The sensor gain"""
        return (self.lmax - self.lmin)/(self.qcalmax - self.qcalmin)

    @property
    def bias(self):
        """Sensor bias"""
        return self.lmax - (self.gain * self.qcalmax)

    @property
    def wavelength(self):
        return gaip.SENSORS[self.spacecraft_id]['sensors'][self.sensor_id]['bands']\
            [str(self.band_num)]['wavelength']

    @property
    def band_desc(self):
        return gaip.SENSORS[self.spacecraft_id]['sensors'][self.sensor_id]['bands']\
            [str(self.band_num)]['desc']

    @property
    def resolution(self):
        return gaip.SENSORS[self.spacecraft_id]['sensors'][self.sensor_id]['bands']\
            [str(self.band_num)]['resolution']

    @property
    def band_type_desc(self):
        return gaip.SENSORS[self.spacecraft_id]['sensors'][self.sensor_id]['bands']\
            [str(self.band_num)]['type_desc']

class Landsat5Acquisition(LandsatAcquisition):

    """ Landsat 5 acquisition. """

    def __init__(self, metadata):
        super(Landsat5Acquisition, self).__init__(metadata)

    @property
    def scene_center_time(self):
        """The acquisition time."""
        return self.scene_center_scan_time

    @property
    def date_acquired(self):
        """The acquisition time."""
        return self.acquisition_date

    @property
    def path(self):
        """The acquisitions path."""
        return self.wrs_path

    @property
    def row(self):
        """The acquisition row."""
        return self.starting_row

class Landsat7Acquisition(LandsatAcquisition):

    """ Landsat 7 acquisition. """

    def __init__(self, metadata):
        super(Landsat7Acquisition, self).__init__(metadata)

    def sortkey(self):
        return self.band_name.replace('band', '')

    @property
    def scene_center_time(self):
        """The acquisition time."""
        return self.scene_center_scan_time

    @property
    def date_acquired(self):
        """The acquisition time."""
        return self.acquisition_date

    @property
    def path(self):
        """The acquisitions path."""
        return self.wrs_path

    @property
    def row(self):
        """The acquisition row."""
        return self.starting_row

class Landsat8Acquisition(LandsatAcquisition):

    """ Landsat 8 acquisition. """

    def __init__(self, metadata):
        super(Landsat8Acquisition, self).__init__(metadata)

    @property
    def samples(self):
        """The number of samples (aka. `width`)."""
        if self.band_type == REF:
            return self.reflective_samples
        if self.band_type == ATM:
            return self.reflective_samples
        if self.band_type == BQA:
            return self.reflective_samples
        if self.band_type == PAN:
            return self.panchromatic_samples
        if self.band_type == THM:
            return self.thermal_samples

    @property
    def lines(self):
        """The number of lines (aka. `height`)."""
        if self.band_type == REF:
            return self.reflective_lines
        if self.band_type == ATM:
            return self.reflective_lines
        if self.band_type == BQA:
            return self.reflective_lines
        if self.band_type == PAN:
            return self.panchromatic_lines
        if self.band_type == THM:
            return self.thermal_lines

    @property
    def grid_cell_size(self):
        """The resolution of the cell."""
        if self.band_type == REF:
            return self.grid_cell_size_reflective
        if self.band_type == ATM:
            return self.grid_cell_size_reflective
        if self.band_type == BQA:
            return self.grid_cell_size_reflective
        if self.band_type == PAN:
            return self.grid_cell_size_panchromatic
        if self.band_type == THM:
            return self.grid_cell_size_thermal

    @property
    def acquisition_date(self):
        """The acquisition time."""
        return self.date_acquired

    @property
    def min_radiance(self):
        """The minimum radiance."""
        return getattr(self, 'radiance_minimum')

    @property
    def max_radiance(self):
        """The minimum radiance."""
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
        """The minimum reflectance."""
        return getattr(self, 'reflectance_minimum')

    @property
    def max_reflectance(self):
        """The maximum reflectance."""
        return getattr(self, 'reflectance_maximum')

    @property
    def zone_number(self):
        """The UTM zone number."""
        return getattr(self, 'utm_zone')

    @property
    def gain(self):
        """The sensor gain
        """
        return self.radiance_mult

    @property
    def bias(self):
        """Sensor bias
        Use value from MTL file for consistency with SceneDataset code
        """
        return self.radiance_add

    @property
    def path(self):
        """The acquisitions path."""
        return self.wrs_path

    @property
    def row(self):
        """The acquisition row."""
        return self.wrs_row


ACQUISITION_TYPE = {
    'Landsat5_TM': Landsat5Acquisition,
    'Landsat7_ETM+': Landsat7Acquisition,
    'LANDSAT_5_TM': Landsat5Acquisition,
    'LANDSAT_7_ETM+': Landsat7Acquisition,
    'LANDSAT_8_OLI': Landsat8Acquisition,
    'LANDSAT_8_OLI_TIRS': Landsat8Acquisition
}


def find_in(path, s, suffix='txt'):
    """Search through `path` and its children for the first occurance of a
    file with `s` in its name. Returns the path of the file or `None`. """
    for root, _, files in os.walk(path):
        for f in files:
            if s in f and f.endswith(suffix):
                return os.path.join(root, f)
    return None

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


def acquisitions(path):
    """
    Return a list of Acquisition objects from `path`. The argument `path`
    can be a MTL file or a directory name. 

    If `path` is a directory will be searched to find an MTL file. If this 
    search fails the directory will be searched for a collection of GeoTiff 
    files.

    The search will include the `path` directory and its children.
    """

    try:
        acqs = acquisitions_via_mtl(path)
    except OSError:
        acqs = acquisitions_via_geotiff(path)

    return acqs

def acquisitions_via_geotiff(path):
    """
    Collect all the GeoTiffs in the supplied directory path and return as
    a list of Acquisitions. Acquisition properties are extracted from the
    filename.
    """
    name_pattern = r'(?P<spacecraft_id>LS\d)_(?P<sensor_id>\w+)_' \
                   r'(?P<product_type>\w+)_(?P<product_id>P\d+)_' \
                   r'GA(?P<product_code>.*)-(?P<station_id>\d+)_' \
                   r'(?P<wrs_path>\d+)_(?P<wrs_row>\d+)_(?P<acqu' \
                   r'isition_date>\d{8})_B(?P<band_num>\d+)\.tif' 

    acqs = []
    if isdir(path):
        p = re.compile(name_pattern)
        for tif_path in find_all_in(path, 'tif'):
    
            new = {}
 
            dir_name, file_name = os.path.split(tif_path)
            match_obj = p.match(file_name)
            if match_obj is not None:
                md = match_obj.groupdict()
                new['FILENAME_FIELDS'] = md

                # find the spacecraft based on the tag
                tag = fixname(md['spacecraft_id'])
                for k, v in SENSORS.iteritems():
                    if v['tag'] == tag:
                        md['spacecraft_id'] = fixname(k)
                        break

                # get spacecraft info from SENSORS
                spacecraft = md['spacecraft_id']
                new['SPACECRAFT'] = {}
                db = SENSORS[spacecraft]
                for k, v in db.iteritems():
                    if k is not 'sensors':
                        try:
                            new['SPACECRAFT'][k.encode('ascii')] = v.encode('ascii')
                        except AttributeError:
                            new['SPACECRAFT'][k.encode('ascii')] = v

            
                # map sensor_id for consistency with SENSOR keys
                if md['sensor_id'] == 'ETM':
                    md['sensor_id'] = 'ETM+'

                # get sensor info from SENSORS

                new['SENSOR_INFO'] = {}
                sensor =  md['sensor_id']
                db = db['sensors'][sensor]
                for k, v in db.iteritems():
                    if k is not 'bands':
                        new['SENSOR_INFO'][k.encode('ascii')] = v
 
                # normalise the band number

                bn = int(md['band_num'])
                if bn % 10 == 0:
                    bn = bn / 10
                md['band_num'] = bn

                # get band info from SENSORS


                bandname = str(bn)
                new['BAND_INFO'] = {}
                db = db['bands'][bandname]
                for k, v in db.iteritems():
                    new['BAND_INFO'][k.encode('ascii')] = v
                band_type = db['type_desc']
                new['BAND_INFO']['band_type'] = BAND_TYPE[band_type]

                # convert acquisition_date to a datetime

                ad = md['acquisition_date']
                md['acquisition_date'] = datetime.datetime(int(ad[0:4]),
                                                           int(ad[4:6]),
                                                           int(ad[6:8]))

                # band_name is required

                md['band_name'] = 'band%d' % (bn, )

                # file and directory name

                md['dir_name'] = dir_name
                md['file_name'] = file_name

                try:
                    acqtype = ACQUISITION_TYPE[spacecraft + '_' + sensor]
                except KeyError:
                    acqtype = LandsatAcquisition


                # create the Acquisition

                acqs.append(acqtype(new))

    return sorted(acqs)

def acquisitions_via_mtl(path):
    """Obtain a list of Acquisition objects from `path`. The argument `path`
    can be a MTL file or a directory name. If `path` is a directory then the 
    MTL file will be search for in the directory and its children."""

    if isdir(path):
        filename = find_in(path, 'MTL')
    else:
        filename = path

    if filename is None:
        raise OSError("Cannot find MTL file in %s" % path)


    data = gaip.load_mtl(filename)
    bandfiles = [k for k in data['PRODUCT_METADATA'].keys() if 'band' in k
                 and 'file_name' in k]
    bands_ = [b.replace('file_name', '').strip('_') for b in bandfiles]

    # The new MTL version for LS7 has 'vcid' in some sections
    # So the following is account for and remove such instances
    bands = []
    for band in bands_:
        if 'vcid' in band:
            band = band.replace('_vcid_', '')
        bands.append(band)

    # We now create an acquisition object for each band and make the
    # parameters names nice.

    acqs = []
    for band in bands:
        bandparts = set(band.split('_'))
        # create a new copy
        new = copy.deepcopy(data)

        # remove unnecessary values
        for kv in new.values():
            for k in kv.keys():
                if 'vcid' in k:
                    nk = k.replace('_vcid_', '')
                    kv[nk] = kv.pop(k)
                else:
                    nk = k
                if bandparts.issubset(set(nk.split('_'))): 
                    # remove the values for the other bands
                    rm = [nk.replace(band, b) for b in bands if b != band]
                    for r in rm:
                        try:
                            del kv[r]
                        except KeyError:
                            pass
                    # rename old key to remove band information
                    newkey = nk.replace(band, '').strip('_')
                    # print "band=%s, k=%s, newkey=%s" % (band, k, newkey)
                    kv[newkey] = kv[nk]
                    del kv[nk]

        # set path
        dir_name = os.path.dirname(os.path.abspath(filename))
        new['PRODUCT_METADATA']['dir_name'] = dir_name

        # set band name and number
        new['PRODUCT_METADATA']['band_name'] = band
        try:
            band_num = int(band.replace('band', '').strip('_'))
            new['PRODUCT_METADATA']['band_num'] = band_num
        except ValueError:
            new['PRODUCT_METADATA']['band_num'] = None

        product = new['PRODUCT_METADATA']
        spacecraft = fixname(product['spacecraft_id'])
        sensor = product['sensor_id']

        # Account for a change in the new MTL files
        if 'acquisition_date' not in product:
            product['acquisition_date'] = product['date_acquired']
        if 'scene_center_scan_time' not in product:
            product['scene_center_scan_time'] = product['scene_center_time']
        if  'product_samples_ref' not in product:
            product['product_samples_ref'] = product['reflective_samples']
        if  'product_lines_ref' not in product:
            product['product_lines_ref'] = product['reflective_lines']
        if  'product_samples_thm' not in product:
            product['product_samples_thm'] = product['thermal_samples']
        if  'product_lines_thm' not in product:
            product['product_lines_thm'] = product['thermal_lines']

        new['SPACECRAFT'] = {}
        db = SENSORS[spacecraft]
        for k, v in db.iteritems():
            if k is not 'sensors':
                try:
                    new['SPACECRAFT'][k.encode('ascii')] = v.encode('ascii')
                except AttributeError:
                    new['SPACECRAFT'][k.encode('ascii')] = v

        new['SENSOR_INFO'] = {}
        try:
            db = db['sensors'][sensor]
        except KeyError, e:
            msg = 'No Match for {} sensor; Trying ETM+!'.format(e.message)
            print msg
            db = db['sensors']['ETM+']

        for k, v in db.iteritems():
            if k is not 'bands':
                new['SENSOR_INFO'][k.encode('ascii')] = v

        bandname = band.replace('band', '').strip('_')
        new['BAND_INFO'] = {}
        db = db['bands'][bandname]

        for k, v in db.iteritems():
            new['BAND_INFO'][k.encode('ascii')] = v
        band_type = db['type_desc']
        new['BAND_INFO']['band_type'] = BAND_TYPE[band_type]

        try:
            acqtype = ACQUISITION_TYPE[spacecraft + '_' + sensor]
        except KeyError:
            acqtype = LandsatAcquisition

        acqs.append(acqtype(new))

    return sorted(acqs)
