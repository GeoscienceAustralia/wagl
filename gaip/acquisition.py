"""
Core code.
"""
from __future__ import absolute_import, print_function
import os
from os.path import isdir, join as pjoin, dirname, basename, exists
import re
import copy
import json
import datetime
import glob
import pandas
from pkg_resources import resource_stream

from functools import total_ordering
import rasterio
from gaip.data import data, data_and_box, no_data
from gaip.geobox import GriddedGeoBox
from gaip.modtran import read_spectral_response
from gaip.mtl import load_mtl


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
    handles 'LANDSAT8' to 'LANDSAT_8'"""
    return re.sub(r'([a-zA-Z]+)_?(\d)',
                  lambda m: m.group(1).upper() + '_' + m.group(2), s)


class AcquisitionsContainer(object):

    """
    A container for dealing with a hierarchial structure
    of acquisitions from different groups, granules, but
    all part of the same geospatial area or scene.
    
    Note: Assuming that each granule contains the same groups.

    The `AcquisitionsContainer.tiled` property indicates whether or
    not a scene is partitioned into several tiles referred to as
    granules.
    """

    def __init__(self, groups=None, granules=None):
        self._tiled = False if granules is None else True
        self._groups = groups
        self._granules = granules

    def __repr__(self):
        fmt = ("****Tiled scene****:\n{tiled}\n"
               "****Granules****:\n{granules}\n"
               "****Groups****:\n{groups}")
        granules = "\n".join(self.granules)
        groups = "\n".join(self.groups)
        return fmt.format(tiled=self.tiled, granules=granules, groups=groups)

    @property
    def tiled(self):
        """
        Indicates whether or not a scene is partitioned into several
        tiles referred to as granules.
        """
        return self._tiled

    @property
    def granules(self):
        """
        Lists the available granules within a scene.
        If `AcquisitionsContainer.tiled` is False, then [None] is
        returned.
        """
        return list(self._granules.keys()) if self.tiled else [None]

    @property
    def groups(self):
        """
        Lists the available groups within a scene.
        """
        if self.tiled:
            grps = list(self._granules.get(self.granules[0]).keys())
        else:
            grps = list(self._groups.keys())
        return grps

    def get_acquisitions(self, group=None, granule=None):
        """
        Return a list of acquisitions for a given granule and group.

        :param group:
            A `str` defining the group layer from which to retrieve
            the acquisitions from. If `None` (default), return the
            acquisitions from the first group in the
            `AcquisitionsContainer.groups` list.

        :param granule:
            A `str` defining the granule layer from which to retrieve
            the acquisitions from. If `None` (default), return the
            acquisitions from the the first granule in the
            `AcquisitionsContainer.granule` list.

        :return:
            A `list` of `Acquisition` objects.
        """
        if self.tiled:
            groups = self.get_granule(granule=granule)
            if group is None:
                acqs = groups[list(groups.keys())[0]]
            else:
                acqs = groups[group]
        else:
            if group is None:
                acqs = self._groups[self.groups[0]]
            else:
                acqs = self._groups[group]

        return acqs

    def get_granule(self, granule=None):
        """
        Return a granule containing groups of `Acquisition` objects.

        :param granule:
            A `str` defining the granule layer from which to retrieve
            groups of `Acquisition` objects. Default is `None`, which
            returns the the first granule in the
            `AcquisitionsContainer.granule` list.

        :return:
            A `dict` containing the groups of `Acquisition` objects
            for a given scene.
        """
        if not self.tiled:
            return self._groups
        if granule is None:
            return self._granules[self.granules[0]]
        else:
            return self._granules[granule]

    def get_root(self, path='', group=None, granule=None):
        """
        Get the root level file system path for a granule and/or group
        within the `AcquisitionContainer` object.

        :param path:
            A `str` containing the root path on which to join the
            granule and/or group layers onto.

        :param group:
            A `str` containing the group layer to be joined onto
            `path`. If group is `None` (default), or `not in`
             `AcquisitionContainer.groups` then no path join occurs.

        :param granule:
            A `str` containing the granule layer to be joined onto
            `path`. If granule is `None` (default), or `not in`
            `AcquisitionContainer.granules` then no path join occurs.

        :return:
            A `str` representing the combined path for group and/or
            granule layers.
        """
        if (granule is None) or (granule not in self.granules):
            root = path
        else:
            root = pjoin(path, granule)

        if (group is not None) and (group in self.groups):
            root = pjoin(root, group)

        return root


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


    def data(self, out=None, window=None, masked=False,
             apply_gain_offset=False, out_no_data=-999):
        """
        Return `numpy.array` of the data for this acquisition.
        If `out` is supplied, it must be a numpy.array into which
        the Acquisition's data will be read.
        """
        return data(self, out=out, window=window, masked=masked,
                    apply_gain_offset=apply_gain_offset,
                    out_no_data=out_no_data)

    def data_and_box(self, out=None, window=None, masked=False):
        """
        Return a tuple comprising the `numpy.array` of the data for this
        Acquisition and the `GriddedGeoBox` describing the spatial extent.
        If `out` is supplied, it must be a numpy.array into which
        the Acquisition's data will be read.
        for this acquisition.
        """
        return data_and_box(self, out=out, window=window, masked=masked)

    def gridded_geo_box(self):
        """Return the `GriddedGeoBox` for this acquisition."""
        with rasterio.open(pjoin(self.dir_name, self.file_name), 'r') as src:
            return GriddedGeoBox.from_dataset(src)

    @property
    def no_data(self):
        """
        Return the no_data value for this acquisition.
        """
        return no_data(self)

    @property
    def gps_file(self):
        return False

    @property
    def scaled_radiance(self):
        """
        Do we have a scaled "at sensor radiance" unit?
        If `True`, then this property needs to be overridden, and define
        the bias and gain properties.
        """
        return False

    def spectral_response(self, as_list=False):
        """
        Reads the spectral response for the sensor.
        """
        fname = 'spectral_response/%s' % self.spectral_filter_file
        spectral_range = range(*self.spectral_range)
        with resource_stream(__name__, fname) as src:
            df = read_spectral_response(src, as_list, spectral_range)
        return df


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
        if self.band_type == BQA:
            return self.product_samples_ref

    @property
    def lines(self):
        """The number of lines (aka. `height`)."""
        if self.band_type == REF:
            return self.product_lines_ref
        if self.band_type == THM:
            return self.product_lines_thm
        if self.band_type == PAN:
            return self.product_lines_pan
        if self.band_type == BQA:
            return self.product_lines_ref

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
        try:
            lmin = self.lmin
        except AttributeError:
            lmin = self.radiance_minimum
        return lmin

    @property
    def max_radiance(self):
        """The maximum radiance (aka. `lmax`)."""
        try:
            lmax = self.lmax
        except AttributeError:
            lmax = self.radiance_maximum
        return lmax

    @property
    def min_quantize(self):
        """THe minimum quantize calibration (aka. `qcal_min`)."""
        try:
            qcal_min = self.qcalmin
        except AttributeError:
            qcal_min = self.quantize_cal_min
        return qcal_min

    @property
    def max_quantize(self):
        """THe maximum quantize calibration (aka. `qcal_max`)."""
        try:
            qcal_max = self.qcalmax
        except AttributeError:
            qcal_max = self.quantize_cal_max
        return qcal_max

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
        return ((self.max_radiance - self.min_radiance) /
                (self.max_quantize - self.min_quantize))

    @property
    def bias(self):
        """Sensor bias"""
        return self.max_radiance - (self.gain * self.max_quantize)

    @property
    def wavelength(self):
        return SENSORS[self.spacecraft_id]['sensors'][self.sensor_id]['bands']\
            [str(self.band_num)]['wavelength']

    @property
    def band_desc(self):
        return SENSORS[self.spacecraft_id]['sensors'][self.sensor_id]['bands']\
            [str(self.band_num)]['desc']

    @property
    def resolution(self):
        return SENSORS[self.spacecraft_id]['sensors'][self.sensor_id]['bands']\
            [str(self.band_num)]['resolution']

    @property
    def band_type_desc(self):
        return SENSORS[self.spacecraft_id]['sensors'][self.sensor_id]['bands']\
            [str(self.band_num)]['type_desc']

    @property
    def scaled_radiance(self):
        """
        Do we have a scaled "at sensor radiance" unit?
        """
        return True


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
        return self.wrs_row


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
        return self.wrs_row


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

    # try:
    #     acqs = acquisitions_via_mtl(path)
    # except OSError:
    #     acqs = acquisitions_via_geotiff(path)

    try:
        acqs = acquisitions_via_safe(path)
    except IOError:
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
                #tag = fixname(md['spacecraft_id'])
                tag = md['spacecraft_id']
                for k, v in SENSORS.items():
                    if v['tag'] == tag:
                        md['spacecraft_id'] = fixname(k)
                        #md['spacecraft_id'] = k
                        break

                # get spacecraft info from SENSORS
                spacecraft = md['spacecraft_id']
                new['SPACECRAFT'] = {}
                db = SENSORS[spacecraft]
                for k, v in db.items():
                    if k is not 'sensors':
                        try:
                            new['SPACECRAFT'][k] = v
                        except AttributeError:
                            new['SPACECRAFT'][k] = v

            
                # map sensor_id for consistency with SENSOR keys
                if md['sensor_id'] == 'ETM':
                    md['sensor_id'] = 'ETM+'

                if md['sensor_id'] == 'OLITIRS':
                    md['sensor_id'] = 'OLI_TIRS'

                # get sensor info from SENSORS

                new['SENSOR_INFO'] = {}
                sensor = md['sensor_id']
                db = db['sensors'][sensor]
                for k, v in db.items():
                    if k is not 'bands':
                        new['SENSOR_INFO'][k] = v
 
                # normalise the band number

                bn = int(md['band_num'])
                if bn % 10 == 0:
                    bn = bn / 10
                md['band_num'] = bn

                # get band info from SENSORS


                bandname = str(bn)
                new['BAND_INFO'] = {}
                db = db['bands'][bandname]
                for k, v in db.items():
                    new['BAND_INFO'][k] = v
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

    return AcquisitionsContainer(groups={'product': sorted(acqs)})


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


    data = load_mtl(filename)
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
        for kv in list(new.values()):
            for k in list(kv.keys()):
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
        band_num = band.replace('band', '').strip('_')
        new['PRODUCT_METADATA']['band_num'] = band_num

        product = new['PRODUCT_METADATA']
        spacecraft = fixname(product['spacecraft_id'])
        product['spacecraft_id'] = spacecraft
        if product['sensor_id'] == 'ETM':
            product['sensor_id'] = 'ETM+'
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

        # Account for missing thermal bands in OLI only products
        if product['sensor_id'] != 'OLI':
            if  'product_samples_thm' not in product:
                product['product_samples_thm'] = product['thermal_samples']
            if  'product_lines_thm' not in product:
                product['product_lines_thm'] = product['thermal_lines']

        new['SPACECRAFT'] = {}
        db = SENSORS[spacecraft]
        for k, v in db.items():
            if k is not 'sensors':
                try:
                    new['SPACECRAFT'][k] = v
                except AttributeError:
                    new['SPACECRAFT'][k] = v

        new['SENSOR_INFO'] = {}
        db = db['sensors'][sensor]

        for k, v in db.items():
            if k is not 'bands':
                new['SENSOR_INFO'][k] = v

        bandname = band.replace('band', '').strip('_')
        new['BAND_INFO'] = {}
        db = db['bands'][bandname]

        for k, v in db.items():
            new['BAND_INFO'][k] = v
        band_type = db['type_desc']
        new['BAND_INFO']['band_type'] = BAND_TYPE[band_type]
  
        try:
            acqtype = ACQUISITION_TYPE[spacecraft + '_' + sensor]
        except KeyError:
            acqtype = LandsatAcquisition

        acqs.append(acqtype(new))

    return AcquisitionsContainer(groups={'product': sorted(acqs)})


def acquisitions_via_safe(path):
    """
    Collect the TOA Radiance images for each granule within a scene.
    Multi-granule & multi-resolution hierarchy format.
    Returns a dict of granules, each granule contains a dict of resolutions,
    and each resolution contains a list of acquisition objects.

    Example:
    {'GRANULE_1': {'R10m': [`acquisition_1`,...,`acquisition_n`],
                   'R20m': [`acquisition_1`,...,`acquisition_n`],
                   'R60m': [`acquisition_1`,...,`acquisition_n`]},
                .
                .
                .
     'GRANULE_N': {'R10m': [`acquisition_1`,...,`acquisition_n`],
                   'R20m': [`acquisition_1`,...,`acquisition_n`],
                   'R60m': [`acquisition_1`,...,`acquisition_n`]}}
    """
    granules = {}
    spacecraft = "SENTINEL-2A"
    sensor = "MSI"

    gps_fname = pjoin(path, "GPS_points")

    granule_dir = pjoin(path, 'GRANULE')
    if not exists(granule_dir):
        raise IOError("GRANULE directory not found: {}".format(granule_dir))

    res_dirs = ["R10m", "R20m", "R60m"]
    for granule in os.listdir(granule_dir):
        resolutions = {}

        img_dir = pjoin(pjoin(granule_dir, granule), 'IMG_DATA')
        if not exists(img_dir):
            continue

        for res_dir in res_dirs:
            acqs = []
            data_dir = pjoin(img_dir, res_dir)
            cwd = os.getcwd()
            os.chdir(data_dir)
            bands = glob.glob('L_output*.hdr')
            os.chdir(cwd)
            bands = [pjoin(data_dir, b.replace('.hdr', '')) for b in bands]

            tile_metadata_fname = 'tile_metadata_{}'.format(res_dir[1:3])
            with open(pjoin(data_dir, tile_metadata_fname), 'r') as src:
                md = src.readlines()

            metadata = {}
            for m in md:
                k, v = m.split('=')
                metadata[k.strip()] = v.strip()

            metadata['EPSG'] = int(metadata['EPSG'])

            dt = pandas.to_datetime(metadata['SENSING_TIME']).to_pydatetime()
            metadata['SENSING_TIME'] = dt

            clon, clat = metadata['centre_lon_lat'].split()
            metadata['centre_lon_lat'] = (float(clon), float(clat))

            ll_lon, ll_lat = metadata['ll_lon_and_ll_lat'].split()
            metadata['ll_lon_and_ll_lat'] = (float(ll_lon), float(ll_lat))

            lr_lon, lr_lat = metadata['lr_lon_and_lr_lat'].split()
            metadata['lr_lon_and_lr_lat'] = (float(lr_lon), float(lr_lat))

            ul_lon, ul_lat = metadata['ul_lon_and_ul_lat'].split()
            metadata['ul_lon_and_ul_lat'] = (float(ul_lon), float(ul_lat))

            ur_lon, ur_lat = metadata['ur_lon_and_ur_lat'].split()
            metadata['ur_lon_and_ur_lat'] = (float(ur_lon), float(ur_lat))

            metadata['ncols'] = int(metadata['ncols'])
            metadata['nrows'] = int(metadata['nrows'])

            metadata['ulx'] = float(metadata['ulx'])
            metadata['uly'] = float(metadata['uly'])

            metadata['resolution'] = float(metadata['resolution'])

            metadata['GPS_Filename'] = gps_fname

            metadata['sensor_id'] = "MSI"

            data = {}
            data['PRODUCT_METADATA'] = copy.deepcopy(metadata)

            #metadata['SPACECRAFT'] = {}
            data['SPACECRAFT'] = {}
            db = SENSORS[spacecraft]
            for k, v in db.items():
                if k is not 'sensors':
                    try:
                        data['SPACECRAFT'][k] = v
                    except AttributeError:
                        data['SPACECRAFT'][k] = v

            data['SENSOR_INFO'] = {}
            db = db['sensors'][sensor]

            for k, v in db.items():
                if k is not 'bands':
                    data['SENSOR_INFO'][k] = v

            for band in bands:
                dname = dirname(band)
                fname = basename(band)
                md = copy.deepcopy(data)
                band_md = md['PRODUCT_METADATA']
                band_md['dir_name'] = dname
                band_md['file_name'] = fname
                bnum = fname.split('_')[2]
                band_name = bnum

                md['BAND_INFO'] = {}
                db_copy = copy.deepcopy(db)
                db_copy = db_copy['bands'][band_name]

                for k, v in db_copy.items():
                    md['BAND_INFO'][k] = v

                band_type = db_copy['type_desc']
                md['BAND_INFO']['band_type'] = BAND_TYPE[band_type]
                band_md['band_type'] = BAND_TYPE[band_type]

                band_md['band_name'] = 'band_{}'.format(bnum)
                band_md['band_num'] = bnum

                acqs.append(Sentinel2aAcquisition(md))

            resolutions[res_dir] = sorted(acqs)
        granules[granule] = resolutions

    return AcquisitionsContainer(granules=granules)


class Sentinel2aAcquisition(Acquisition):

    def __init__(self, metadata):
        super(Sentinel2aAcquisition, self).__init__(metadata)

    @property
    def samples(self):
        """The number of samples (aka. `width`)."""
        return self.ncols

    @property
    def lines(self):
        """The number of lines (aka. `height`)."""
        return self.nrows

    @property
    def width(self):
        """The width of the acquisition (aka. `samples`)."""
        return self.ncols

    @property
    def height(self):
        """The height of the acquisition (aka. `lines`)."""
        return self.nrows

    @property
    def grid_cell_size(self):
        """The resolution of the cell."""
        return self.resolution

    @property
    def acquisition_date(self):
        """The acquisition time."""
        return self.SENSING_TIME

    @property
    def scene_centre_datetime(self):
        # TODO: check that the tile centre times are different
        return self.SENSING_TIME

    @property
    def scene_center_datetime(self):
        return self.SENSING_TIME

    @property
    def scene_centre_date(self):
        return self.SENSING_TIME

    @property
    def decimal_hour(self):
        """The time in decimal."""
        time = self.SENSING_TIME
        return (time.hour + (time.minute + (time.second
                                            + time.microsecond / 1000000.0)
                             / 60.0) / 60.0)

    # TODO: update nbar_constants for new sensor
    @property
    def spacecraft_id(self):
        return "Sentinel2A"

    @property
    def satellite_name(self):
        return "Sentinel2A"

    # @property
    # def sensor_id(self):
    #     return "MSI"

    @property
    def gps_file(self):
        return True

        #geobox = self.gridded_geo_box()
        #ymin = numpy.min([geobox.ul_lonlat[1], geobox.ur_lonlat[1],
        #                  geobox.lr_lonlat[1], geobox.ll_lonlat[1]])
        #ymax= numpy.max([geobox.ul_lonlat[1], geobox.ur_lonlat[1],
        #                 geobox.lr_lonlat[1], geobox.ll_lonlat[1]])

        #subs = df[(df.lat >= ymin) & (df.lat <= ymax)]
        #idx = subs.shape[0] // 2 - 1
        #alon0 = subs.iloc[idx].lon
        #alat0 = subs.iloc[idx].lat

    def read_gps_file(self):
        df = pandas.read_csv(self.GPS_Filename, sep=' ', names=['lon', 'lat'],
                             header=None)
        return df

    @property
    def gain(self):
        return 0.01

    @property
    def bias(self):
        return 0.0

    @property
    def scaled_radiance(self):
        """
        Do we have a scaled "at sensor radiance" unit?
        """
        return True
