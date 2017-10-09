"""
Core code.
"""
from __future__ import absolute_import, print_function
import os
from os.path import isdir, join as pjoin, dirname, basename, exists, splitext
import re
import copy
import json
import datetime
import glob
from xml.etree import ElementTree
from functools import total_ordering
from dateutil import parser
import pandas
from pkg_resources import resource_stream

import rasterio
from gaip.geobox import GriddedGeoBox
from gaip.modtran import read_spectral_response
from gaip.mtl import load_mtl
from gaip.constants import BandType


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

    def __init__(self, label, groups=None, granules=None):
        self._tiled = False if granules is None else True
        self._groups = groups
        self._granules = granules
        self._label = label

    def __repr__(self):
        fmt = ("****Tiled scene****:\n{tiled}\n"
               "****Granules****:\n{granules}\n"
               "****Groups****:\n{groups}")
        granules = "\n".join(self.granules) if self.tiled else ""
        groups = "\n".join(self.groups)
        return fmt.format(tiled=self.tiled, granules=granules, groups=groups)

    @property
    def label(self):
        """
        Return the scene label.
        """
        return self._label

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
        return sorted(list(self._granules.keys())) if self.tiled else [None]

    @property
    def groups(self):
        """
        Lists the available groups within a scene.
        """
        if self.tiled:
            grps = sorted(list(self._granules.get(self.granules[0]).keys()))
        else:
            grps = sorted(list(self._groups.keys()))
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

    def get_root(self, path='/', group=None, granule=None):
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

    # def __init__(self, metadata):
    #     for v in metadata.values():
    #         self.__dict__.update(v)
    def __init__(self, pathname, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        self._pathname = pathname
        self._acquisition_datetime = acquisition_datetime
        self._band_name = band_name
        self._band_id = band_id

        self._gps_file = False
        self._scaled_radiance = False

        self._gain = 1.0
        self._bias = 0.0

        if metadata is not None:
            for key, value in metadata.items():
                setattr(self, key, value)

        self._open()

    def _open(self):
        with rasterio.open(self.pathname) as ds:
            self._samples = ds.width
            self._lines = ds.height

    @property
    def pathname(self):
        return self._pathname

    @property
    def acquisition_datetime(self):
        return self._acquisition_datetime

    @property
    def band_name(self):
        return self._band_name

    @property
    def band_id(self):
        return self._band_id

    @property
    def samples(self):
        """The number of samples (aka. `width`)."""
        return self._samples

    @property
    def lines(self):
        """The number of lines (aka. `height`)."""
        return self._lines

    @property
    def gps_file(self):
        return self._gps_file

    @property
    def scaled_radiance(self):
        """
        Do we have a scaled "at sensor radiance" unit?
        If `True`, then this property needs to be overridden, and define
        the bias and gain properties.
        """
        return self._scaled_radiance

    @property
    def gain(self):
        return self._gain

    @property
    def bias(self):
        return self._bias

    def __eq__(self, other):
        return self.band_name == other.band_name

    def __lt__(self, other):
        return self.sortkey() < other.sortkey()

    def __repr__(self):
        return 'Acquisition(band_name=' + self.band_name + ')'

    def sortkey(self):
        """Representation used for sorting objects."""
        return self.band_name

    def data(self, out=None, window=None, masked=False,
             apply_gain_offset=False, out_no_data=-999):
        """
        Return `numpy.array` of the data for this acquisition.
        If `out` is supplied, it must be a numpy.array into which
        the Acquisition's data will be read.
        """
        with rasterio.open(self.pathname) as ds:
            # convert to at sensor radiance as required
            if apply_gain_offset:
                if out is None:
                    out = ds.read(1, window=window, masked=masked)
                else:
                    ds.read(1, out=out, window=window, masked=masked)

                # check for no data
                no_data = acq.no_data if acq.no_data is not None else 0
                nulls = out == no_data

                # gain & offset; y = mx + b
                data = acq.gain * out + acq.bias

                # set the out_no_data value inplace of the input no data value
                data[nulls] = out_no_data

                return data
            return ds.read(1, out=out, window=window, masked=masked)

    def data_and_box(self, out=None, window=None, masked=False):
        """
        Return a tuple comprising the `numpy.array` of the data for this
        Acquisition and the `GriddedGeoBox` describing the spatial extent.
        If `out` is supplied, it must be a numpy.array into which
        the Acquisition's data will be read.
        for this acquisition.
        """
        return data_and_box(self, out=out, window=window, masked=masked)
        with rasterio.open(self.pathname) as ds:
            box = GriddedGeoBox.from_dataset(ds)
            if window is not None:
                rows = window[0][1] - window[0][0]
                cols = window[1][1] - window[1][0]
                prj = ds.crs.wkt
                res = ds.res

                # Get the new UL co-ordinates of the array
                ul_x, ul_y = ds.transform * (window[1][0], window[0][0])
                box = GriddedGeoBox(shape=(rows, cols), origin=(ul_x, ul_y),
                                    pixelsize=res, crs=prj)
            return (fo.read(1, out=out, window=window, masked=masked), box)

    def gridded_geo_box(self):
        """Return the `GriddedGeoBox` for this acquisition."""
        with rasterio.open(self.pathname) as src:
            return GriddedGeoBox.from_dataset(src)

    def decimal_hour(self):
        """The time in decimal."""
        time = self.acquisition_datetime
        dec_hour = (time.hour + (time.minute + (time.second
                                                + time.microsecond / 1000000.0)
                                 / 60.0) / 60.0)
        return dec_hour

    def julian_day(self):
        """
        Return the Juilan Day of the acquisition_datetime.
        """
        return int(self.acquisition_datetime.strftime('%j'))

    @property
    def no_data(self):
        """
        Return the no_data value for this acquisition.
        Assumes that the acquisition is a single band file.
        """
        with rasterio.open(self.pathname) as ds:
            nodata_list = ds.nodatavals
            return nodata_list[0]

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

    def __init__(self, pathname, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        self.min_radiance = 0
        self.max_radiance = 1
        self.min_quantize = 0
        self.max_quantize = 1

        super(LandsatAcquisition, self).__init__(pathname,
                                                 acquisition_datetime,
                                                 band_name=band_name,
                                                 band_id=band_id,
                                                 metadata=metadata)

        self._gain = ((self.max_radiance - self.min_radiance) /
                      (self.max_quantize - self.min_quantize))
        self._bias = self.max_radiance - (self.gain * self.max_quantize)

        self._scaled_radiance = True


class Landsat5Acquisition(LandsatAcquisition):

    """ Landsat 5 acquisition. """

    def __init__(self, pathname, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Landsat5Acquisition, self).__init__(pathname,
                                                  acquisition_datetime,
                                                  band_name=band_name,
                                                  band_id=band_id,
                                                  metadata=metadata)

        self.platform_id = 'LANDSAT-5'
        self.sensor_id = 'TM'
        self.tle_format = 'l5_%4d%s_norad.txt'
        self.tag = 'LS5'
        self.altitude = 705000.0
        self.inclination = 1.7139133254584316445390643346558
        self.omega = 0.001059
        self.radius = 7285600.0
        self.semi_major_axis = 7083160.0
        self.maximum_view_angle = 9.0


class Landsat7Acquisition(LandsatAcquisition):

    """ Landsat 7 acquisition. """

    def __init__(self, pathname, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Landsat7Acquisition, self).__init__(pathname,
                                                  acquisition_datetime,
                                                  band_name=band_name,
                                                  band_id=band_id,
                                                  metadata=metadata)

        self.platform = 'LANDSAT-7'
        self.sensor_id = 'ETM+'
        self.tle_format = 'L7%4d%sASNNOR.S00'
        self.tag = 'LS7'
        self.altitude = 705000.0
        self.inclination = 1.7139133254584316445390643346558
        self.omega = 0.001059
        self.radius = 7285600.0
        self.semi_major_axis = 7083160.0
        self.maximum_view_angle = 9.0


class Landsat8Acquisition(LandsatAcquisition):

    """ Landsat 8 acquisition. """

    def __init__(self, pathname, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Landsat8Acquisition, self).__init__(pathname,
                                                  acquisition_datetime,
                                                  band_name=band_name,
                                                  band_id=band_id,
                                                  metadata=metadata)

        self.platform_id = 'LANDSAT-8'
        self.sensor_id = 'OLI'
        self.tle_format = 'L8%4d%sASNNOR.S00'
        self.tag = 'LS8'
        self.altitude = 705000.0
        self.inclination = 1.7139133254584316445390643346558
        self.omega = 0.001059
        self.radius = 7285600.0
        self.semi_major_axis = 7083160.0
        self.maximum_view_angle = 9.0


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

def acquisitions2(pathname):
    """
    Return an acquisitions container.
    """
    if '.zip' in splitext(pathname):
        # assume it is the SAFE format from ESA and open it
        container = acquisitions_via_safe2(pathname)

def acquisitions_via_safe2(pathname):
    """
    Read the SAFE format and return an acquisitions container.
    """
    archive = zipfile.ZipFile(pathname)
    xmlfiles = [s for s in archive.namelist() if "MTD_MSIL1C.xml" in s]

    if len(xmlfiles) == 0:
        pattern = pathname.replace('PRD_MSIL1C', 'MTD_SAFL1C')
        pattern = pattern.replace('.zip', '.xml')
        xmlzipfiles = [s for s in archive.namelist() if pattern in s]

    mtd_xml = archive.read(xmlfiles[0])
    xml_root = ElementTree.XML(mtd_xml)

    # what do we do about the 'scene_centre_time' ???
    # DATATAKE_SENSING_START is another potential field to read...
    product_start_time = parser.parse(xml_root.findall('./*/Product_Info/PRODUCT_START_TIME')[0].text)
    product_stop_time = parser.parse(xml_root.findall('./*/Product_Info/PRODUCT_STOP_TIME')[0].text)

    platform = xml_root.findall('./*/Product_Info/*/SPACECRAFT_NAME')[0].text
    # need sensor name


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

    # set path
    dir_name = os.path.dirname(os.path.abspath(filename))

    data = load_mtl(filename)
    bandfiles = [k for k in data['PRODUCT_METADATA'].keys() if 'band' in k
                 and 'file_name' in k]
    bands_ = [b.replace('file_name', '').strip('_') for b in bandfiles]

    # create an acquisition object for each band and attach
    # some appropriate metadata/attributes

    # shortcuts to the requried levels
    prod_md = data['PRODUCT_METADATA']
    rad_md = data['MIN_MAX_RADIANCE']
    quant_md = data['MIN_MAX_PIXEL_VALUE']

    # acquisition datetime
    acq_date = prod_md.get('acquisition_date', prod_md['date_acquired'])
    centre_time = prod_md.get('scene_center_scan_time',
                              prod_md['scene_center_time'])
    acq_datetime = datetime.datetime.combine(acq_date, centre_time)

    # platform and sensor id's
    platform_id = fixname(prod_md['spacecraft_id'])
    sensor_id = prod_md['sensor_id']
    if sensor_id == 'ETM':
        sensor_id = 'ETM+'

    # get the appropriate landsat acquisition class
    try:
        acqtype = ACQUISITION_TYPE['_'.join([platform_id, sensor_id])]
    except KeyError:
        acqtype = LandsatAcquisition

    # solar angles
    if 'PRODUCT_PARAMETERS' in data:
        solar_azimuth = data['PRODUCT_PARAMETERS']['sun_azimuth']
        solar_elevation = data['PRODUCT_PARAMETERS']['sun_elevation']
    else:
        solar_azimuth = data['IMAGE_ATTRIBUTES']['sun_azimuth']
        solar_elevation = data['IMAGE_ATTRIBUTES']['sun_elevation']

    # bands to ignore
    ignore = ['band_quality']

    # supported bands for the given platform & sensor id's
    supported_bands = SENSORS[platform_id][sensor_id]['bands']

    acqs = []
    for band in bands_:
        if band in ignore:
            continue

        # band id
        if 'vcid' in band:
            band_id = band.replace('_vcid_', '').replace('band', '').strip('_')
        else:
            band_id = band.replace('band', '').strip('_')

        if band_id not in supported_bands:
            continue

        # band info stored in sensors.json
        sensor_band_info = supported_bands[band_id]

        # band id name, band filename, band full file pathname
        band_name = 'BAND {}'.format(band_id)
        band_fname = prod_md.get('{}_file_name'.format(band),
                                 prod_md['file_name_{}'.format(band)])
        fname = pjoin(dir_name, band_fname)

        min_rad = rad_md.get('lmin_{}'.format(band),
                             rad_md['radiance_minimum_{}'.format(band)])
        max_rad = rad_md.get('lmax_{}'.format(band),
                             rad_md['radiance_maximum_{}'.format(band)])

        min_quant = quant_md.get('qcalmin_{}'.format(band),
                                 quant_md['quantize_cal_min_{}'.format(band)])
        max_quant = quant_md.get('qcalmax_{}'.format(band),
                                 quant_md['quantize_cal_max_{}'.format(band)])

        # metadata
        attrs = {k: v for k, v in sensor_band_info.items()}
        attrs['solar_azimuth'] = solar_azimuth
        attrs['solar_elevation'] = solar_elevation
        attrs['min_radiance'] = min_rad
        attrs['max_radiance'] = max_rad
        attrs['min_quantize'] = min_quant
        attrs['max_quantize'] = max_quant

        acqs.append(acqtype(fname, acq_datetime, band_name, band_id, attrs))

    return AcquisitionsContainer(label=basename(path),
                                 groups={'product': sorted(acqs)})


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
                md['BAND_INFO']['band_type'] = BandType[band_type]
                band_md['band_type'] = BandType[band_type]

                band_md['band_name'] = 'band_{}'.format(bnum)
                band_md['band_num'] = bnum

                acqs.append(Sentinel2aAcquisition(md))

            resolutions[res_dir] = sorted(acqs)
        granules[granule] = resolutions

    return AcquisitionsContainer(label=basename(path), granules=granules)


class Sentinel2aAcquisition(Acquisition):

    def __init__(self, pathname, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Sentinel2aAcquisition, self).__init__(pathname,
                                                    acquisition_datetime,
                                                    band_name=band_name,
                                                    band_id=band_id,
                                                    metadata=metadata)

        self.platform_id: 'SENTINEL-2A'
        self.sensor_id = 'MSI'
        self.tle_format = 'S2A%4d%sASNNOR.S00'
        self.tag = 'S2A'
        self.altitude: 786000.0
        self.inclination: 1.721243708316808
        self.omega: 0.001039918
        self.semi_major_axis: 7164137.0
        self.maximum_view_angle: 20.0

        self._gps_file = True
        self._gain = 0.01
        self._bias = 0.0

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
