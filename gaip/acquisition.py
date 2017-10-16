"""
The acquisition is a core object passed throughout gaip,
which provides a lightweight `Acquisition` class.
An `Acquisition` instance provides data read methods,
and contains various attributes/metadata regarding the band
as well as the sensor, and sensor's platform.

An `AcquisitionsContainer` provides a container like structure
for handling scenes consisting of multiple Granules/Tiles,
and of differing resolutions.
"""
from __future__ import absolute_import, print_function
import os
from os.path import isdir, join as pjoin, dirname, basename
import re
import json
import datetime
from xml.etree import ElementTree
from functools import total_ordering
import zipfile
from dateutil import parser
import numpy
import numexpr
from scipy import interpolate
import pandas
from pkg_resources import resource_stream
from nested_lookup import nested_lookup
import pyproj

import rasterio
from gaip.geobox import GriddedGeoBox
from gaip.modtran import read_spectral_response
from gaip.mtl import load_mtl
from gaip.constants import BandType


with open(pjoin(dirname(__file__), 'sensors.json')) as fo:
    SENSORS = json.load(fo)

def fixname(s):
    """Fix satellite name.
       Performs 'Landsat7' to 'LANDSAT_7', 'LANDSAT8' to 'LANDSAT_8',
       'Landsat-5' to 'LANDSAT_5'.
    """
    return re.sub(r'([a-zA-Z]+)[_-]?(\d)',
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

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        self._pathname = pathname
        self._uri = uri
        self._acquisition_datetime = acquisition_datetime
        self._band_name = band_name
        self._band_id = band_id

        self._gps_file = False

        if metadata is not None:
            for key, value in metadata.items():
                if key == 'band_type':
                    value = BandType[value]
                setattr(self, key, value)

        self._open()

    def _open(self):
        """
        A private method for opening the dataset and
        retrieving the dimensional information.
        """
        with rasterio.open(self.uri) as ds:
            self._samples = ds.width
            self._lines = ds.height

    @property
    def pathname(self):
        """
        The pathname of the level1 dataset.
        """
        return self._pathname

    @property
    def uri(self):
        """
        The uri of the acquisition.
        """
        return self._uri

    @property
    def acquisition_datetime(self):
        """
        The acquisitions centre scantime.
        """
        return self._acquisition_datetime

    @property
    def band_name(self):
        """
        The band name, which goes by the format `BAND {}`.
        """
        return self._band_name

    @property
    def band_id(self):
        """
        The band id as given in the `sensors.json` file.
        """
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
        """
        Does the acquisition have an associated GPS file?
        """
        return self._gps_file

    def __eq__(self, other):
        return self.band_name == other.band_name

    def __lt__(self, other):
        return self.sortkey() < other.sortkey()

    def __repr__(self):
        return 'Acquisition(band_name=' + self.band_name + ')'

    def sortkey(self):
        """Representation used for sorting objects."""
        return self.band_name

    def data(self, out=None, window=None, masked=False):
        """
        Return `numpy.array` of the data for this acquisition.
        If `out` is supplied, it must be a numpy.array into which
        the Acquisition's data will be read.
        """
        with rasterio.open(self.uri) as ds:
            data = ds.read(1, out=out, window=window, masked=masked)

        return data

    def radiance_data(self, window=None, out_no_data=-999):
        """
        Return the data as radiance in watts/(m^2*micrometre).
        Override with a custom version for a specific sensor.
        """
        raise NotImplementedError

    def data_and_box(self, out=None, window=None, masked=False):
        """
        Return a tuple comprising the `numpy.array` of the data for this
        Acquisition and the `GriddedGeoBox` describing the spatial extent.
        If `out` is supplied, it must be a numpy.array into which
        the Acquisition's data will be read.
        for this acquisition.
        """
        with rasterio.open(self.uri) as ds:
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
        with rasterio.open(self.uri) as src:
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
        with rasterio.open(self.uri) as ds:
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

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
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

    def radiance_data(self, window=None, out_no_data=-999):
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


class Landsat5Acquisition(LandsatAcquisition):

    """ Landsat 5 acquisition. """

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Landsat5Acquisition, self).__init__(pathname, uri,
                                                  acquisition_datetime,
                                                  band_name=band_name,
                                                  band_id=band_id,
                                                  metadata=metadata)

        self.platform_id = 'LANDSAT_5'
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

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Landsat7Acquisition, self).__init__(pathname, uri,
                                                  acquisition_datetime,
                                                  band_name=band_name,
                                                  band_id=band_id,
                                                  metadata=metadata)

        self.platform_id = 'LANDSAT_7'
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

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Landsat8Acquisition, self).__init__(pathname, uri,
                                                  acquisition_datetime,
                                                  band_name=band_name,
                                                  band_id=band_id,
                                                  metadata=metadata)

        self.platform_id = 'LANDSAT_8'
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
    Return an instance of `AcquisitionsContainer` containing the
    acquisitions for each Granule Group (if applicable) and
    each sub Group.
    """
    try:
        container = acquisitions_via_safe(path)
    except zipfile.BadZipFile:
        try:
            container = acquisitions_via_mtl(path)
        except OSError:
            raise IOError("No acquisitions found in: {}".format(path))

    return container


def acquisitions_via_mtl(pathname):
    """Obtain a list of Acquisition objects from `path`. The argument `path`
    can be a MTL file or a directory name. If `path` is a directory then the 
    MTL file will be search for in the directory and its children."""

    if isdir(pathname):
        filename = find_in(pathname, 'MTL')
    else:
        filename = pathname

    if filename is None:
        raise OSError("Cannot find MTL file in %s" % pathname)

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
    solar_azimuth = nested_lookup('sun_azimuth', data)[0]
    solar_elevation = nested_lookup('sun_elevation', data)[0]

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
        band_name = attrs.pop('band_name')

        acqs.append(acqtype(pathname, fname, acq_datetime, band_name, band_id,
                            attrs))

    return AcquisitionsContainer(label=basename(pathname),
                                 groups={'product': sorted(acqs)})


def acquisitions_via_safe(pathname):
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
    def group_helper(fname, resolution_groups):
        """
        A helper function to find the resolution group
        and the band_id, and bail rather than loop over
        everything.
        """
        for key, group in resolution_groups.items():
            for item in group:
                if item in fname:
                    band_id = re.sub(r'B[0]?', '', item)
                    return key, band_id
        return None, None

    archive = zipfile.ZipFile(pathname)
    xmlfiles = [s for s in archive.namelist() if "MTD_MSIL1C.xml" in s]

    if not xmlfiles:
        pattern = basename(pathname.replace('PRD_MSIL1C', 'MTD_SAFL1C'))
        pattern = pattern.replace('.zip', '.xml')
        xmlfiles = [s for s in archive.namelist() if pattern in s]

    mtd_xml = archive.read(xmlfiles[0])
    xml_root = ElementTree.XML(mtd_xml)

    # platform id, TODO: sensor name
    search_term = './*/Product_Info/*/SPACECRAFT_NAME'
    platform_id = fixname(xml_root.findall(search_term)[0].text)

    # alternate lookup for different granule versions
    search_term = './*/Product_Info/PROCESSING_BASELINE'
    processing_baseline = xml_root.findall(search_term)[0].text

    # supported bands for this sensor
    supported_bands = SENSORS[platform_id]['MSI']['bands']

    # TODO: extend this to incorporate a S2b selection
    acqtype = Sentinel2aAcquisition

    # safe archive for S2a has two band name mappings, and we need to map ours
    band_map = {'0': '1', '1': '2', '2': '3', '3': '4', '4': '5', '5': '6',
                '6': '7', '7': '8', '8': '8a', '9': '9', '10': '10',
                '11': '11', '12': '12'}

    # conversion factors
    c1 = {band_map[str(i)]: float(v) for i, v in
          enumerate(['0.001'] * 10 + ['0.0005'] + ['0.0002'] + ['0.00005'])}
    
    # earth -> sun distance correction factor; d2 =  1/ U
    search_term = './*/Product_Image_Characteristics/Reflectance_Conversion/U'
    u = float(xml_root.findall(search_term)[0].text)

    # quantification value
    search_term = './*/Product_Image_Characteristics/QUANTIFICATION_VALUE'
    qv = float(xml_root.findall(search_term)[0].text)

    # exoatmospheric solar irradiance
    solar_irradiance = {}
    for irradiance in xml_root.iter('SOLAR_IRRADIANCE'):
        band_irradiance = irradiance.attrib
        mapped_band_id = band_map[band_irradiance['bandId']]
        solar_irradiance[mapped_band_id] = float(irradiance.text)

    # assume multiple granules
    single_granule_archive = False

    search_term = './*/Product_Info/Product_Organisation/Granule_List/Granules'
    granules = {granule.get('granuleIdentifier'):
                    [imid.text for imid in granule.findall('IMAGE_ID')]
                for granule in xml_root.findall(search_term)}

    if not granules:
        single_granule_archive = True
        granules = {granule.get('granuleIdentifier'):
                        [imid.text for imid in granule.findall('IMAGE_FILE')]
                    for granule in xml_root.findall(search_term[:-1])}

    # resolution groups
    band_groups = {'R10m': ['B02', 'B03', 'B04', 'B08'],
                   'R20m': ['B05', 'B06', 'B07', 'B11', 'B12', 'B8A'],
                   'R60m': ['B01', 'B09', 'B10']}

    granule_groups = {}
    for granule_id, images in granules.items():
        res_groups = {'R10m': [],
                      'R20m': [],
                      'R60m': []}

        granule_xmls = [s for s in archive.namelist() if 'MTD_TL.xml' in s]
        if not granule_xmls:
            pattern = granule_id.replace('MSI', 'MTD')
            pattern = pattern.replace(''.join(['_N', processing_baseline]),
                                      '.xml')

            granule_xmls = [s for s in archive.namelist() if pattern in s]

        granule_xml = archive.read(granule_xmls[0])
        granule_root = ElementTree.XML(granule_xml)

        img_data_path = ''.join(['zip:', pathname, '!', archive.namelist()[0]])

        if not single_granule_archive: 
            img_data_path = ''.join([img_data_path,
                                     pjoin('GRANULE', granule_id, 'IMG_DATA')])

        # acquisition centre datetime
        search_term = './*/SENSING_TIME'
        acq_time = parser.parse(granule_root.findall(search_term)[0].text)

        for image in images:
            # image filename
            img_fname = ''.join([pjoin(img_data_path, image), '.jp2'])

            # band id
            group, band_id = group_helper(img_fname, band_groups)
            if band_id not in supported_bands:
                continue

            # band info stored in sensors.json
            sensor_band_info = supported_bands[band_id]

            # image attributes/metadata
            attrs = {k: v for k, v in sensor_band_info.items()}
            attrs['solar_irradiance'] = solar_irradiance[band_id]
            attrs['d2'] = 1 / u
            attrs['qv'] = qv
            attrs['c1'] = c1[band_id]
            attrs['granule_xml'] = granule_xmls[0]
            band_name = attrs.pop('band_name')

            res_groups[group].append(acqtype(pathname, img_fname, acq_time,
                                             band_name, band_id, attrs))

        granule_groups[granule_id] = {k: sorted(v) for k, v in
                                      res_groups.items()}

    return AcquisitionsContainer(label=basename(pathname),
                                 granules=granule_groups)


class Sentinel2aAcquisition(Acquisition):

    """ Sentinel-2a acquisition. """

    def __init__(self, pathname, uri, acquisition_datetime, band_name='BAND 1',
                 band_id='1', metadata=None):
        super(Sentinel2aAcquisition, self).__init__(pathname, uri,
                                                    acquisition_datetime,
                                                    band_name=band_name,
                                                    band_id=band_id,
                                                    metadata=metadata)

        self.platform_id = 'SENTINEL_2A'
        self.sensor_id = 'MSI'
        self.tle_format = 'S2A%4d%sASNNOR.S00'
        self.tag = 'S2A'
        self.altitude = 786000.0
        self.inclination = 1.721243708316808
        self.omega = 0.001039918
        self.semi_major_axis = 7164137.0
        self.maximum_view_angle = 20.0

        self._gps_file = True
        self._solar_zenith = None

    def read_gps_file(self):
        """
        Returns the recorded gps data as a `pandas.DataFrame`.
        """
        # open the zip archive and get the xml root
        archive = zipfile.ZipFile(self.pathname)
        xml_file = [s for s in archive.namelist() if
                    ("DATASTRIP" in s) & (".xml" in s)][0]
        xml_root = ElementTree.XML(archive.read(xml_file))

        # coordinate transform
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

        # output container
        data = {'longitude': [],
                'latitude': [],
                'altitude': [],
                'timestamp': []}

        # there are a few columns of data that could be of use
        # but for now, just get the location and timestamp from the
        # gps points list
        gps = xml_root.findall('./*/Ephemeris/GPS_Points_List')[0]
        for point in gps.iter('GPS_Point'):
            x, y, z = [float(i) / 1000 for i in
                       point.findtext('POSITION_VALUES').split()]
            gps_time = parser.parse(point.findtext('GPS_TIME'))

            # coordinate transformation
            lon, lat, alt = pyproj.transform(ecef, lla, x, y, z, radians=False)

            data['longitude'].append(lon)
            data['latitude'].append(lat)
            data['altitude'].append(alt)
            data['timestamp'].append(gps_time)
        
        return pandas.DataFrame(data)

    def _retrieve_solar_zenith(self):
        """
        If radiance is not available, calculate it into a
        temporary file.
        Can we use tempfile, that way it is removed from disk
        upon object deletion?
        Or simply keep it in memory?
        Code adapted from https://github.com/umwilm/SEN2COR.
        """
        def rbspline(y_coords, x_coords, zdata):
            """
            A wrapper providing the call signiture for RectBivariateSpline.
            Code adapted from https://github.com/umwilm/SEN2COR.
            """
            y = numpy.arange(zdata.shape[0], dtype=numpy.float32)
            x = numpy.arange(zdata.shape[1], dtype=numpy.float32)
            func = interpolate.RectBivariateSpline(x, y, zdata)

            return func(y_coords, x_coords)

        archive = zipfile.ZipFile(self.pathname)
        xml_root = ElementTree.XML(archive.read(self.granule_xml))
        
        # read the low res solar zenith data
        search_term = './*/Tile_Angles/Sun_Angles_Grid/Zenith/Values_List'
        values = xml_root.findall(search_term)[0]

        dims = (len(values), len(values[0].text.split()))
        data = numpy.zeros(dims, dtype='float32')
        for i, val in enumerate(values.iter('VALUES')):
            data[i] = val.text.split()

        # correct solar_zenith dimensions
        if self.lines < self.samples:
            last_row = int(data[0].size * float(self.lines) /
                           float(self.samples) + 0.5)
            solar_zenith = data[0:last_row, :]
        elif self.samples < self.lines:
            last_col = int(data[1].size * float(self.samples) /
                           float(self.lines) + 0.5)
            solar_zenith = data[:, 0:last_col]
        else:
            solar_zenith = data

        numpy.absolute(solar_zenith, out=solar_zenith)
        numpy.clip(solar_zenith, 0, 70.0, out=solar_zenith)

        # interpolate across the full dimensions of the acquisition
        dims = solar_zenith.shape
        y = numpy.arange(self.lines) / (self.lines - 1) * dims[0]
        x = numpy.arange(self.samples) / (self.samples - 1) * dims[1]

        solar_zenith = numpy.float32(rbspline(y, x, solar_zenith))
        self._solar_zenith = numpy.radians(solar_zenith, out=solar_zenith)

    def radiance_data(self, window=None, out_no_data=-999):
        """
        Return the data as radiance in watts/(m^2*micrometre).
        """
        # retrieve the solar zenith if we haven't already done so
        if self._solar_zenith is None:
            self._retrieve_solar_zenith()

        # Python style index
        if window is None:
            idx = (slice(None, None), slice(None, None))
        else:
            idx = (slice(window[0][0], window[0][1]),
                   slice(window[1][0], window[1][1]))

        # coefficients
        #pylint: disable=unused-argument
        sf = numpy.float32(1 / (self.c1 * self.qv))
        pi_d2 = numpy.float32(numpy.pi * self.d2)
        esun = numpy.float32(self.solar_irradiance / 10)
        solar_zenith = self._solar_zenith[idx]

        # toa reflectance
        data = self.data(window=window)

        # check for no data
        no_data = self.no_data if self.no_data is not None else 0
        nulls = data == no_data

        # inversion
        expr = "((data * esun * cos(solar_zenith) * sf) / pi_d2) * 0.01"
        radiance = numexpr.evaluate(expr)
        radiance[nulls] = out_no_data

        return radiance
