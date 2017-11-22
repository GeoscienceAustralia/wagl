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
from os.path import isdir, join as pjoin, dirname, basename, splitext
import re
import json
import datetime
from xml.etree import ElementTree
import zipfile
from dateutil import parser
from nested_lookup import nested_lookup


from .base import Acquisition, AcquisitionsContainer
from .sentinel import Sentinel2aAcquisition, Sentinel2aAcquisitionOnDisk
from .landsat import ACQUISITION_TYPE, LandsatAcquisition

from ..mtl import load_mtl


with open(pjoin(dirname(__file__), 'sensors.json')) as fo:
    SENSORS = json.load(fo)

def fixname(s):
    """Fix satellite name.
       Performs 'Landsat7' to 'LANDSAT_7', 'LANDSAT8' to 'LANDSAT_8',
       'Landsat-5' to 'LANDSAT_5'.
    """
    return re.sub(r'([a-zA-Z]+)[_-]?(\d)',
                  lambda m: m.group(1).upper() + '_' + m.group(2), s)


def find_in(path, s, suffix='txt'):
    """Search through `path` and its children for the first occurance of a
    file with `s` in its name. Returns the path of the file or `None`. """
    for root, _, files in os.walk(path):
        for f in files:
            if s in f and f.endswith(suffix):
                return os.path.join(root, f)
    return None


def acquisitions(path, hint=None):
    """
    Return an instance of `AcquisitionsContainer` containing the
    acquisitions for each Granule Group (if applicable) and
    each sub Group.
    """

    if hint == 's2_directory':
        container = acquisitions_s2_directory(path)
    elif splitext(path)[1] == '.zip':
        container = acquisitions_via_safe(path)
    else:
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
    supported_bands = SENSORS[platform_id][sensor_id]['band_ids']

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


class Sentinel2AConstants(object):
    """
    Defines shared constants between the two sentinel2 metadata 
    retrieval patterns
    """

    def __init__(self):
        # pylint: disable=invalid-name
        # safe archive for S2a has two band name mappings, and we need to map ours
        self.band_map = {
            '0': '1', '1': '2', '2': '3', '3': '4', '4': '5', '5': '6',
            '6': '7', '7': '8', '8': '8A', '9': '9', '10': '10',
            '11': '11', '12': '12'
        }
        # conversion factors
        self.c1 = {self.band_map[str(i)]: v for i, v in
                   enumerate([0.001] * 10 + [0.0005] + [0.0002] + [0.00005])}

        # radiance (w/m^2)  scale factors
        self.rsf = {self.band_map[str(i)]: v for i, v in
                    enumerate([0.01] * 10 + [None] + [0.002] + [0.0005])}



def acquisitions_s2_directory(pathname):
    """
    Collect the TOA Radiance images for each granule within a scene.
    Multi-granule & multi-resolution hierarchy format.
    Returns a dict of granules, each granule contains a dict of resolutions,
    and each resolution contains a list of acquisition objects.

    it is assumed that the pathname points to a directory containing
    information pulled from AWS S3 granules (s3://sentinel-s2-l1c)
    with additional information retrieved from the productInfo.json
    sitting in a subfolder

    Example:
    {'GRANULE_1': {'R10m': [`acquisition_1`,...,`acquisition_n`],
                   'R20m': [`acquisition_1`,...,`acquisition_n`],
                   'R60m': [`acquisition_1`,...,`acquisition_n`]},
     'GRANULE_N': {'R10m': [`acquisition_1`,...,`acquisition_n`],
                   'R20m': [`acquisition_1`,...,`acquisition_n`],
                   'R60m': [`acquisition_1`,...,`acquisition_n`]}}
    """

    granule_xml = pathname + '/metadata.xml'
    acqtype = Sentinel2aAcquisitionOnDisk
    s2a_const = Sentinel2AConstants()

    search_paths = {
        'datastrip/metadata.xml': [
            {
                'key': 'platform_id',
                'search_path': ".//*/SPACECRAFT_NAME",
                'parse': lambda x: fixname(x[0].text) if x else None,
            },
            {
                'key': 'processing_baseline',
                'search_path': ".//*/PROCESSING_BASELINE",
                'parse': lambda x: x[0].text if x else None,
            },
            {
                'key': 'u',
                'search_path': './/*/Reflectance_Conversion/U',
                'parse': lambda x: float(x[0].text) if x else None,
            },
            {
                'key': 'qv',
                'search_path': './/*/QUANTIFICATION_VALUE',
                'parse': lambda x: float(x[0].text) if x else None,
            },
            {
                'key': 'solar_irradiance_list',
                'search_path': './/*/SOLAR_IRRADIANCE',
                'parse': lambda x: {
                    si.attrib['bandId']: float(si.text) for si in x
                },
            },
        ],
        'metadata.xml': [
            {
                'key': 'granule_id',
                'search_path': './/*/TILE_ID',
                'parse': lambda x: x[0].text if x else None
            },
            {
                'key': 'acq_time',
                'search_path': './/*/SENSING_TIME',
                'parse': lambda x: parser.parse(x[0].text) if x else None,
            },
        ]
    }

    acquisition_data = {}
    for fn, terms in search_paths.items():
        xml_root = None
        with open('/home/ubuntu/granule/' + fn, 'r') as fd:
            xml_root = ElementTree.XML(fd.read())
        for term in terms:
            acquisition_data[term['key']] = term['parse'](xml_root.findall(term['search_path']))

    supported_bands = SENSORS[acquisition_data['platform_id']]['MSI']['band_ids']
    # resolution groups
    band_groups = {'R10m': ['B02', 'B03', 'B04', 'B08'],
                   'R20m': ['B05', 'B06', 'B07', 'B11', 'B12', 'B8A'],
                   'R60m': ['B01', 'B09', 'B10']}

    res_groups = {'R10m': [],
                  'R20m': [],
                  'R60m': []}

    for band_id in s2a_const.band_map:
        if s2a_const.band_map[band_id] not in supported_bands:
            continue

        band_name = 'B{}'.format(s2a_const.band_map[band_id].zfill(2))
        img_fname = pathname + '/' + band_name + '.jp2'

        if not os.path.isfile(img_fname):
            continue

        group = None
        for group, bands in band_groups.items():
            if band_name in bands:
                break
        else:
            continue  # group not found

        attrs = {k: v for k, v in supported_bands[s2a_const.band_map[band_id]].items()}
        attrs['solar_irradiance'] = acquisition_data['solar_irradiance_list'][band_id]
        attrs['d2'] = 1 / acquisition_data['u']
        attrs['qv'] = acquisition_data['qv']
        attrs['c1'] = s2a_const.c1[s2a_const.band_map[band_id]]
        attrs['radiance_scale_factor'] = s2a_const.rsf[s2a_const.band_map[band_id]]
        attrs['granule_xml'] = granule_xml
        band_name = attrs.pop('band_name')

        acq_time = acquisition_data['acq_time']


        res_groups[group].append(acqtype(pathname, img_fname, acq_time,
                                         band_name, band_id, attrs))

    granule_id = acquisition_data['granule_id']

    granule_groups = {}
    granule_groups[granule_id] = {k: sorted(v) for k, v in res_groups.items()}

    return AcquisitionsContainer(label=basename(pathname),
                                 granules=granule_groups)


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

    s2a_const = Sentinel2AConstants()

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
    supported_bands = SENSORS[platform_id]['MSI']['band_ids']

    # TODO: extend this to incorporate a S2b selection
    acqtype = Sentinel2aAcquisition

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
        mapped_band_id = s2a_const.band_map[band_irradiance['bandId']]
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
            attrs['c1'] = s2a_const.c1[band_id]
            attrs['radiance_scale_factor'] = s2a_const.rsf[band_id]
            attrs['granule_xml'] = granule_xmls[0]
            band_name = attrs.pop('band_name')

            res_groups[group].append(acqtype(pathname, img_fname, acq_time,
                                             band_name, band_id, attrs))

        granule_groups[granule_id] = {k: sorted(v) for k, v in
                                      res_groups.items()}

    return AcquisitionsContainer(label=basename(pathname),
                                 granules=granule_groups)
