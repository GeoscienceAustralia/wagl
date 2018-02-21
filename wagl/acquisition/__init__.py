"""
The acquisition is a core object passed throughout wagl,
which provides a lightweight `Acquisition` class.
An `Acquisition` instance provides data read methods,
and contains various attributes/metadata regarding the band
as well as the sensor, and sensor's platform.

An `AcquisitionsContainer` provides a container like structure
for handling scenes consisting of multiple Granules/Tiles,
and of differing resolutions.
"""
from __future__ import absolute_import, print_function
from collections import OrderedDict
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
from .sentinel import Sentinel2aAcquisition, Sentinel2bAcquisition, s2_index_to_band_id
from .sentinel import Sentinel2aSinergiseAcquisition, Sentinel2bSinergiseAcquisition
from .landsat import ACQUISITION_TYPE, LandsatAcquisition

from ..mtl import load_mtl


# resolution group format
RESG_FMT = "RES-GROUP-{}"


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

    if hint == 's2_sinergise':
        container = acquisitions_s2_sinergise(path)
    elif splitext(path)[1] == '.zip':
        container = acquisitions_via_safe(path)
    else:
        try:
            container = acquisitions_via_mtl(path)
        except OSError:
            raise IOError("No acquisitions found in: {}".format(path))

    return container


def create_resolution_groups(acqs):
    """
    Given a list of acquisitions, return an OrderedDict containing
    groups of acquisitions based on the resolution.
    The order of groups is highest to lowest resolution.
    Resolution group names are given as:

    * RES-GROUP-0, RES-GROUP-1, ..., RES-GROUP-N

    :param acqs:
        A `list` of acquisition objects.

    :return:
        An OrderedDict detailed as follows:

        {'RES-GROUP-0: [acquisition, acquisition],
         'RES-GROUP-1: [acquisition, acquisition],
         'RES-GROUP-N: [acquisition, acquisition]}
    """
    fmt = "RES-GROUP-{}"
    # 0 -> n resolution sets (higest res to lowest res)
    resolutions = sorted(set([acq.resolution for acq in acqs]))
    res_groups = OrderedDict([(fmt.format(i), []) for i, _ in
                              enumerate(resolutions)])

    for acq in sorted(acqs):
        group = fmt.format(resolutions.index(acq.resolution))
        res_groups[group].append(acq)

    return res_groups


def acquisitions_via_mtl(pathname):
    """
    Obtain a list of Acquisition objects from `pathname`.
    The argument `pathname` can be a MTL file or a directory name.
    If `pathname` is a directory then the MTL file will be search
    for in the directory and its children.
    Returns an instance of `AcquisitionsContainer`.
    """

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

    # granule id
    granule_id = nested_lookup('landsat_scene_id', data)[0]

    # bands to ignore
    ignore = ['band_quality']

    # supported bands for the given platform & sensor id's
    band_configurations = SENSORS[platform_id][sensor_id]['band_ids']

    acqs = []
    for band in bands_:
        if band in ignore:
            continue

        # band id
        if 'vcid' in band:
            band_id = band.replace('_vcid_', '').replace('band', '').strip('_')
        else:
            band_id = band.replace('band', '').strip('_')

        # band info stored in sensors.json
        sensor_band_info = band_configurations.get(band_id, {})

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
        if attrs.get('supported_band'):
            attrs['solar_azimuth'] = solar_azimuth
            attrs['solar_elevation'] = solar_elevation
            attrs['min_radiance'] = min_rad
            attrs['max_radiance'] = max_rad
            attrs['min_quantize'] = min_quant
            attrs['max_quantize'] = max_quant

        # band_name is an internal property of acquisitions class
        band_name = attrs.pop('band_name', band_id)

        acqs.append(acqtype(pathname, fname, acq_datetime, band_name, band_id,
                            attrs))

    # resolution groups dict
    res_groups = create_resolution_groups(acqs)

    return AcquisitionsContainer(label=basename(pathname),
                                 granules={granule_id: res_groups})


def acquisitions_s2_sinergise(pathname):
    """
    Collect the TOA Radiance images for each granule within a scene.
    Multi-granule & multi-resolution hierarchy format.
    Returns an instance of `AcquisitionsContainer`.

    it is assumed that the pathname points to a directory containing
    information pulled from AWS S3 granules (s3://sentinel-s2-l1c)
    with additional information retrieved from the productInfo.json
    sitting in a subfolder
    """

    granule_xml = pathname + '/metadata.xml'

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
                    s2_index_to_band_id(si.attrib['bandId']): float(si.text) for si in x
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
        with open(pathname + '/' + fn, 'rb') as fd:
            xml_root = ElementTree.XML(fd.read())
        for term in terms:
            acquisition_data[term['key']] = term['parse'](xml_root.findall(term['search_path']))

    band_configurations = SENSORS[acquisition_data['platform_id']]['MSI']['band_ids']
    esa_ids = ['B02', 'B03', 'B04', 'B08', 'B05', 'B06', 'B07', 'B11',
               'B12', 'B8A', 'B01', 'B09', 'B10']

    if 'S2A' in acquisition_data['granule_id']:
        acqtype = Sentinel2aSinergiseAcquisition
    else:
        # assume it is S2B
        acqtype = Sentinel2bSinergiseAcquisition

    acqs = []
    for band_id in band_configurations:

        # If it is a configured B-format transform it to the correct format
        if re.match('[0-9].?', band_id):
            band_name = 'B{}'.format(band_id.zfill(2))

        img_fname = pathname + '/' + band_name + '.jp2'

        if not os.path.isfile(img_fname):
            continue

        if band_name not in esa_ids:
            continue

        attrs = {k: v for k, v in band_configurations[band_id].items()}
        if attrs.get('supported_band'):
            attrs['solar_irradiance'] = acquisition_data['solar_irradiance_list'][band_id]
            attrs['d2'] = 1 / acquisition_data['u']
            attrs['qv'] = acquisition_data['qv']

        # Required attribute for packaging
        attrs['granule_xml'] = granule_xml

        # band_name is an internal property of acquisitions class
        band_name = attrs.pop('band_name', band_id)

        acq_time = acquisition_data['acq_time']


        acqs.append(acqtype(pathname, img_fname, acq_time, band_name, band_id,
                            attrs))

    granule_id = acquisition_data['granule_id']

    # resolution groups dict
    granule_groups = {}
    granule_groups[granule_id] = create_resolution_groups(acqs)

    return AcquisitionsContainer(label=basename(pathname),
                                 granules=granule_groups)


def acquisitions_via_safe(pathname):
    """
    Collect the TOA Radiance images for each granule within a scene.
    Multi-granule & multi-resolution hierarchy format.
    Returns an instance of `AcquisitionsContainer`.
    """
    def band_id_helper(fname, esa_ids):
        """
        A helper function to find the band_id.
        """
        # TODO: do we need this func any more as res groups are now
        # derived rather pre-determined
        for esa_id in esa_ids:
            if esa_id in fname:
                band_id = re.sub(r'B[0]?', '', esa_id)
                return band_id
        return None

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
    band_configurations = SENSORS[platform_id]['MSI']['band_ids']

    if basename(pathname)[0:3] == 'S2A':
        acqtype = Sentinel2aAcquisition
    else:
        # assume it is S2B
        acqtype = Sentinel2bAcquisition

    # earth -> sun distance correction factor; d2 =  1/ U
    search_term = './*/Product_Image_Characteristics/Reflectance_Conversion/U'
    u = float(xml_root.findall(search_term)[0].text)

    # quantification value
    search_term = './*/Product_Image_Characteristics/QUANTIFICATION_VALUE'
    qv = float(xml_root.findall(search_term)[0].text)

    # exoatmospheric solar irradiance
    solar_irradiance = {}
    for irradiance in xml_root.iter('SOLAR_IRRADIANCE'):
        band_id = s2_index_to_band_id(irradiance.attrib['bandId'])
        solar_irradiance[band_id] = float(irradiance.text)

    search_term = './*/Product_Info/Product_Organisation/Granule_List/Granules'
    grn_elements = xml_root.findall(search_term)

    # handling multi vs single granules + variants of each type
    if not grn_elements:
        grn_elements = xml_root.findall(search_term[:-1])

    if grn_elements[0].findtext('IMAGE_ID'):
        search_term = 'IMAGE_ID'
    else:
        search_term = 'IMAGE_FILE'

    granules = {granule.get('granuleIdentifier'):
                    [imid.text for imid in granule.findall(search_term)]
                for granule in grn_elements}

    # ESA image ids
    esa_ids = ['B02', 'B03', 'B04', 'B08', 'TCI', 'B05', 'B06', 'B07', 'B11',
               'B12', 'B8A', 'B01', 'B09', 'B10']

    granule_groups = {}
    for granule_id, images in granules.items():

        granule_xmls = [s for s in archive.namelist() if 'MTD_TL.xml' in s]
        if not granule_xmls:
            pattern = granule_id.replace('MSI', 'MTD')
            pattern = pattern.replace(''.join(['_N', processing_baseline]),
                                      '.xml')

            granule_xmls = [s for s in archive.namelist() if pattern in s]

        granule_xml = archive.read(granule_xmls[0])
        granule_root = ElementTree.XML(granule_xml)

        # handling different metadata versions for image paths
        img_data_path = ''.join(['zip:', pathname, '!', archive.namelist()[0]])
        if basename(images[0]) == images[0]:
            img_data_path = ''.join([img_data_path,
                                     pjoin('GRANULE', granule_id, 'IMG_DATA')])

        # acquisition centre datetime
        search_term = './*/SENSING_TIME'
        acq_time = parser.parse(granule_root.findall(search_term)[0].text)

        acqs = []
        for image in images:
            # image filename
            img_fname = ''.join([pjoin(img_data_path, image), '.jp2'])

            # band id
            band_id = band_id_helper(img_fname, esa_ids)

            # band info stored in sensors.json
            sensor_band_info = band_configurations.get(band_id)
            attrs = {k: v for k, v in sensor_band_info.items()}

            # image attributes/metadata
            if sensor_band_info.get('supported_band'):
                attrs['solar_irradiance'] = solar_irradiance[band_id]
                attrs['d2'] = 1 / u
                attrs['qv'] = qv

            # Required attribute for packaging
            attrs['granule_xml'] = granule_xmls[0]

            # band_name is an internal property of acquisitions class
            band_name = attrs.pop('band_name', band_id)

            acqs.append(acqtype(pathname, img_fname, acq_time, band_name,
                                band_id, attrs))

        # resolution groups dict
        granule_groups[granule_id] = create_resolution_groups(acqs)

    return AcquisitionsContainer(label=basename(pathname),
                                 granules=granule_groups)
