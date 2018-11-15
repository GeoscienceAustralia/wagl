#!/usr/bin/env python

"""
Various metadata extraction and creation, and writing tools.
"""

from __future__ import absolute_import, print_function
from datetime import datetime as dtime, timezone as dtz
import os
from os.path import dirname
import pwd
import socket
import uuid
import numpy
import pandas
import rasterio
import yaml
from yaml.representer import Representer
import wagl
from wagl.constants import BrdfParameters, DatasetName, POINT_FMT, GroupName
from wagl.constants import BandType
from wagl.hdf5 import write_scalar, read_h5_table, read_scalar

yaml.add_representer(numpy.int8, Representer.represent_int)
yaml.add_representer(numpy.uint8, Representer.represent_int)
yaml.add_representer(numpy.int16, Representer.represent_int)
yaml.add_representer(numpy.uint16, Representer.represent_int)
yaml.add_representer(numpy.int32, Representer.represent_int)
yaml.add_representer(numpy.uint32, Representer.represent_int)
yaml.add_representer(numpy.int, Representer.represent_int)
yaml.add_representer(numpy.int64, Representer.represent_int)
yaml.add_representer(numpy.uint64, Representer.represent_int)
yaml.add_representer(numpy.float, Representer.represent_float)
yaml.add_representer(numpy.float32, Representer.represent_float)
yaml.add_representer(numpy.float64, Representer.represent_float)
yaml.add_representer(numpy.ndarray, Representer.represent_list)


def extract_ancillary_metadata(fname):
    """
    Extracts the change (last metadata change), modified,
    accessed, and owner user id.

    :param fname:
        A string containing the full file pathname to a file
        on disk.

    :return:
        A `dictionary` with keys `ctime`, `mtime`, `atime`,
        and `owner`.
    """

    def _get_utc_datetime(timestamp):
        return dtime.utcfromtimestamp(timestamp).replace(tzinfo=dtz.utc)

    res = {}
    fstat = os.stat(fname)
    res['ctime'] = _get_utc_datetime(fstat.st_ctime)
    res['mtime'] = _get_utc_datetime(fstat.st_mtime)
    res['atime'] = _get_utc_datetime(fstat.st_atime)
    res['owner_id'] = str(fstat.st_uid)
    return res


def get_system_information():
    utc_now = dtime.utcnow().replace(tzinfo=dtz.utc).isoformat()
    system_info = {'uname': ' '.join(os.uname()),
                   'hostname': socket.getfqdn(),
                   'runtime_id': str(uuid.uuid1()),
                   'time_processed': utc_now}
    return system_info


def read_metadata_tags(fname, bands):
    """
    Retrieves the metadata tags for a list of bands from a `GDAL`
    compliant dataset.

    :param fname:
        A string containing the full file pathname to a file
        on disk.

    :param bands:
        A `list` containing the band numbers (1 -> n) from which
	to retreive the metadata tags.

    :return:
        A `pandas.DataFrame`.
    """
    with rasterio.open(fname) as ds:
        tag_data = {k: [] for k in ds.tags(1).keys()}
        for band in bands:
            tags = ds.tags(band)
            for tag in tags:
                tag_data[tag].append(tags[tag])

    return pandas.DataFrame(tag_data)


def create_ard_yaml(acquisitions, ancillary_group, out_group, normalization_angle, sbt=False):
    """
    Write the NBAR metadata captured during the entire workflow to a
    HDF5 SCALAR dataset using the yaml document format.

    :param acquisitions:
        A `list` of `Acquisition` instances.

    :param ancillary_group:
        The root HDF5 `Group` that contains the ancillary data
        collected via wagl.ancillary.collect_ancillary>

    :param out_group:
        A `h5py.Group` object opened for write access.

    :param sbt:
        A `bool` indicating whether to create the
        Surface Brightness Temperature yaml dataset.
        Default is False.

    :return:
        None; The yaml document is written to the HDF5 file.
    """
    def load_sbt_ancillary(group):
        """
        Load the sbt ancillary data retrieved during the worlflow.
        """
        point_data = {DatasetName.DEWPOINT_TEMPERATURE.value: {},
                      DatasetName.SURFACE_GEOPOTENTIAL.value: {},
                      DatasetName.TEMPERATURE_2M.value: {},
                      DatasetName.SURFACE_RELATIVE_HUMIDITY.value: {},
                      DatasetName.GEOPOTENTIAL.value: {},
                      DatasetName.RELATIVE_HUMIDITY.value: {},
                      DatasetName.TEMPERATURE.value: {}}

        npoints = group[DatasetName.COORDINATOR.value].shape[0]
        for point in range(npoints):
            pnt_grp = group[POINT_FMT.format(p=point)]
            lonlat = tuple(pnt_grp.attrs['lonlat'])

            # scalars
            dname = DatasetName.DEWPOINT_TEMPERATURE.value
            point_data[dname][lonlat] = read_scalar(pnt_grp, dname)

            dname = DatasetName.SURFACE_GEOPOTENTIAL.value
            point_data[dname][lonlat] = read_scalar(pnt_grp, dname)

            dname = DatasetName.TEMPERATURE_2M.value
            point_data[dname][lonlat] = read_scalar(pnt_grp, dname)

            dname = DatasetName.SURFACE_RELATIVE_HUMIDITY.value
            point_data[dname][lonlat] = read_scalar(pnt_grp, dname)

            # tables
            dname = DatasetName.GEOPOTENTIAL.value
            dset = pnt_grp[dname]
            attrs = {k: v for k, v in dset.attrs.items()}
            df = read_h5_table(pnt_grp, dname)
            for column in df.columns:
                attrs[column] = df[column].values
            point_data[dname][lonlat] = attrs

            dname = DatasetName.RELATIVE_HUMIDITY.value
            dset = pnt_grp[dname]
            attrs = {k: v for k, v in dset.attrs.items()}
            df = read_h5_table(pnt_grp, dname)
            for column in df.columns:
                attrs[column] = df[column].values
            point_data[dname][lonlat] = attrs
            
            dname = DatasetName.TEMPERATURE.value
            dset = pnt_grp[dname]
            attrs = {k: v for k, v in dset.attrs.items()}
            df = read_h5_table(pnt_grp, dname)
            for column in df.columns:
                attrs[column] = df[column].values
            point_data[dname][lonlat] = attrs

        return point_data

    def load_ancillary(acquisitions, fid, sbt=False):
        """
        Load the ancillary data retrieved during the workflow.
        """
        # retrieve the averaged ancillary if available
        anc_grp = fid.get(GroupName.ANCILLARY_AVG_GROUP.value)
        if anc_grp is None:
            anc_grp = fid

        dname = DatasetName.AEROSOL.value
        aerosol_data = read_scalar(anc_grp, dname)
        dname = DatasetName.WATER_VAPOUR.value
        water_vapour_data = read_scalar(anc_grp, dname)
        dname = DatasetName.OZONE.value
        ozone_data = read_scalar(anc_grp, dname)
        dname = DatasetName.ELEVATION.value
        elevation_data = read_scalar(anc_grp, dname)

        ancillary = {'aerosol': aerosol_data,
                     'water_vapour': water_vapour_data,
                     'ozone': ozone_data,
                     'elevation': elevation_data,
                     'normalization_angle': normalization_angle}

        if sbt:
            sbt_ancillary = load_sbt_ancillary(fid)
            for key in sbt_ancillary:
                ancillary[key] = sbt_ancillary[key]
        else:
            for acq in acquisitions:
                if acq.band_type == BandType.THERMAL:
                    continue

                bn = acq.band_name
                for param in BrdfParameters:
                    fmt = DatasetName.BRDF_FMT.value
                    dname = fmt.format(band_name=bn, parameter=param.value)
                    dset = fid[dname]
                    key = dname.lower().replace('-', '_')
                    ancillary[key] = {k: v for k, v in dset.attrs.items()}
                    ancillary[key]['value'] = dset[()]
                    ancillary[key]['type'] = key

        return ancillary

    acquisition = acquisitions[0]
    level1_path = acquisition.pathname
    acq_datetime = (
        acquisition.acquisition_datetime
        .replace(tzinfo=dtz.utc)
        .isoformat()
    )
    source_info = {'source_level1': level1_path,
                   'acquisition_datetime': acq_datetime,
                   'platform_id': acquisition.platform_id,
                   'sensor_id': acquisition.sensor_id}

    # ancillary metadata tracking
    for key, value in extract_ancillary_metadata(level1_path).items():
        if isinstance(value, dtime):
            source_info[key] = value.isoformat()
        else:
            source_info[key] = value

    # load the ancillary and remove fields not of use to ODC
    ancillary = load_ancillary(acquisitions, ancillary_group, sbt)
    for item in ancillary:
        for remove in ['CLASS', 'VERSION', 'query_date', 'data_source']:
            ancillary[item].pop(remove, None)

    software_versions = {'wagl': {'version': wagl.__version__,
                                  'repo_url': 'https://github.com/GeoscienceAustralia/wagl.git'}, # pylint: disable=line-too-long
                         'modtran': {'version': '5.2.1',
                                     'repo_url': 'http://www.ontar.com/software/productdetails.aspx?item=modtran'} # pylint: disable=line-too-long
                        }

    algorithm = {}
    if sbt:
        dname = DatasetName.SBT_YAML.value
        algorithm['sbt_doi'] = 'TODO'
    else:
        dname = DatasetName.NBAR_YAML.value
        algorithm['algorithm_version'] = 2.0
        algorithm['arg25_doi'] = 'http://dx.doi.org/10.4225/25/5487CC0D4F40B'
        algorithm['nbar_doi'] = 'http://dx.doi.org/10.1109/JSTARS.2010.2042281'
        algorithm['nbar_terrain_corrected_doi'] = 'http://dx.doi.org/10.1016/j.rse.2012.06.018' # pylint: disable=line-too-long

    metadata = {'system_information': get_system_information(),
                'source_datasets': source_info,
                'ancillary': ancillary,
                'algorithm_information': algorithm,
                'software_versions': software_versions}
    
    # output
    yml_data = yaml.dump(metadata, default_flow_style=False)
    write_scalar(yml_data, dname, out_group, attrs={'file_format': 'yaml'})


def create_pq_yaml(acquisition, ancillary, tests_run, out_group):
    """
    Write the PQ metadata captured during the entire workflow to a
    HDF5 SCALAR dataset using the yaml document format.

    :param acquisition:
        An instance of `acquisition`.

    :param ancillary:
        A dict containing the ancillary information.

    :param test_run:
        A dict containing the key/value pairs of tests and whether
        or not a given test was run.

    :param out_group:
        A `h5py.Group` object opened for write access.

    :return:
        None; The yaml document is written to the HDF5 file.
    """

    source_info = {'source_l1t': dirname(acquisition.dir_name),
                   'source_reflectance': 'NBAR'}

    algorithm = {'software_version': wagl.__version__,
                 'software_repository': 'https://github.com/GeoscienceAustralia/wagl.git', # pylint: disable=line-too-long
                 'pq_doi': 'http://dx.doi.org/10.1109/IGARSS.2013.6723746'}
    
    metadata = {'system_information': get_system_information(),
                'source_data': source_info,
                'algorithm_information': algorithm,
                'ancillary': ancillary,
                'tests_run': tests_run}

    # output
    dname = DatasetName.PQ_YAML.value
    yml_data = yaml.dump(metadata, default_flow_style=False)
    write_scalar(yml_data, dname, out_group, attrs={'file_format': 'yaml'})
