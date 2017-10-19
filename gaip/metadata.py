#!/usr/bin/env python

"""
Various metadata extraction and creation, and writing tools.
"""

from __future__ import absolute_import, print_function
from datetime import datetime as dtime
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
import gaip
from gaip.constants import BrdfParameters, DatasetName, POINT_FMT, GroupName
from gaip.constants import BandType
from gaip.hdf5 import write_scalar, read_h5_table, read_scalar

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
        A `dictionary` with keys `change`, `modified`, `accessed`,
        and `user`.
    """
    res = {}
    fstat = os.stat(fname)
    res['change'] = dtime.utcfromtimestamp(fstat.st_ctime)
    res['modified'] = dtime.utcfromtimestamp(fstat.st_mtime)
    res['accessed'] = dtime.utcfromtimestamp(fstat.st_atime)
    res['user'] = pwd.getpwuid(fstat.st_uid).pw_gecos
    return res


def read_meatadata_tags(fname, bands):
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


def create_ard_yaml(acquisitions, ancillary_group, out_group, sbt=False):
    """
    Write the NBAR metadata captured during the entire workflow to a
    HDF5 SCALAR dataset using the yaml document format.

    :param acquisitions:
        A `list` of `Acquisition` instances.

    :param ancillary_group:
        The root HDF5 `Group` that contains the ancillary data
        collected via gaip.ancillary.collect_ancillary>

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
        point_data = {DatasetName.dewpoint_temperature.value: {},
                      DatasetName.surface_geopotential.value: {},
                      DatasetName.temperature_2m.value: {},
                      DatasetName.surface_relative_humidity.value: {},
                      DatasetName.geopotential.value: {},
                      DatasetName.relative_humidity.value: {},
                      DatasetName.temperature.value: {}}

        npoints = group[DatasetName.coordinator.value].shape[0]
        for point in range(npoints):
            pnt_grp = group[POINT_FMT.format(p=point)]
            lonlat = tuple(pnt_grp.attrs['lonlat'])

            # scalars
            dname = DatasetName.dewpoint_temperature.value
            point_data[dname][lonlat] = read_scalar(pnt_grp, dname)

            dname = DatasetName.surface_geopotential.value
            point_data[dname][lonlat] = read_scalar(pnt_grp, dname)

            dname = DatasetName.temperature_2m.value
            point_data[dname][lonlat] = read_scalar(pnt_grp, dname)

            dname = DatasetName.surface_relative_humidity.value
            point_data[dname][lonlat] = read_scalar(pnt_grp, dname)

            # tables
            dname = DatasetName.geopotential.value
            dset = pnt_grp[dname]
            attrs = {k: v for k, v in dset.attrs.items()}
            df = read_h5_table(pnt_grp, dname)
            for column in df.columns:
                attrs[column] = df[column].values
            point_data[dname][lonlat] = attrs

            dname = DatasetName.relative_humidity.value
            dset = pnt_grp[dname]
            attrs = {k: v for k, v in dset.attrs.items()}
            df = read_h5_table(pnt_grp, dname)
            for column in df.columns:
                attrs[column] = df[column].values
            point_data[dname][lonlat] = attrs
            
            dname = DatasetName.temperature.value
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
        anc_grp = fid.get(GroupName.ancillary_avg_group.value)
        if anc_grp is None:
            anc_grp = fid

        dname = DatasetName.aerosol.value
        aerosol_data = read_scalar(anc_grp, dname)
        dname = DatasetName.water_vapour.value
        water_vapour_data = read_scalar(anc_grp, dname)
        dname = DatasetName.ozone.value
        ozone_data = read_scalar(anc_grp, dname)
        dname = DatasetName.elevation.value
        elevation_data = read_scalar(anc_grp, dname)

        ancillary = {'aerosol': aerosol_data,
                     'water_vapour': water_vapour_data,
                     'ozone': ozone_data,
                     'elevation': elevation_data}

        if sbt:
            sbt_ancillary = load_sbt_ancillary(fid)
            for key in sbt_ancillary:
                ancillary[key] = sbt_ancillary[key]
        else:
            brdf_data = {}
            for acq in acquisitions:
                if acq.band_type == BandType.Thermal:
                    continue

                brdf = {}
                bn = acq.band_name
                for param in BrdfParameters:
                    fmt = DatasetName.brdf_fmt.value
                    dset = fid[fmt.format(band_name=bn, parameter=param.name)]
                    brdf[param.name] = {k: v for k, v in dset.attrs.items()}
                    brdf[param.name]['value'] = dset[()]

                brdf_data[bn] = brdf

            ancillary['brdf'] = brdf_data

        return ancillary

    acquisition = acquisitions[0]
    level1_path = acquisition.pathname
    source_info = {'source_level1': level1_path,
                   'acquisition_datetime': acquisition.acquisition_datetime,
                   'platform_id': acquisition.platform_id,
                   'sensor_id': acquisition.sensor_id}

    # ancillary metadata tracking
    for key, value in extract_ancillary_metadata(level1_path).items():
        source_info[key] = value

    ancillary = load_ancillary(acquisitions, ancillary_group, sbt)

    software_versions = {'gaip': {'version': gaip.__version__,
                                  'repo_url': 'https://github.com/GeoscienceAustralia/gaip.git'}, # pylint: disable=line-too-long
                         'modtran': {'version': '5.2.1',
                                     'repo_url': 'http://www.ontar.com/software/productdetails.aspx?item=modtran'} # pylint: disable=line-too-long
                        }

    algorithm = {}
    if sbt:
        dname = DatasetName.sbt_yaml.value
        algorithm['sbt_doi'] = 'TODO'
    else:
        dname = DatasetName.nbar_yaml.value
        algorithm['algorithm_version'] = 2.0
        algorithm['arg25_doi'] = 'http://dx.doi.org/10.4225/25/5487CC0D4F40B'
        algorithm['nbar_doi'] = 'http://dx.doi.org/10.1109/JSTARS.2010.2042281'
        algorithm['nbar_terrain_corrected_doi'] = 'http://dx.doi.org/10.1016/j.rse.2012.06.018' # pylint: disable=line-too-long

    system_info = {'uname': ' '.join(os.uname()),
                   'hostname': socket.getfqdn(),
                   'runtime_id': uuid.uuid1(),
                   'time_processed': dtime.utcnow().isoformat()}

    metadata = {'system_information': system_info,
                'source_data': source_info,
                'ancillary_data': ancillary,
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
    system_info = {'uname': ' '.join(os.uname()),
                   'hostname': socket.getfqdn(),
                   'runtime_id': uuid.uuid1(),
                   'time_processed': dtime.utcnow().isoformat()}

    source_info = {'source_l1t': dirname(acquisition.dir_name),
                   'source_reflectance': 'NBAR'}

    algorithm = {'software_version': gaip.__version__,
                 'software_repository': 'https://github.com/GeoscienceAustralia/gaip.git', # pylint: disable=line-too-long
                 'pq_doi': 'http://dx.doi.org/10.1109/IGARSS.2013.6723746'}
    
    metadata = {'system_information': system_info,
                'source_data': source_info,
                'algorithm_information': algorithm,
                'ancillary': ancillary,
                'tests_run': tests_run}

    # output
    dname = DatasetName.pq_yaml.value
    yml_data = yaml.dump(metadata, default_flow_style=False)
    write_scalar(yml_data, dname, out_group, attrs={'file_format': 'yaml'})
