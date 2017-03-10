#!/usr/bin/env python

from __future__ import absolute_import
from datetime import datetime as dtime
import os
import pwd
import subprocess
import numpy
import yaml
from yaml.representer import Representer
from gaip import constants


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


def write_nbar_yaml(acquisition, level1_path, ozone_data, aerosol_data,
                    water_vapour_data, elevation_data, brdf_data, out_fname):
    """
    Write the NBAR metadata captured during the entire workflow to disk
    using the yaml document format.
    """

    source_info = {}
    source_info['source_scene'] = level1_path
    source_info['scene_centre_datetime'] = acquisition.scene_centre_datetime
    source_info['platform'] = acquisition.spacecraft_id
    source_info['sensor'] = acquisition.sensor_id
    source_info['path'] = acquisition.path
    source_info['row'] = acquisition.row

    # ancillary metadata tracking
    md = extract_ancillary_metadata(level1_path)
    for key in md:
        source_info[key] = md[key]

    ancillary = {}
    ancillary['aerosol'] = aerosol_data
    ancillary['water_vapour'] = water_vapour_data
    ancillary['ozone'] = ozone_data
    ancillary['elevation'] = elevation_data

    # Get the required BRDF LUT & factors list
    nbar_constants = constants.NBARConstants(acquisition.spacecraft_id,
                                                  acquisition.sensor_id)

    bands = nbar_constants.get_brdf_lut()
    brdf_factors = nbar_constants.get_brdf_factors()

    brdf = {}
    band_fmt = 'band_{}'
    for band in bands:
        data = {}
        for factor in brdf_factors:
            data[factor] = brdf_data[(band, factor)]
        brdf[band_fmt.format(band)] = data

    ancillary['brdf'] = brdf

    # TODO (a) retrieve software version from git once deployed
    algorithm = {}
    algorithm['algorithm_version'] = 2.0 # hardcode for now see TODO (a)
    algorithm['software_repository'] = ('https://github.com/'
                                        'GeoscienceAustralia/'
                                        'ga-neo-landsat-processor.git')
    algorithm['arg25_doi'] = 'http://dx.doi.org/10.4225/25/5487CC0D4F40B'
    algorithm['nbar_doi'] = 'http://dx.doi.org/10.1109/JSTARS.2010.2042281'
    algorithm['nbar_terrain_corrected_doi'] = ('http://dx.doi.org/10.1016/'
                                               'j.rse.2012.06.018')

    system_info = {}
    proc = subprocess.Popen(['uname', '-a'], stdout=subprocess.PIPE)
    system_info['node'] = proc.stdout.read()
    system_info['time_processed'] = dtime.utcnow()

    metadata = {}
    metadata['system_information'] = system_info
    metadata['source_data'] = source_info
    metadata['ancillary_data'] = ancillary
    metadata['algorithm_information'] = algorithm
    
    # Account for NumPy dtypes
    yaml.add_representer(numpy.float, Representer.represent_float)
    yaml.add_representer(numpy.float32, Representer.represent_float)
    yaml.add_representer(numpy.float64, Representer.represent_float)

    # output
    with open(out_fname, 'w') as src:
        yaml.dump(metadata, src, default_flow_style=False)
