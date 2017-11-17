#!/usr/bin/env python

"""
A recursive utility that extracts all SCALAR, IMAGE & TABLE datasets
from a HDF5 file. Any hierarchial structure will be replicated on
disk as directories.
"""

import os
from os.path import join as pjoin, normpath, dirname, exists, basename
from posixpath import join as ppjoin
from posixpath import basename as pbasename
from functools import partial
import argparse
import numpy
import h5py
import yaml
from yaml.representer import Representer

from gaip.data import write_img
from gaip.geobox import GriddedGeoBox
from gaip.hdf5 import read_h5_table

IGNORE = ['crs_wkt', 'geotransform']

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


def convert_image(dataset, output_directory):
    """
    Converts a HDF5 `IMAGE` Class dataset to a compressed GeoTiff,
    with deflate zlevel 1 compression.
    Any attributes stored with the image will be written as dataset
    level metadata tags, and not band level tags.
    All attributes will also be written to a yaml file.

    :param dataset:
        A HDF5 `IMAGE` Class dataset.

    :param output_directory:
        A filesystem path to the directory that will be the root
        directory for any images extracted.

    :return:
        None, outputs are written directly to disk.
    """
    geobox = GriddedGeoBox.from_dataset(dataset)
    tags = {k: v for k, v in dataset.attrs.items() if k not in IGNORE}
    if 'no_data_value' in tags:
        no_data = tags.pop('no_data_value')
    else:
        no_data = None

    tags['history'] = "Converted from HDF5 IMAGE to GeoTiff."

    # TODO: get x & y chunks from 3D images
    kwargs = {'fmt': 'GTiff',
              'geobox': geobox,
              'compress': 'deflate',
              'options': {'zlevel': 1},
              'tags': tags,
              'nodata': no_data}

    base_fname = pjoin(output_directory, normpath(dataset.name.strip('/')))
    out_fname = ''.join([base_fname, '.tif'])

    if not exists(dirname(out_fname)):
        os.makedirs(dirname(out_fname))

    write_img(dataset, out_fname, **kwargs)

    out_fname = ''.join([base_fname, '.yaml'])
    tags = {k: v for k, v in dataset.attrs.items()}
    with open(out_fname, 'w') as src:
        yaml.dump(tags, src, default_flow_style=False, indent=4)


def convert_table(group, dataset_name, output_directory):
    """
    Converts a HDF5 `TABLE` Class dataset to a csv file.
    All attributes will be written to a yaml file.

    :param dataset_name:
        A `str` containing the pathname to the HDF5 `Table` Class
        dataset.

    :param output_directory:
        A filesystem path to the directory that will be the root
        directory for any tables extracted.

    :return:
        None, outputs are written directly to disk.
    """
    df = read_h5_table(group, dataset_name)
    dataset = group[dataset_name]
    tags = {k: v for k, v in dataset.attrs.items()}

    base_fname = pjoin(output_directory, normpath(dataset_name.strip('/')))
    out_fname = ''.join([base_fname, '.csv'])

    if not exists(dirname(out_fname)):
        os.makedirs(dirname(out_fname))

    df.to_csv(out_fname)

    out_fname = ''.join([base_fname, '.yaml'])
    with open(out_fname, 'w') as src:
        yaml.dump(tags, src, default_flow_style=False, indent=4)


def convert_scalar(dataset, output_directory):
    """
    Converts a HDF5 scalar dataset to a yaml file.
    All attributes will be included in the yaml file.

    :param dataset:
        A HDF5 scalar dataset.

    :param output_directory:
        A filesystem path to the directory that will be the root
        directory for any images extracted.

    :return:
        None, outputs are written directly to disk.
    """
    if dataset.attrs.get('file_format') == 'yaml':
        tags = yaml.load(dataset[()])
    else:
        tags = {k: v for k, v in dataset.attrs.items()}
        data = dataset[()]
        if isinstance(data, numpy.bytes_):
            data = data.decode('utf-8')
        tags[basename(dataset.name)] = data

    base_fname = pjoin(output_directory, normpath(dataset.name.strip('/')))
    out_fname = ''.join([base_fname, '.yaml'])

    if not exists(dirname(out_fname)):
        os.makedirs(dirname(out_fname))

    with open(out_fname, 'w') as src:
        yaml.dump(tags, src, default_flow_style=False, indent=4)


def extract(output_directory, group, name):
    """
    A simple utility that sends an object to the appropriate
    extraction utility.
    """
    dataset_name = ppjoin(group.name, name.decode('utf-8'))
    obj = group[dataset_name]
    obj_class = obj.attrs.get('CLASS')

    if obj_class == 'IMAGE':
        convert_image(obj, output_directory)
    elif obj_class == 'TABLE':
        convert_table(group, dataset_name, output_directory)
    elif obj_class == 'SCALAR':
        convert_scalar(obj, output_directory)
    else:
        return None


def run(fname, outdir, pathname):
    """ Run dataset conversion tree. """
    # note: lower level h5py access is required in order to visit links
    with h5py.File(fname, 'r') as fid:
        if pathname in fid:
            obj = fid[pathname]
            if isinstance(obj, h5py.Group):
                root = h5py.h5g.open(obj.id, b'.')
                root.links.visit(partial(extract, outdir, obj))
            elif isinstance(obj, h5py.Dataset):
                name = pbasename(obj.name).encode('utf-8')
                extract(outdir, obj.parent, name)
            else:
                return
        else:
            msg = "{pathname} not found in {fname}"
            print(msg.format(pathname=pathname, fname=fname))
            return


def _parser():
    """ Argument parser. """
    description = "Extracts HDF5 datasets to either GeoTiff or CSV."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--filename", required=True,
                        help="The file from which to extract datasets from.")
    parser.add_argument("--outdir", required=True,
                        help=("The output directory that will contain the "
                              "extracted datasets."))
    parser.add_argument("--pathname", default="/",
                        help=("The pathname to either a Dataset or a Group. "
                              "Default is '/', which is root level."))

    return parser


def main():
    """ Main execution. """
    parser = _parser()
    args = parser.parse_args()
    run(args.filename, args.outdir, args.pathname)
