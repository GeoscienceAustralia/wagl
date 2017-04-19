#!/usr/bin/env python

import os
from os.path import join as pjoin, normpath, dirname, exists, basename
from posixpath import join as ppjoin
from functools import partial
import argparse
import h5py
import yaml

from gaip.data import write_img
from gaip.geobox import GriddedGeoBox
from gaip.hdf5 import read_table

IGNORE = ['crs_wkt', 'geotransform']


def convert_image(dataset, output_directory):
    """
    Converts a HDF5 `IMAGE` Class dataset to a compressed GeoTiff,
    with deflate zlevel 1 compression.
    Any attributes stored with the image will be written as dataset
    level metadata tags, and not band level tags.

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

    kwargs = {'fmt': 'GTiff',
              'geobox': geobox,
              'compress': 'deflate',
              'options': {'zlevel': 1},
              'tags': tags,
              'nodata': no_data}

    out_fname = pjoin(output_directory, normpath(dataset.name.strip('/')))
    out_fname = ''.join([out_fname, '.tif'])

    if not exists(dirname(out_fname)):
        os.makedirs(dirname(out_fname))

    write_img(dataset, out_fname, **kwargs)


def convert_table(group, dataset_name, output_directory):
    """
    Converts a HDF5 `TABLE` Class dataset to a csv file.
    Any attributes stored with the table will not be carried across
    with the csv.

    :param dataset_name:
        A `str` containing the pathname to the HDF5 `Table` Class
        dataset.

    :param output_directory:
        A filesystem path to the directory that will be the root
        directory for any tables extracted.

    :return:
        None, outputs are written directly to disk.
    """
    df = read_table(group, dataset_name)

    out_fname = pjoin(output_directory, normpath(dataset_name.strip('/')))
    out_fname = ''.join([out_fname, '.csv'])

    if not exists(dirname(out_fname)):
        os.makedirs(dirname(out_fname))

    df.to_csv(out_fname)


def convert_scalar(dataset, output_directory):
    """
    Converts a HDF5 scalar dataset to a yaml file.
    All attributes will be output.

    :param dataset:
        A HDF5 scalar dataset.

    :param output_directory:
        A filesystem path to the directory that will be the root
        directory for any images extracted.

    :return:
        None, outputs are written directly to disk.
    """
    tags = {k: v for k, v in dataset.attrs.items()}
    tags[basename(dataset.name)] = dataset[()]

    out_fname = pjoin(output_directory, normpath(dataset.name.strip('/')))
    out_fname = ''.join([out_fname, '.yaml'])

    if not exists(dirname(out_fname)):
        os.makedirs(dirname(out_fname))

    with open(out_fname, 'w') as src:
        yaml.dump(tags, src, default_flow_style=False)


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


def run(fname, outdir):
    """ Run dataset conversion tree. """
    # note: lower level h5py access is required in order to visit links
    with h5py.File(fname, 'r') as fid:
        root = h5py.h5g.open(fid.fid, b'/')
        root.links.visit(partial(extract, outdir, fid))


def _parser():
    """ Argument parser. """
    description = "Extracts HDF5 datasets to either GeoTiff or CSV."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--filename", required=True,
                        help="The file from which to extract datasets from.")
    parser.add_argument("--outdir", required=True,
                        help=("The output directory that will contain the "
                              "extracted datasets."))

    return parser


def main():
    """ Main execution. """
    parser = _parser()
    args = parser.parse_args()
    run(args.filename, args.outdir)
