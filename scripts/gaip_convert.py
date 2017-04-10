#!/usr/bin/env python

import os
from os.path import join as pjoin, normpath, dirname, exists
from posixpath import join as ppjoin
from functools import partial
import argparse
import h5py

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

    kwargs = {'fmt': 'GTiff',
              'geobox': geobox,
              'compress': 'deflate',
              'options': {'zlevel': 1},
              'tags': tags,
              'nodata': no_data}

    out_fname = pjoin(output_directory, normpath(dataset.name.strip('/')))
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
    df.to_csv(out_fname)


def extract(output_directory, group, name):
    """
    A simple utility that sends an object to the appropriate
    extraction utility.
    """
    dataset_name = ppjoin(group.name, name.decode('utf-8'))
    obj = group[dataset_name]

    if obj.attrs.get('CLASS') == 'IMAGE':
        convert_image(obj, output_directory)
    elif obj.attrs.get('CLASS') == 'TABLE':
        convert_table(group, dataset_name, output_directory)
    else:
        # TODO: scalars
        return None

def main(fname, output_directory):
    """
    Main execution.
    """
    # note: lower level h5py access is required in order to visit links
    with h5py.File(fname, 'r') as fid:
        root = h5py.h5g.open(fid.fid, b'/')
        root.links.visit(partial(extract, output_directory, fid))


if __name__ == '__main__':
    description = "Extracts HDF5 datasets to either GeoTiff or CSV."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--filename", required=True,
                        help="The file from which to extract datasets from.")
    parser.add_argument("--outdir", required=True,
                        help=("The output directory that will contain the "
                              "extracted datasets."))
    args = parser.parse_args()
    main(args.filename, args.outdir)
