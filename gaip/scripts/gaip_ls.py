#!/usr/bin/env python

"""
A recursive utility that lists all Groups and Datasets contained
within a HDF5.
"""

import argparse
import h5py
from gaip.hdf5 import h5ls


def run(fname):
    """ List the contents of a HDF5 file. """
    # note: lower level h5py access is required in order to visit links
    with h5py.File(fname, 'r') as fid:
        h5ls(fid)


def _parser():
    """ Argument parser. """
    description = "Lists the contents of a HDF5 file."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--filename", required=True,
                        help="The file from which to list the contents.")

    return parser


def main():
    """ Main execution. """
    parser = _parser()
    args = parser.parse_args()
    run(args.filename)
