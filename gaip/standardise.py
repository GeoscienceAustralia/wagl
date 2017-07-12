#!/usr/bin/env python

from os.path import join as pjoin
from posixpath import join as ppjoin
import h5py

from gaip.acquisition import acquisitions
from gaip.longitude_latitude_arrays import create_lon_lat_grids


def card4l(level1, model, vertices, method, pixel_quality, landsea, out_fname):
    """
    CEOS Analysis Ready Data for Land.
    A workflow for producing standardised products that meet the
    CARD4L specification.

    TODO: modtran path, ancillary paths, tle path, compression, ytile
    """
    container = acquisitions(level1)

    with h5py.File(out_fname, 'w') as fid:
        for granule in container.granules:
            granule_root = container.get_root(granule=granule)

            for cgroup in container.groups:
                acqs = container.get_acquisitions(granule=granule, group=cgroup)
                root = container.get_root(granule=granule, group=cgroup)
                root_grp = fid.create_group(root)

                 # longitude and latitude
                 create_lon_lat_grids(acqs[0].gridded_geo_box(), root_group,
                                      compression=compression, y_tile=y_tile)
