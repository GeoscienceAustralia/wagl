#!/usr/bin/env python

"""
Various interpolation methods.
"""

from __future__ import absolute_import
from os.path import basename, splitext
from posixpath import join as ppjoin
import math
from scipy.interpolate import Rbf
import numpy as np
import h5py
from gaip.constants import DatasetName, Model, GroupName
from gaip.hdf5 import dataset_compression_kwargs
from gaip.hdf5 import write_h5_image
from gaip.hdf5 import read_h5_table

DEFAULT_ORIGIN = (0, 0)
DEFAULT_SHAPE = (8, 8)


def bilinear(shape, fUL, fUR, fLR, fLL, dtype=np.float64):
    """
    Bilinear interpolation of four scalar values.

    :param shape:
        Shape of interpolated grid (nrows, ncols).

    :param fUL:
        Data value at upper-left (NW) corner.

    :param fUR:
        Data value at upper-right (NE) corner.

    :param fLR:
        Data value at lower-right (SE) corner.

    :param fLL:
        Data value at lower-left (SW) corner.

    :param dtype:
        Data type (numpy I presume?).

    :return:
        Array of data values interpolated between corners.
    """

    s, t = [a.astype(dtype) for a in np.ogrid[0:shape[0], 0:shape[1]]]

    s /= (shape[0] - 1.0)
    t /= (shape[1] - 1.0)

    result = (s * (t * fLR + (1.0 - t) * fLL) + (1.0 - s) *
              (t * fUR + (1.0 - t) * fUL))

    return result


def indices(origin=DEFAULT_ORIGIN, shape=DEFAULT_SHAPE):
    """
    Generate corner indices for a grid block.

    :param origin:
        Block origin (2-tuple).

    :param shape:
        Block shape (2-tuple: nrows, ncols).

    :return:
        Corner indices: (xmin, xmax, ymin, ymax).
    """
    return (origin[0], origin[0] + shape[0] - 1,
            origin[1], origin[1] + shape[1] - 1)


def subdivide(origin=DEFAULT_ORIGIN, shape=DEFAULT_SHAPE):
    """
    Generate indices for grid sub-blocks.

    :param origin:
        Block origin (2-tuple).

    :param shape:
        Block shape (nrows, ncols).

    :return:
        Dictionary containing sub-block corner indices:
            { 'UL': <list of 2-tuples>,
              'UR': <list of 2-tuples>,
              'LL': <list of 2-tuples>,
              'LR': <list of 2-tuples> }
    """
    i0, ie, j0, je = indices(origin, shape)
    ic = origin[0] + shape[0] // 2
    jc = origin[1] + shape[1] // 2

    return {
        'UL': [(i0, j0), (i0, jc), (ic, j0), (ic, jc)],
        'LL': [(ic, j0), (ic, jc), (ie, j0), (ie, jc)],
        'UR': [(i0, jc), (i0, je), (ic, jc), (ic, je)],
        'LR': [(ic, jc), (ic, je), (ie, jc), (ie, je)],
    }


def interpolate_block(origin=DEFAULT_ORIGIN, shape=DEFAULT_SHAPE,
                      eval_func=None, grid=None):
    """
    Interpolate a grid block.

    :param origin:
        Block origin (2-tuple).

    :param shape:
        Block shape (nrows, ncols).

    :param eval_func:
        Evaluator function.
    :type eval_func:
        callable; accepts grid indices i, j and returns a scalar value.

    :param grid:
        Grid array.
    :type grid:
        :py:class:`numpy.array`.

    :return:
        Interpolated block array if grid argument is None. If grid argument
        is supplied its elements are modified in place and this function
        does not return a value.
    """
    i0, i1, j0, j1 = indices(origin, shape)

    fUL = eval_func(i0, j0)
    fLL = eval_func(i1, j0)
    fUR = eval_func(i0, j1)
    fLR = eval_func(i1, j1)

    if grid is None:
        return bilinear(shape, fUL, fUR, fLR, fLL)

    grid[i0:i1 + 1, j0:j1 + 1] = bilinear(shape, fUL, fUR, fLR, fLL)


def interpolate_grid(depth=0, origin=DEFAULT_ORIGIN, shape=DEFAULT_SHAPE,
                     eval_func=None, grid=None):
    """
    Interpolate a data grid.

    :param depth:
        Recursive bisection depth.
    :type depth:
        :py:class:`int`

    :param origin:
        Block origin,
    :type origin:
        :py:class:`tuple` of length 2.

    :param shape:
        Block shape.
    :type shape:
        :py:class:`tuple` of length 2 ``(nrows, ncols)``.

    :param eval_func:
        Evaluator function.
    :type eval_func:
        callable; accepts grid indices i, j and returns a scalar value.

    :param grid:
        Grid array.
    :type grid:
        :py:class:`numpy.array`.

    :todo:
        Move arguments ``eval_func`` and ``grid`` to positions 1 and 2, and remove
        defaults (and the check that they are not ``None`` at the top of the function
        body).
    """
    assert eval_func is not None
    assert grid is not None

    if depth == 0:
        interpolate_block(origin, shape, eval_func, grid)
    else:
        blocks = subdivide(origin, shape)
        for (kUL, kUR, kLL, kLR) in blocks.values():
            block_shape = (kLR[0] - kUL[0] + 1, kLR[1] - kUL[1] + 1)
            interpolate_grid(depth - 1, kUL, block_shape, eval_func, grid)


def fortran_bilinear_interpolate(cols, rows, locations, samples,
                                 row_start, row_end, row_centre):
    """
    Original NBAR interpolation scheme.
    Sheared 4-cell bilinear, implemented in fortran.
    """
    from gaip.__bilinear_interpolation import bilinear_interpolation as fortran

    assert len(samples) == 3*3
    assert len(locations) == len(samples)

    s1 = samples[[0, 1, 3, 4]]
    s2 = samples[[1, 2, 4, 5]]
    s3 = samples[[3, 4, 6, 7]]
    s4 = samples[[4, 5, 7, 8]]

    output = np.empty((rows, cols), dtype=np.float32)

    fortran(cols, rows, locations, s1, s2, s3, s4,
            row_start, row_end, row_centre, output.T)

    return output


def rbf_interpolate(cols, rows, locations, samples, *_,
                    chunking=True, kernel='gaussian'):
    """scipy radial basis function interpolation"""

    rbf = Rbf(locations[:, 1], locations[:, 0], samples - samples.mean(),
              function=kernel)

    if not chunking:
        return rbf(*np.mgrid[:cols, :rows]).astype(np.float32) + samples.mean()
    else:
        xchunks, ychunks = cols//100 + 1, rows//100 + 1
        raster = np.empty((rows, cols), dtype=np.float32)
        for y in np.array_split(np.arange(rows), ychunks):
            for x in np.array_split(np.arange(cols), xchunks):
                xy = np.meshgrid(x, y) # note sparse not supported by scipy
                r = rbf(*xy).astype(np.float32) + samples.mean()
                raster[y[0]:y[-1]+1, x[0]:x[-1]+1] = r
        return raster


def sheared_bilinear_interpolate(cols, rows, locations, samples,
                                 row_start, row_end, row_centre):
    """
    Generalisation of the original NBAR interpolation scheme
    """

    n = len(samples)
    grid_size = int(math.sqrt(n)) - 1

    assert (grid_size+1)**2 == n and not (grid_size & 1) and not grid_size % 1
    # Assume count of samples is 9 or 25, 49, 81.. (Grid size is 2, 4, 6, ..)

    vertex_shape = (grid_size+1,)*2
    locations = locations.reshape(vertex_shape+(2,))
    samples = samples.reshape(vertex_shape)

    once = lambda f: f() # only for code structure
    @once
    def lines():
        """Place parcel boundaries alongside track (by 1D linear interpolation)"""
        L = np.empty((grid_size+1, rows), dtype=np.uint64)

        middle_vertex = grid_size//2

        L[0] = row_start
        L[middle_vertex] = row_centre
        L[-1] = row_end

        for i in range(1, middle_vertex):
            L[i] = row_start + (row_centre - row_start) * (i/middle_vertex)
            L[i+middle_vertex] = row_centre + (row_end - row_centre) * (i/middle_vertex)

        return L.reshape(grid_size+1, rows, 1) # enable broadcast

    y, x = np.ogrid[:rows, :cols]

    @once
    def parcellation(masking=True):
        """Generate parcellation map"""
        # Would use e.g. a matplotlib poly drawing routine,
        # except if curved sides are required.

        zones = np.full((rows, cols), 0, dtype=np.int8)

        # first axis
        vlines = lines
        for line in vlines[1:-1]:
            zones += (x >= line)
            # note, needn't be concerned with how edges are zoned
            # because they will be masked out downstream

        # second axis
        hlines = locations[:, 0, 0] # y (row) component of left-edge samples
        for line in hlines[1:-1]:
            zones += (y >= line) * grid_size

        if masking:
            zones[(x < vlines[0]) | (x > vlines[-1]) |
                  (y < hlines[0]) | (y > hlines[-1])] = -1

        return zones

    def shear(i, j, both_sides=False):
        """Warp to straighten edges of trapezoid"""
        if not both_sides:
            return x - lines[0] # original style NBAR algorithm
        else:
            left = lines[j]
            width = lines[j+1] - left

            xx = x - left

            return xx / width # * width.mean() if matrix nearly singular

    def patch(i, j, x=x, with_shear=True):
        """bilinear cell"""
        vertices = locations[i:i+2, j:j+2].reshape(4, 2)
        values = samples[i:i+2, j:j+2].reshape(4)

        if with_shear: # then re-map coordinates
            x = shear(i, j)
            old = vertices
            vertices = old.astype(np.float32, copy=True)
            vertices[:, 1] = x[list(old.T)]

        matrix = np.ones((4, 4))
        matrix[:, 1:3] = vertices
        matrix[:, 3] = vertices[:, 0] * vertices[:, 1]

        a = np.linalg.solve(matrix, values) # determine coefficients

        return a[0] + a[1]*y + a[2]*x + a[3]*x*y # broadcast as raster

    result = np.full((rows, cols), np.nan, dtype=np.float32)
    for i in range(grid_size):
        for j in range(grid_size):
            subset = parcellation == i*grid_size + j
            result[subset] = patch(i, j)[subset]

    return result


def _interpolate(acq, factor, sat_sol_angles_fname, coefficients_fname,
                 ancillary_fname, out_fname, compression, y_tile, method):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(sat_sol_angles_fname, 'r') as sat_sol,\
        h5py.File(coefficients_fname, 'r') as coef,\
        h5py.File(ancillary_fname, 'r') as anc,\
        h5py.File(out_fname, 'w') as out_fid:

        grp1 = anc[GroupName.ancillary_group.value]
        grp2 = sat_sol[GroupName.sat_sol_group.value]
        grp3 = coef[GroupName.coefficients_group.value]
        interpolate(acq, factor, grp1, grp2, grp3, out_fid, compression,
                    y_tile, method)


def interpolate(acq, factor, ancillary_group, satellite_solar_group,
                coefficients_group, out_group=None, compression='lzf',
                y_tile=100, method=None):
    # TODO: more docstrings
    """Perform interpolation."""
    geobox = acq.gridded_geo_box()
    cols, rows = geobox.get_shape_xy()

    # read the relevant tables into DataFrames
    coordinator = read_h5_table(ancillary_group, DatasetName.coordinator.value)
    boxline = read_h5_table(satellite_solar_group, DatasetName.boxline.value)

    if factor in Model.nbar.factors:
        dataset_name = DatasetName.nbar_coefficients.value
    elif factor in Model.sbt.factors:
        dataset_name = DatasetName.sbt_coefficients.value
    else:
        msg = "Factor name not found in available factors: {}"
        raise ValueError(msg.format(Model.standard.factors))

    coefficients = read_h5_table(coefficients_group, dataset_name)

    coord = np.zeros((coordinator.shape[0], 2), dtype='int')
    map_x = coordinator.map_x.values
    map_y = coordinator.map_y.values
    coord[:, 1], coord[:, 0] = (map_x, map_y) * ~geobox.transform
    centre = boxline.bisection_index.values + 1
    start = boxline.start_index.values + 1
    end = boxline.end_index.values + 1

    band = acq.band_num
    band_records = coefficients.band_id == 'BAND {}'.format(band)
    samples = coefficients[factor][band_records].values

    if method is None:
        if len(samples) == 9:
            method = 'linear'
        else:
            method = 'shear'

    func_map = {'linear': fortran_bilinear_interpolate,
                'shear': sheared_bilinear_interpolate,
                'rbf': rbf_interpolate}
    assert method in func_map

    result = func_map[method](cols, rows, coord, samples, start, end, centre)

    # setup the output file/group as needed
    if out_group is None:
        fid = h5py.File('interpolated-coefficients.h5', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    group = fid.create_group(GroupName.interp_group.value)

    fmt = DatasetName.interpolation_fmt.value
    dset_name = fmt.format(factor=factor, band=band)
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))
    no_data = -999
    kwargs['fillvalue'] = no_data
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal(),
             'no_data_value': no_data,
             'interpolation_method': method}
    desc = ("Contains the interpolated result of factor {}"
            "for band {} from sensor {}.")
    attrs['Description'] = desc.format(factor, band, acq.satellite_name)
    write_h5_image(result, dset_name, group, attrs, **kwargs)

    if out_group is None:
        return fid


def link_interpolated_data(data, out_fname):
    """
    Links the individual interpolated results into a
    single file for easier access.
    """
    group_path = GroupName.interp_group.value
    for key in data:
        fname = data[key]
        with h5py.File(fname, 'r') as fid:
            dataset_names = list(fid[group_path].keys())

        with h5py.File(out_fname, 'a') as fid:
            for dname in dataset_names:
                dataset_name = ppjoin(group_path, dname)
                fid[dataset_name] = h5py.ExternalLink(fname, dataset_name)
