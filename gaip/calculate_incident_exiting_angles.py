#!/usr/bin/env python

"""
Calculates 2D grids of incident, exiting and relative azimuthal angles.
"""

import numpy
import h5py

from eotools.geobox import GriddedGeoBox
from eotools import tiling
from gaip import as_array
from gaip import GriddedGeoBox
from gaip import exiting_angle
from gaip import incident_angle
from gaip import dataset_compression_kwargs
from gaip import attach_image_attributes


def _incident_angles_wrapper(satellite_solar_fname, slope_aspect_fname,
                             out_fname, x_tile=None, y_tile=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(satellite_solar_fname, 'r') as sat_sol,\
        h5py.File(slope_aspect_fname, 'r') as slp_asp:

        solar_zen_dset = sat_sol['solar-zenith-angles']
        solar_azi_dset = sat_sol['solar-azimuthi-angles']
        slope_dset = slp_asp['slope']
        aspect_dset = slp_asp['aspect']

        shape = solar_zen_dset.shape
        transform = solar_zen_dset.attrs['geotransform']
        origin = (transform[0], transform[3])
        crs = solar_zen_dset.attrs['crs_wkt']
        pixelsize = (abs(transform[1]), abs(transform[5]))
        geobox = GriddedGeoBox(shape, origin, pixelsize, crs)

        fid = incident_angles(solar_zen_dset, solar_azi_dset, slope_dset,
                              aspect_dset, geobox, out_fname, x_tile, y_tile)

    fid.close()

    return


def incident_angles(solar_zenith_dataset, solar_azimuth_dataset, slope_datset,
                    aspect_dataset, geobox, out_fname=None, compression='lzf',
                    x_tile=None, y_tile=None):
    """
    Calculates the incident angle and the azimuthal incident angle.

    :param solar_zenith_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the solar zenith
        angles when index/sliced.

    :param solar_azimuth_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the solar azimuth
        angles when index/sliced.

    :param slope_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the slope values
        when index/sliced.

    :param aspect_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the aspect angles
        when index/sliced.

    :param geobox:
        An instance of a GriddedGeoBox object.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset names will be as follows:

        * incident-angles
        * azimuthal-incident-angles

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param x_tile:
        Defines the tile size along the x-axis. Default is None which
        equates to all elements along the x-axis.

    :param y_tile:
        Defines the tile size along the y-axis. Default is None which
        equates to all elements along the y-axis.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    shape = geobox.get_shape_yx()
    rows, cols = shape
    crs = geobox.crs.ExportToWkt()

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('incident-angles.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))
    no_data = -999
    kwargs['shape'] = shape
    kwargs['fillvalue'] = no_data
    kwargs['dtype'] = 'float32'

    # output datasets
    incident_dset = fid.create_dataset('incident-angles', **kwargs)
    azi_inc_dset = fid.create_dataset('azimuthal-incident-angles', **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': crs,
             'geotransform': geobox.affine.to_gdal(),
             'no_data_value': no_data}
    desc = "The incident angles in degrees."
    attrs['Description'] = desc
    attach_image_attributes(incident_dset, attrs)

    desc = "The azimuthal incident angles in degrees."
    attrs['Description'] = desc
    attach_image_attributes(azi_inc_dset, attrs)

    # Initialise the tiling scheme for processing
    if x_tile is None:
        x_tile = cols
    if y_tile is None:
        y_tile = rows
    tiles = tiling.generate_tiles(cols, rows, x_tile, y_tile, generator=False)

    # Loop over each tile
    for tile in tiles:
        # Row and column start and end locations
        ystart = tile[0][0]
        xstart = tile[1][0]
        yend = tile[0][1]
        xend = tile[1][1]
        idx = (slice(ystart, yend, slice(xstart, xend)))

        # Tile size
        ysize = yend - ystart
        xsize = xend - xstart

        # Read the data for the current tile
        # Convert to required datatype and transpose
        sol_zen = as_array(solar_zenith_dataset[idx],
                           dtype=numpy.float32, transpose=True)
        sol_azi = as_array(solar_azimuth_dataset[idx],
                           dtype=numpy.float32, transpose=True)
        slope = as_array(slope_datset[idx],
                         dtype=numpy.float32, transpose=True)
        aspect = as_array(aspect_dataset[idx],
                          dtype=numpy.float32, transpose=True)

        # Initialise the work arrays
        incident = numpy.zeros((ysize, xsize), dtype='float32')
        azi_incident = numpy.zeros((ysize, xsize), dtype='float32')

        # Process the current tile
        incident_angle(xsize, ysize, sol_zen, sol_azi, slope, aspect,
                       incident.transpose(), azi_incident.transpose())

        # Write the current tile to disk
        incident_dset[idx] = incident
        azi_inc_dset[idx] = azi_incident

    fid.flush()
    return fid


def _exiting_angles_wrapper(satellite_solar_fname, slope_aspect_fname,
                            out_fname, x_tile=None, y_tile=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(satellite_solar_fname, 'r') as sat_sol,\
        h5py.File(slope_aspect_fname, 'r') as slp_asp:

        sat_view_dset = sat_sol['satellite-view-angles']
        sat_azi_dset = sat_sol['satellite-azimuth-angles']
        slope_dset = slp_asp['slope']
        aspect_dset = slp_asp['aspect']

        shape = sat_view_dset.shape
        transform = sat_view_dset.attrs['geotransform']
        origin = (transform[0], transform[3])
        crs = sat_view_dset.attrs['crs_wkt']
        pixelsize = (abs(transform[1]), abs(transform[5]))
        geobox = GriddedGeoBox(shape, origin, pixelsize, crs)

        fid = incident_angles(sat_view_dset, sat_azi_dset, slope_dset,
                              aspect_dset, geobox, out_fname, x_tile, y_tile)

    fid.close()

    return


def exiting_angles(satellite_view_dataset, satellite_azimuth_dataset,
                   slope_dataset, aspect_dataset, geobox, out_fname=None,
                   compression='lzf', x_tile=None, y_tile=None):
    """
    Calculates the exiting angle and the azimuthal exiting angle.

    :param satellite_view_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the satellite view
        angles when index/sliced.

    :param satellite_azimuth_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the satellite
        azimuth angles when index/sliced.

    :param slope_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the slope values
        when index/sliced.

    :param aspect_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the aspect angles
        when index/sliced.

    :param geobox:
        An instance of a GriddedGeoBox object.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset names will be as follows:

        * exiting-angles
        * azimuthal-exiting-angles

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param x_tile:
        Defines the tile size along the x-axis. Default is None which
        equates to all elements along the x-axis.

    :param y_tile:
        Defines the tile size along the y-axis. Default is None which
        equates to all elements along the y-axis.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    shape = geobox.get_shape_yx()
    rows, cols = shape
    crs = geobox.crs.ExportToWkt()

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('incident-angles.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))
    no_data = -999
    kwargs['shape'] = shape
    kwargs['fillvalue'] = no_data
    kwargs['dtype'] = 'float32'

    # output datasets
    exiting_dset = fid.create_dataset('exiting-angles', **kwargs)
    azi_exit_dset = fid.create_dataset('azimuthal-exiting-angles', **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': crs,
             'geotransform': geobox.affine.to_gdal(),
             'no_data_value': no_data}
    desc = "The exiting angles in degrees."
    attrs['Description'] = desc
    attach_image_attributes(exiting_dset, attrs)

    desc = "The azimuthal exiting angles in degrees."
    attrs['Description'] = desc
    attach_image_attributes(azi_exit_dset, attrs)

    # Initialise the tiling scheme for processing
    if x_tile is None:
        x_tile = cols
    if y_tile is None:
        y_tile = rows
    tiles = tiling.generate_tiles(cols, rows, x_tile, y_tile, generator=False)

    # Loop over each tile
    for tile in tiles:
        # Row and column start and end locations
        ystart = tile[0][0]
        xstart = tile[1][0]
        yend = tile[0][1]
        xend = tile[1][1]
        idx = (slice(ystart, yend, slice(xstart, xend)))

        # Tile size
        ysize = yend - ystart
        xsize = xend - xstart

        # Read the data for the current tile
        # Convert to required datatype and transpose
        sat_view = as_array(satellite_view_dataset[idx],
                            dtype=numpy.float32, transpose=True)
        sat_azi = as_array(satellite_azimuth_dataset[idx],
                           dtype=numpy.float32, transpose=True)
        slope = as_array(slope_dataset[idx],
                         dtype=numpy.float32, transpose=True)
        aspect = as_array(aspect_dataset[idx],
                          dtype=numpy.float32, transpose=True)

        # Initialise the work arrays
        exiting = numpy.zeros((ysize, xsize), dtype='float32')
        azi_exiting = numpy.zeros((ysize, xsize), dtype='float32')

        # Process the current tile
        exiting_angle(xsize, ysize, sat_view, sat_azi, slope, aspect,
                      exiting.transpose(), azi_exiting.transpose())

        # Write the current to disk
        exiting_dset[idx] = exiting
        azi_exit_dset[idx] = azi_exiting

    fid.flush()
    return fid


def _relative_azimuth_slope_wrapper(incident_angles_fname,
                                    exiting_angles_fname, out_fname,
                                    compression='lzf', x_tile=None,
                                    y_tile=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(incident_angles_fname, 'r') as inci_angles,\
        h5py.File(exiting_angles_fname, 'r') as exit_angles:

        azi_inci_dset = inci_angles['azimuthal-incident-angles']
        azi_exit_dset = exit_angles['azimuthal-exiting-angles']

        shape = azi_inci_dset.shape
        transform = azi_inci_dset.attrs['geotransform']
        origin = (transform[0], transform[3])
        crs = azi_inci_dset.attrs['crs_wkt']
        pixelsize = (abs(transform[1]), abs(transform[5]))
        geobox = GriddedGeoBox(shape, origin, pixelsize, crs)

        fid = relative_azimuth_slope(azi_inci_dset, azi_exit_dset, geobox,
                                     out_fname, x_tile, y_tile)

    fid.close()

    return


def relative_azimuth_slope(azimuth_incident_dataset,
                           azimuth_exiting_dataset, geobox, out_fname=None,
                           compression='lzf', x_tile=None, y_tile=None):
    """
    Calculates the relative azimuth angle on the slope surface.

    :param azimuth_incident_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the azimuthal
        incident angles when index/sliced.

    :param azimuth_exiting_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the azimuthal
        exiting angles when index/sliced.

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset names will be as follows:

        * relative-azimuth-angles

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param x_tile:
        Defines the tile size along the x-axis. Default is None which
        equates to all elements along the x-axis.

    :param y_tile:
        Defines the tile size along the y-axis. Default is None which
        equates to all elements along the y-axis.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    shape = geobox.get_shape_yx()
    rows, cols = shape
    crs = geobox.crs.ExportToWkt()

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('relative-azimuth-angles.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, geobox.x_size()))
    no_data = -999
    kwargs['shape'] = shape
    kwargs['fillvalue'] = no_data
    kwargs['dtype'] = 'float32'

    # output datasets
    out_dset = fid.create_dataset('relative-azimuth-angles', **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': crs,
             'geotransform': geobox.affine.to_gdal(),
             'no_data_value': no_data}
    desc = "The relative azimuth angles on the slope surface in degrees."
    attrs['Description'] = desc
    attach_image_attributes(out_dset, attrs)

    # Initialise the tiling scheme for processing
    if x_tile is None:
        x_tile = cols
    if y_tile is None:
        y_tile = rows
    tiles = tiling.generate_tiles(cols, rows, x_tile, y_tile, generator=False)

    # Loop over each tile
    for tile in tiles:
        # Row and column start and end locations
        ystart, yend = tile[0]
        xstart, xend = tile[1]
        idx = (slice(ystart, yend, slice(xstart, xend)))

        # Read the data for the current tile
        azi_inc = azimuth_incident_dataset[idx]
        azi_exi = azimuth_exiting_dataset[idx]

        # Process the tile
        rel_azi = azi_inc - azi_exi
        rel_azi[rel_azi <= -180.0] += 360.0
        rel_azi[rel_azi > 180.0] -= 360.0

        # Write the current tile to disk
        out_dset[idx] = rel_azi

    fid.flush()
    return fid
