#!/usr/bin/env python

"""
Calculates 2D grids of incident, exiting and relative azimuthal angles.
"""

from __future__ import absolute_import, print_function
import numpy
import h5py

from wagl.constants import DatasetName, GroupName
from wagl.geobox import GriddedGeoBox
from wagl.tiling import generate_tiles
from wagl.data import as_array
from wagl.hdf5 import H5CompressionFilter, attach_image_attributes
from wagl.__exiting_angle import exiting_angle
from wagl.__incident_angle import incident_angle


def _incident_exiting_angles(satellite_solar_fname, slope_aspect_fname,
                             out_fname, compression=H5CompressionFilter.LZF,
                             filter_opts=None, incident=True):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(satellite_solar_fname, 'r') as sat_sol,\
        h5py.File(slope_aspect_fname, 'r') as slp_asp,\
        h5py.File(out_fname, 'w') as out_fid:

        grp1 = sat_sol[GroupName.SAT_SOL_GROUP.value]
        grp2 = slp_asp[GroupName.SLP_ASP_GROUP.value]
        if incident:
            incident_angles(grp1, grp2, out_fid, compression, filter_opts)
        else:
            exiting_angles(grp1, grp2, out_fid, compression, filter_opts)


def incident_angles(satellite_solar_group, slope_aspect_group, out_group=None,
                    compression=H5CompressionFilter.LZF, filter_opts=None):
    """
    Calculates the incident angle and the azimuthal incident angle.

    :param satellite_solar_group:
        The root HDF5 `Group` that contains the solar zenith and
        solar azimuth datasets specified by the pathnames given by:

        * DatasetName.SOLAR_ZENITH
        * DatasetName.SOLAR_AZIMUTH
        
    :param slope_aspect_group:
        The root HDF5 `Group` that contains the slope and aspect
        datasets specified by the pathnames given by:

        * DatasetName.SLOPE
        * DatasetName.ASPECT

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be as follows:

        * DatasetName.INCIDENT
        * DatasetName.AZIMUTHAL_INCIDENT

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF 

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    # dataset arrays
    dname = DatasetName.SOLAR_ZENITH.value
    solar_zenith_dataset = satellite_solar_group[dname]
    dname = DatasetName.SOLAR_AZIMUTH.value
    solar_azimuth_dataset = satellite_solar_group[dname]
    slope_dataset = slope_aspect_group[DatasetName.SLOPE.value]
    aspect_dataset = slope_aspect_group[DatasetName.ASPECT.value]

    geobox = GriddedGeoBox.from_dataset(solar_zenith_dataset)
    shape = geobox.get_shape_yx()
    rows, cols = shape
    crs = geobox.crs.ExportToWkt()

    # Initialise the output files
    if out_group is None:
        fid = h5py.File('incident-angles.h5', 'w', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    if GroupName.INCIDENT_GROUP.value not in fid:
        fid.create_group(GroupName.INCIDENT_GROUP.value)

    if filter_opts is None:
        filter_opts = {}

    grp = fid[GroupName.INCIDENT_GROUP.value]
    tile_size = solar_zenith_dataset.chunks
    filter_opts['chunks'] = tile_size
    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    no_data = numpy.nan
    kwargs['shape'] = shape
    kwargs['fillvalue'] = no_data
    kwargs['dtype'] = 'float32'

    # output datasets
    dataset_name = DatasetName.INCIDENT.value
    incident_dset = grp.create_dataset(dataset_name, **kwargs)
    dataset_name = DatasetName.AZIMUTHAL_INCIDENT.value
    azi_inc_dset = grp.create_dataset(dataset_name, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': crs,
             'geotransform': geobox.transform.to_gdal(),
             'no_data_value': no_data}
    desc = "Contains the incident angles in degrees."
    attrs['description'] = desc
    attrs['alias'] = 'incident'
    attach_image_attributes(incident_dset, attrs)

    desc = "Contains the azimuthal incident angles in degrees."
    attrs['description'] = desc
    attrs['alias'] = 'azimuthal-incident'
    attach_image_attributes(azi_inc_dset, attrs)

    # process by tile
    for tile in generate_tiles(cols, rows, tile_size[1], tile_size[0]):
        # Tile size
        ysize = tile[0].stop - tile[0].start
        xsize = tile[1].stop - tile[1].start

        # Read the data for the current tile
        # Convert to required datatype and transpose
        sol_zen = as_array(solar_zenith_dataset[tile],
                           dtype=numpy.float32, transpose=True)
        sol_azi = as_array(solar_azimuth_dataset[tile],
                           dtype=numpy.float32, transpose=True)
        slope = as_array(slope_dataset[tile],
                         dtype=numpy.float32, transpose=True)
        aspect = as_array(aspect_dataset[tile],
                          dtype=numpy.float32, transpose=True)

        # Initialise the work arrays
        incident = numpy.zeros((ysize, xsize), dtype='float32')
        azi_incident = numpy.zeros((ysize, xsize), dtype='float32')

        # Process the current tile
        incident_angle(xsize, ysize, sol_zen, sol_azi, slope, aspect,
                       incident.transpose(), azi_incident.transpose())

        # Write the current tile to disk
        incident_dset[tile] = incident
        azi_inc_dset[tile] = azi_incident

    if out_group is None:
        return fid


def exiting_angles(satellite_solar_group, slope_aspect_group, out_group=None,
                   compression=H5CompressionFilter.LZF, filter_opts=None):
    """
    Calculates the exiting angle and the azimuthal exiting angle.

    :param satellite_solar_group:
        The root HDF5 `Group` that contains the satellite view and
        satellite azimuth datasets specified by the pathnames given by:

        * DatasetName.SATELLITE_VIEW
        * DatasetName.SATELLITE_AZIMUTH
        
    :param slope_aspect_group:
        The root HDF5 `Group` that contains the slope and aspect
        datasets specified by the pathnames given by:

        * DatasetName.SLOPE
        * DatasetName.ASPECT

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be as follows:

        * DatasetName.EXITING
        * DatasetName.AZIMUTHAL_EXITING

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF 

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    # dataset arrays
    dname = DatasetName.SATELLITE_VIEW.value
    satellite_view_dataset = satellite_solar_group[dname]
    dname = DatasetName.SATELLITE_AZIMUTH.value
    satellite_azimuth_dataset = satellite_solar_group[dname]
    slope_dataset = slope_aspect_group[DatasetName.SLOPE.value]
    aspect_dataset = slope_aspect_group[DatasetName.ASPECT.value]

    geobox = GriddedGeoBox.from_dataset(satellite_view_dataset)
    shape = geobox.get_shape_yx()
    rows, cols = shape
    crs = geobox.crs.ExportToWkt()

    # Initialise the output files
    if out_group is None:
        fid = h5py.File('exiting-angles.h5', 'w', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    if GroupName.EXITING_GROUP.value not in fid:
        fid.create_group(GroupName.EXITING_GROUP.value)

    if filter_opts is None:
        filter_opts = {}

    grp = fid[GroupName.EXITING_GROUP.value]
    tile_size = satellite_view_dataset.chunks
    filter_opts['chunks'] = tile_size
    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    no_data = numpy.nan
    kwargs['shape'] = shape
    kwargs['fillvalue'] = no_data
    kwargs['dtype'] = 'float32'

    # output datasets
    dataset_name = DatasetName.EXITING.value
    exiting_dset = grp.create_dataset(dataset_name, **kwargs)
    dataset_name = DatasetName.AZIMUTHAL_EXITING.value
    azi_exit_dset = grp.create_dataset(dataset_name, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': crs,
             'geotransform': geobox.transform.to_gdal(),
             'no_data_value': no_data}
    desc = "Contains the exiting angles in degrees."
    attrs['description'] = desc
    attrs['alias'] = 'exiting'
    attach_image_attributes(exiting_dset, attrs)

    desc = "Contains the azimuthal exiting angles in degrees."
    attrs['description'] = desc
    attrs['alias'] = 'azimuthal-exiting'
    attach_image_attributes(azi_exit_dset, attrs)

    # process by tile
    for tile in generate_tiles(cols, rows, tile_size[1], tile_size[0]):
        # Tile size
        ysize = tile[0].stop - tile[0].start
        xsize = tile[1].stop - tile[1].start

        # Read the data for the current tile
        # Convert to required datatype and transpose
        sat_view = as_array(satellite_view_dataset[tile],
                            dtype=numpy.float32, transpose=True)
        sat_azi = as_array(satellite_azimuth_dataset[tile],
                           dtype=numpy.float32, transpose=True)
        slope = as_array(slope_dataset[tile],
                         dtype=numpy.float32, transpose=True)
        aspect = as_array(aspect_dataset[tile],
                          dtype=numpy.float32, transpose=True)

        # Initialise the work arrays
        exiting = numpy.zeros((ysize, xsize), dtype='float32')
        azi_exiting = numpy.zeros((ysize, xsize), dtype='float32')

        # Process the current tile
        exiting_angle(xsize, ysize, sat_view, sat_azi, slope, aspect,
                      exiting.transpose(), azi_exiting.transpose())

        # Write the current to disk
        exiting_dset[tile] = exiting
        azi_exit_dset[tile] = azi_exiting

    if out_group is None:
        return fid


def _relative_azimuth_slope(incident_angles_fname, exiting_angles_fname,
                            out_fname, compression=H5CompressionFilter.LZF,
                            filter_opts=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(incident_angles_fname, 'r') as inci_fid,\
        h5py.File(exiting_angles_fname, 'r') as exit_fid,\
        h5py.File(out_fname, 'w') as out_fid:

        grp1 = inci_fid[GroupName.INCIDENT_GROUP.value]
        grp2 = exit_fid[GroupName.EXITING_GROUP.value]
        relative_azimuth_slope(grp1, grp2, out_fid, compression, filter_opts)


def relative_azimuth_slope(incident_angles_group, exiting_angles_group,
                           out_group=None, compression=H5CompressionFilter.LZF,
                           filter_opts=None):
    """
    Calculates the relative azimuth angle on the slope surface.

    :param incident_angles_group:
        The root HDF5 `Group` that contains the azimuthal incident
        angle dataset specified by the pathname given by:

        * DatasetName.AZIMUTHAL_INCIDENT

    :param exiting_angles_group:
        The root HDF5 `Group` that contains the azimuthal exiting
        angle dataset specified by the pathname given by:

        * DatasetName.AZIMUTHAL_EXITING

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be as follows:

        * DatasetName.RELATIVE_SLOPE

    :param compression:
        The compression filter to use.
        Default is H5CompressionFilter.LZF 

    :filter_opts:
        A dict of key value pairs available to the given configuration
        instance of H5CompressionFilter. For example
        H5CompressionFilter.LZF has the keywords *chunks* and *shuffle*
        available.
        Default is None, which will use the default settings for the
        chosen H5CompressionFilter instance.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    # dataset arrays
    dname = DatasetName.AZIMUTHAL_INCIDENT.value
    azimuth_incident_dataset = incident_angles_group[dname]
    dname = DatasetName.AZIMUTHAL_EXITING.value
    azimuth_exiting_dataset = exiting_angles_group[dname]

    geobox = GriddedGeoBox.from_dataset(azimuth_incident_dataset)
    shape = geobox.get_shape_yx()
    rows, cols = shape
    crs = geobox.crs.ExportToWkt()

    # Initialise the output files
    if out_group is None:
        fid = h5py.File('relative-azimuth-angles.h5', 'w', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    if GroupName.REL_SLP_GROUP.value not in fid:
        fid.create_group(GroupName.REL_SLP_GROUP.value)

    if filter_opts is None:
        filter_opts = {}

    grp = fid[GroupName.REL_SLP_GROUP.value]
    tile_size = azimuth_incident_dataset.chunks
    filter_opts['chunks'] = tile_size
    kwargs = compression.config(**filter_opts).dataset_compression_kwargs()
    no_data = numpy.nan
    kwargs['shape'] = shape
    kwargs['fillvalue'] = no_data
    kwargs['dtype'] = 'float32'

    # output datasets
    out_dset = grp.create_dataset(DatasetName.RELATIVE_SLOPE.value, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': crs,
             'geotransform': geobox.transform.to_gdal(),
             'no_data_value': no_data}
    desc = ("Contains the relative azimuth angles on the slope surface in "
            "degrees.")
    attrs['description'] = desc
    attrs['alias'] = 'relative-slope'
    attach_image_attributes(out_dset, attrs)

    # process by tile
    for tile in generate_tiles(cols, rows, tile_size[1], tile_size[0]):
        # Read the data for the current tile
        azi_inc = azimuth_incident_dataset[tile]
        azi_exi = azimuth_exiting_dataset[tile]

        # Process the tile
        rel_azi = azi_inc - azi_exi
        rel_azi[rel_azi <= -180.0] += 360.0
        rel_azi[rel_azi > 180.0] -= 360.0

        # Write the current tile to disk
        out_dset[tile] = rel_azi

    if out_group is None:
        return fid
