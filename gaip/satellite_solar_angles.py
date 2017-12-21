#!/usr/bin/env python

"""
Satellite and Solar angle claculations over a 2D grid.
"""

from __future__ import absolute_import, print_function
import math
import ephem
import numpy as np
import h5py

from osgeo import osr
from gaip.constants import DatasetName, GroupName
from gaip.hdf5 import dataset_compression_kwargs, attach_image_attributes
from gaip.hdf5 import attach_table_attributes, attach_attributes
from gaip.tle import load_tle
from gaip.__sat_sol_angles import angle
from gaip.__satellite_model import set_satmod
from gaip.__track_time_info import set_times

CRS = "EPSG:4326"


def convert_to_lonlat(geobox, col_index, row_index):
    """
    Converts arrays of row and column indices into latitude and
    longitude (WGS84 datum).

    :param geobox:
        An instance of a GriddedGeoBox object.

    :param col_index:
        A 1D `NumPy` array of integers representing the column indices.

    :param row_index:
        A 1D `NumPy` array of integers representing the row indices.

    :return:
        2x 1D `NumPy` arrays, of type float64, containing the
        (longitude, latitude) coordinates.
    """
    # Define the TO_CRS for lon & lat outputs
    sr = osr.SpatialReference()
    sr.SetFromUserInput(CRS)

    lon = np.zeros(row_index.shape, dtype='float64')
    lat = np.zeros(row_index.shape, dtype='float64')

    for i, coord in enumerate(zip(col_index, row_index)):
        map_xy = geobox.convert_coordinates(coord)
        lonlat = geobox.transform_coordinates(map_xy, to_crs=sr)
        lon[i] = lonlat[0]
        lat[i] = lonlat[1]

    return lon, lat
    

def create_centreline_dataset(geobox, x, n, out_group):
    """
    Creates the centre line dataset.

    :param geobox:
        An instance of a GriddedGeoBox object.

    :param x:
        A 1D np array of type int with the same shape as
        `geobox.shape[0]`.
        Details the column number starting at 0.

    :param n:
        A 1D np array of type int with the same shape as x.
        Details whether or not the track point coordinate is
        averaged.

    :param out_group:
        A writeable HDF5 `Group` object.

    :return:
        None, results are written into the H5Group defined by
        out_group.
    """
    # centreline
    # if more than one pixel in a line was a track point the coordinates
    # are averaged
    wh = n > 1.5
    x[wh] = x[wh] / n[wh]

    # check whether there is no centre pixel in the line. It is assumed that
    # at least the adjacent lines have pixel
    wh = n < 0.5
    temp = x[0:2].copy()
    x[wh] = np.roll(x, 1)[wh]
    # account for first element potentially being changed with the
    # last element
    if wh[0]:
        x[0] = temp[1]

    # convert X centre points to integers (basically array co-ordinates)
    # and correct for FORTRAN offset
    x = np.rint(x) - 1

    rows, _ = geobox.shape
    y = np.arange(rows)

    dtype = np.dtype([('row_index', 'int64'), ('col_index', 'int64'),
                      ('n_pixels', 'float'), ('latitude', 'float64'),
                      ('longitude', 'float64')])
    data = np.zeros(rows, dtype=dtype)
    lon, lat = convert_to_lonlat(geobox, x, y)

    data['row_index'] = y
    data['col_index'] = x
    data['n_pixels'] = n
    data['latitude'] = lat
    data['longitude'] = lon

    kwargs = dataset_compression_kwargs()
    dname = DatasetName.centreline.value
    cent_dset = out_group.create_dataset(dname, data=data, **kwargs)
    desc = ("Contains the array, latitude and longitude coordinates of the "
            "satellite track path.")
    attrs = {'description': desc,
             'array_coordinate_offset': 0}
    attach_table_attributes(cent_dset, title='Centreline', attrs=attrs)


def first_and_last(array):
    """
    Utility to find indices of first and last true value in an array.

    Example:
                            0 1 2 3 4 5 6 7 8 9
                                |       |
        >>> first_and_last([0,0,1,1,0,1,1,0,0,0])
        (2, 6)
        >>> first_and_last([0,0,0]) # if no element of truth:
        (-1, -1)
        >>> first_and_last([1])
        (0, 0)
    """
    i, = np.nonzero(array) # assume array only has one dimension to unpack
    return (i[0], i[-1]) if len(i) else (-1, -1)


def asymetric_linspace(start, stop, num, midpoint):
    """
    Utility like numpy.linspace but with custom midpoint.

        >>> assymetric_linspace(start=10, stop=20, num=5, midpoint=18)
        [10, 14, 18, 19, 20]
    """
    front = np.linspace(start, midpoint, num//2, endpoint=False, dtype='int64')
    back = np.linspace(midpoint, stop, num//2 + 1, dtype='int64')
    return list(front)+list(back)


def swathe_edges(threshold, array):
    """
    Find left and right edges of swathe.

    Takes raster array, compares values against threshold, and returns a pair
    of vectors which represent the indices of the first and last pixel in
    each row of the array.
    """
    start = np.empty(array.shape[0], dtype='int')
    end = np.empty(array.shape[0], dtype='int')
    for i, row in enumerate(array):
        start[i], end[i] = first_and_last(row <= threshold)
    return start, end


def create_boxline(geobox, view_angle_dataset, centreline_dataset, out_group,
                   max_angle=9.0):
    """
    Creates the boxline (satellite track bi-section) dataset.

    :param geobox:
        An instance of a GriddedGeoBox object.

    :param view_angle_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the satellite view
        angles when index/sliced.

    :param centreline_dataset:
        The dataset created by the create_centreline function.

    :param out_group:
        A writeable HDF5 `Group` object.

    :param max_angle:
        The maximum viewing angle. Default is 9.0 degrees.

    :return:
        None, results are written into the H5Group defined by
        out_group.
    """
    rows, _ = view_angle_dataset.shape
    
    # calculate the column start and end indices
    # (for filtering out pixels of the ortho' array where no observations
    # are expected because the sensor look-angle would be too peripheral.)
    # TODO: similar filtering for pixels where the line acquisition time would
    # be outside of the scene aquisition window.
    istart, iend = swathe_edges(max_angle, view_angle_dataset)

    row_index = np.arange(rows)
    col_index = centreline_dataset['col_index'][:]

    # record curves for parcellation (of the raster into interpolation cells)
    boxline_dtype = np.dtype([('row_index', 'int64'),
                              ('bisection_index', 'int64'),
                              ('npoints', 'int64'),
                              ('start_index', 'int64'),
                              ('end_index', 'int64'),
                              ('bisection_longitude', 'float64'),
                              ('bisection_latitude', 'float64'),
                              ('start_longitude', 'float64'),
                              ('start_latitude', 'float64'),
                              ('end_longitude', 'float64'),
                              ('end_latitude', 'float64')])
    boxline = np.empty(rows, dtype=boxline_dtype)
    boxline['row_index'] = row_index
    boxline['bisection_index'] = col_index
    boxline['npoints'] = centreline_dataset['n_pixels'][:]
    boxline['start_index'] = istart
    boxline['end_index'] = iend

    # lon/lat conversions
    lon, lat = convert_to_lonlat(geobox, col_index, row_index)
    boxline['bisection_longitude'] = lon
    boxline['bisection_latitude'] = lat

    lon, lat = convert_to_lonlat(geobox, istart, row_index)
    boxline['start_longitude'] = lon
    boxline['start_latitude'] = lat

    lon, lat = convert_to_lonlat(geobox, iend, row_index)
    boxline['end_longitude'] = lon
    boxline['end_latitude'] = lat

    kwargs = dataset_compression_kwargs()
    desc = ("Contains the bi-section, column start and column end array "
            "coordinates.")
    attrs = {'description': desc,
             'array_coordinate_offset': 0}
    dname = DatasetName.boxline.value
    box_dset = out_group.create_dataset(dname, data=boxline, **kwargs)
    attach_table_attributes(box_dset, title='Boxline', attrs=attrs)


def create_vertices(acquisition, boxline_dataset, vertices=(3, 3)):
    """
    Defines the point locations, and the number of points, across a
    spatial grid.

    :param acquisition:
        An instance of an `Acquisition` object.

    :param boxline:
        The dataset containing the bi-section (satellite track)
        coordinates. The datatype should be the same as that returned
        by the `satellite_solar_angles.create_boxline` function.
    :type boxline:
        [('row_index', 'int64'), ('bisection_index', 'int64'),
         ('npoints', 'int64'), ('start_index', 'int64'),
         ('end_index', 'int64')]

    :param vertices:
        An integer 2-tuple indicating the number of rows and columns
        of sample-locations ("coordinator") to produce.
        The vertex columns should be an odd number.
        Default is (3, 3).

    :return:
        A `NumPy` dataset of the following datatype:

        * np.dtype([('row_index', 'int64'),
                    ('col_index', 'int64'),
                    ('latitude', 'float64'),
                    ('longitude', 'float64'),
                    ('map_y', 'int64'),
                    ('map_x', 'int64')])
    """
    geobox = acquisition.gridded_geo_box()
    rows, cols = geobox.shape

    if rows < vertices[0] | cols < vertices[1]:
        msg = ("Vertices must be >= to the acquisition dimensions! "
               "Acquisition dimensions: {}, Vertices: {}")
        raise ValueError(msg.format((rows, cols), vertices))

    xcentre = boxline_dataset['bisection_index']
    npoints = boxline_dataset['npoints']
    istart = boxline_dataset['start_index']
    iend = boxline_dataset['end_index']

    # special handling if the satellite track does not intersect both ends
    # of the raster
    track_end_rows = set(first_and_last(npoints))
    partial_track = track_end_rows - {0, rows-1, -1}
    mid_row = rows // 2
    if -1 in track_end_rows: # track doesn't intersect raster
        mid_col = cols // 2
    elif partial_track: # track intersects only part of raster
        mid_col = ({xcentre[0], xcentre[-1]} - {0, cols-1, -1}).pop()
        mid_row = partial_track.pop()
    else: # track fully available for deference
        mid_col = None

    # Note, assumes that if track intersects two rows then it also
    # intersects all intervening rows.

    grid_rows = asymetric_linspace(0, rows-1, vertices[0], midpoint=mid_row)

    locations = np.empty((vertices[0], vertices[1], 2), dtype='int64')
    for ig, ir in enumerate(grid_rows): # row indices for sample-grid & raster
        grid_line = asymetric_linspace(
            istart[ir], iend[ir], vertices[1], mid_col or xcentre[ir])
        locations[ig, :, 0] = ir
        locations[ig, :, 1] = grid_line
    locations = locations.reshape(vertices[0] * vertices[1], 2)

    # custom datatype for coordinator
    coordinator_dtype = np.dtype([('row_index', 'int64'),
                                  ('col_index', 'int64'),
                                  ('latitude', 'float64'),
                                  ('longitude', 'float64'),
                                  ('map_y', 'int64'),
                                  ('map_x', 'int64')])
    coordinator = np.empty(locations.shape[0], dtype=coordinator_dtype)
    coordinator['row_index'] = locations[:, 0]
    coordinator['col_index'] = locations[:, 1]

    map_xy = (locations[:, 1], locations[:, 0]) * geobox.transform
    coordinator['map_y'] = map_xy[1]
    coordinator['map_x'] = map_xy[0]

    lon, lat = convert_to_lonlat(geobox, locations[:, 1], locations[:, 0])
    coordinator['latitude'] = lat
    coordinator['longitude'] = lon

    return coordinator


def calculate_julian_century(datetime):
    """
    Given a datetime object return the julian century from the 2000 epoch.

    :param datetime:
        A datetime object containing the date to be converted to a
        Julian centuries since 2000/01/01 12:00.

    :return:
        A floating point value representing the Julian centuries since
        the 2000/01/01 12:00 epoch.
    """

    # Convert the scene timestamp to a julian date
    d = ephem.date(datetime)
    jdate = ephem.julian_date(d)

    # Get the J2000 epoch
    epoch = ephem.date((2000, 1, 1, 12.00))
    j2_epoch = ephem.julian_date(epoch)

    # Note:
    # This differes from online sources such as
    # http://www.pietro.org/astro_util_staticdemo/FDetailDateConversions.htm
    # http://en.wikipedia.org/wiki/Equinox_(celestial_coordinates)
    # which use:
    # 2000 + (jdate - epoch) / 365.25
    century = (jdate - j2_epoch) / 36525

    return century


def setup_spheroid(proj_wkt):
    """
    Given a WKT projection string, determine the spheroid paramaters
    that will be used in calcultating the angle grids.

    :param proj_wkt:
        A string containing valid WKT projection information.

    :return:
        A floating point np array of 4 elements containing
        spheroidal paramaters.

            * Index 0 contains the spheroid Major Axis.
            * Index 1 contains the spheroid Inverse Flattening.
            * Index 2 contains the spheroid Squared Eccentricity.
            * Index 3 contains the Earth rotational angular velocity in
              radians/second.

        Also a np dataset of the following datatype:

            * dtype = [('semi_major_axis', 'float64'),
                       ('inverse_flattening', 'float64'),
                       ('eccentricity_squared', 'float64'),
                       ('earth_rotational_angular_velocity', 'float64')]
    """
    dtype = np.dtype([('semi_major_axis', 'float64'),
                      ('inverse_flattening', 'float64'),
                      ('eccentricity_squared', 'float64'),
                      ('earth_rotational_angular_velocity', 'float64')])
    dset = np.zeros(1, dtype=dtype)

    # Define the spatial reference
    sr = osr.SpatialReference()
    sr.ImportFromWkt(proj_wkt)

    # Spheroid major axis
    dset['semi_major_axis'] = sr.GetSemiMajor()

    # Inverse flattening
    dset['inverse_flattening'] = sr.GetInvFlattening()

    # Eccentricity squared
    dset['eccentricity_squared'] = 1.0 - (1.0 - 1.0 /
                                          dset['inverse_flattening'])**2

    # Earth rotational angular velocity rad/sec
    # Other sources such as:
    # http://www.oosa.unvienna.org/pdf/icg/2012/template/WGS_84.pdf
    # state 0.000072921150 as the mean value
    dset['earth_rotational_angular_velocity'] = 0.000072722052

    return np.array(dset.tolist()).squeeze(), dset


def setup_orbital_elements(acquisition, tle_path):
    """
    Given an ephemeral object and a datetime object, calculate the
    satellite orbital paramaters used for calculating the angle grids.

    :param acquisition:
        An `Acquisition` object.

    :param tle_path:
        A `str` to the directory containing the Two Line Element data.

    :return:
        A floating point np array of 3 elements containing the
        satellite ephemeral bodies orbital paramaters.

            * Index 0 contains the obrital inclination in degrees.
            * Index 1 contains the semi major raidus in metres.
            * Index 2 contains the angular velocity in radians/sec^1.

        Also a np dataset of the following datatype:

            * dtype = [('orbital_inclination', 'float64'),
                       ('semi_major_radius', 'float64'),
                       ('angular_velocity', 'float64')]
    """
    dtype = np.dtype([('orbital_inclination', 'float64'),
                      ('semi_major_radius', 'float64'),
                      ('angular_velocity', 'float64')])
    dset = np.zeros(1, dtype=dtype)

    ephemeral = load_tle(acquisition, tle_path)

    # If we have None, then no suitable TLE was found, so use values gathered
    # by the acquisition object
    if ephemeral is None:
        # orbital inclination (degrees)
        dset['orbital_inclination'] = math.degrees(acquisition.inclination)
        # semi_major radius (m)
        dset['semi_major_radius'] = acquisition.semi_major_axis
        # angular velocity (rad sec-1)
        dset['angular_velocity'] = acquisition.omega
    else:
        ephemeral.compute(acquisition.acquisition_datetime)
        pi = np.pi
        n = ephemeral._n  # number or orbits per day
        s = 24*60*60  # Seconds in a day
        mu = 398600441800000.0  # Earth Gravitational parameter m^3s^-2

        # orbital inclination (degrees)
        dset['orbital_inclination'] = np.rad2deg(ephemeral._inc)

        # semi_major radius (m)
        # http://smallsats.org/2012/12/06/two-line-element-set-tle/
        dset['semi_major_radius'] = (mu / (2*pi*n/s)**2)**(1/3)

        # angular velocity (rad sec-1)
        dset['angular_velocity'] = (2*pi*n) / s

    return np.array(dset.tolist()).squeeze(), dset


def setup_smodel(centre_lon, centre_lat, spheroid, orbital_elements):
    """
    Setup the satellite model.
    A wrapper routine for the `set_satmod` Fortran module built via
    ``F2Py``.

    :param centre_lon:
        The longitude of the scene centre.

    :param centre_lat:
        The lattitude of the scene centre.

    :param spheroid:
        A 4 element floating point array containing the Earth
        spheroidal paramaters.

            * Index 0 contains the spheroid Major Axis.
            * Index 1 contains the spheroid Inverse Flattening.
            * Index 2 contains the spheroid Squared Eccentricity.
            * Index 3 contains the Earth rotational angular velocity in
              radians/second.

    :param orbital_elements:
        A 3 element floating point array containing the satellite
        orbital elements.

            * Index 0 contains the obrital inclination in degrees.
            * Index 1 contains the semi major raidus in metres.
            * Index 2 contains the angular velocity in radians/sec^1.

    :return:
        A floating point np array of 12 elements containing the
        satellite model paramaters.

            * Index 0 contains phi0.
            * Index 1 contains phi0_p.
            * Index 2 contains rho0.
            * Index 3 contains t0.
            * Index 4 contains lam0.
            * Index 5 contains gamm0.
            * Index 6 contains beta0.
            * Index 7 contains rotn0.
            * Index 8 contains hxy0.
            * Index 9 contains N0.
            * Index 10 contains H0.
            * Index 11 contains th_ratio0.

        Also np dataset of the following datatype:

            * dtype = [('phi0', 'f8'), ('phi0_p', 'f8'), ('rho0', 'f8'),
                       ('t0', 'f8'), ('lam0', 'f8'), ('gamm0', 'f8'),
                       ('beta0', 'f8'), ('rotn0', 'f8'), ('hxy0', 'f8'),
                       ('N0', 'f8'), ('H0', 'f8'), ('th_ratio0', 'f8')]
    """
    smodel, _ = set_satmod(centre_lon, centre_lat, spheroid, orbital_elements)

    columns = ['phi0', 'phi0_p', 'rho0', 't0', 'lam0', 'gamm0', 'beta0',
               'rotn0', 'hxy0', 'N0', 'H0', 'th_ratio0']
    dtype = np.dtype([(col, 'float64') for col in columns])
    smodel_dset = np.zeros(1, dtype=dtype)
    smodel_dset[0] = smodel

    return smodel, smodel_dset


def setup_times(ymin, ymax, spheroid, orbital_elements, smodel, ntpoints=12):
    """
    Setup the satellite track times.
    A wrapper routine for the ``set_times`` Fortran module built via
    ``F2Py``.

    :param ymin:
        The minimum lattitude in the array extent.

    :param ymax:
        The maximum lattitude in the array extent.

    :param spheroid:
        A 4 element floating point array containing the Earth
        spheroidal paramaters.

            * Index 0 contains the spheroid Major Axis.
            * Index 1 contains the spheroid Inverse Flattening.
            * Index 2 contains the spheroid Squared Eccentricity.
            * Index 3 contains the Earth rotational angular velocity in
              radians/second.

    :param orbital_elements:
        A 3 element floating point array containing the satellite
        orbital elements.

            * Index 0 contains the obrital inclination in degrees.
            * Index 1 contains the semi major raidus in metres.
            * Index 2 contains the angular velocity in radians/sec^1.

    :param smodel:
        A floating point np array of 12 elements containing the
        satellite model paramaters:

            * Index 0 contains phi0.
            * Index 1 contains phi0_p.
            * Index 2 contains rho0.
            * Index 3 contains t0.
            * Index 4 contains lam0.
            * Index 5 contains gamm0.
            * Index 6 contains beta0.
            * Index 7 contains rotn0.
            * Index 8 contains hxy0.
            * Index 9 contains N0.
            * Index 10 contains H0.
            * Index 11 contains th_ratio0.

    :param ntpoints:
        The number of time sample points to be calculated along the
        satellite track. Default is 12.

    :return:
        A floating point np array of [ntpoints,8] containing the
        satellite track times and other information.

            * Index 0 t.
            * Index 1 rho.
            * Index 2 phi_p.
            * Index 3 lam.
            * Index 4 beta.
            * Index 5 hxy.
            * Index 6 mj.
            * Index 7 skew.

        Also a np dataset of the datatype:

            * dtype = [('t', 'f8'), ('rho', 'f8'), ('phi_p', 'f8'),
                       ('lam', 'f8'), ('beta', 'f8'), ('hxy', 'f8'),
                       ('mj', 'f8'), ('skew', 'f8')]
    """
    track, _ = set_times(ymin, ymax, ntpoints, spheroid, orbital_elements,
                         smodel)

    columns = ['t', 'rho', 'phi_p', 'lam', 'beta', 'hxy', 'mj', 'skew']
    dtype = np.dtype([(col, 'float64') for col in columns])
    track_dset = np.zeros(ntpoints, dtype=dtype)
    track_dset['t'] = track[:, 0]
    track_dset['rho'] = track[:, 1]
    track_dset['phi_p'] = track[:, 2]
    track_dset['lam'] = track[:, 3]
    track_dset['beta'] = track[:, 4]
    track_dset['hxy'] = track[:, 5]
    track_dset['mj'] = track[:, 6]
    track_dset['skew'] = track[:, 7]

    return track, track_dset


def _store_parameter_settings(fid, spheriod, orbital_elements,
                              satellite_model, satellite_track, params):
    """
    An internal function for storing the parameter settings for the
    calculate_angles workflow.
    """
    group = fid.create_group('PARAMETERS')

    # generic parameters
    dname = DatasetName.generic.value
    dset = group.create_dataset(dname, 'GENERIC PARAMETERS')
    attach_attributes(dset, params)

    # sheroid
    desc = "The spheroid used in the satellite and solar angles calculation."
    attrs = {'description': desc}
    dname = DatasetName.spheroid.value
    sph_dset = group.create_dataset(dname, data=spheriod)
    attach_table_attributes(sph_dset, title='Spheroid', attrs=attrs)

    # orbital elements
    desc = ("The satellite orbital parameters used in the satellite and "
            "solar angles calculation.")
    attrs = {'description': desc}
    dname = DatasetName.orbital_elements.value
    orb_dset = group.create_dataset(dname, data=orbital_elements)
    attach_table_attributes(orb_dset, title='Orbital Elements', attrs=attrs)

    # satellite model
    desc = ("The satellite model used in the satellite and solar angles "
            "calculation.")
    attrs = {'description': desc}
    dname = DatasetName.satellite_model.value
    sat_dset = group.create_dataset(dname, data=satellite_model)
    attach_table_attributes(sat_dset, title='Satellite Model', attrs=attrs)

    # satellite track
    desc = ("The satellite track information used in the satellite and solar "
            "angles calculation.")
    attrs = {'description': desc}
    dname = DatasetName.satellite_track.value
    track_dset = group.create_dataset(dname, data=satellite_track)
    attach_table_attributes(track_dset, title='Satellite Track', attrs=attrs)


def _calculate_angles(acquisition, lon_lat_fname, out_fname=None,
                      compression='lzf', tle_path=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    with h5py.File(lon_lat_fname, 'r') as lon_lat_fid,\
        h5py.File(out_fname, 'w') as fid:
        lon_lat_grp = lon_lat_fid[GroupName.lon_lat_group.value]
        calculate_angles(acquisition, lon_lat_grp, fid, compression, tle_path)


def calculate_angles(acquisition, lon_lat_group, out_group=None,
                     compression='lzf', tle_path=None):
    """
    Calculate the satellite view, satellite azimuth, solar zenith,
    solar azimuth, and relative aziumth angle grids, as well as the
    time grid. All grids are output as float32 ENVI files.
    A wrapper routine for the ``angle_all`` Fortran module built via
    ``F2Py``.

    :param acquisition:
        An instance of an `Acquisition` object.

    :param lon_lat_group::
        The root HDF5 `Group` that contains the longitude and
        latitude datasets.
        The dataset pathnames are given by:

        * DatasetName.lon
        * DatasetName.lat

    :param out_group:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver. Otherwise,
        a writeable HDF5 `Group` object.

        The dataset names will be as follows:

        * DatasetName.satellite_view
        * DatasetName.satellite_azimuth
        * DatasetName.solar_zenith
        * DatasetName.solar_azimuth
        * DatasetName.relative_azimuth
        * DatasetName.acquisition_time
        * DatasetName.centreline
        * DatasetName.boxline
        * DatasetName.spheroid
        * DatasetName.orbital_elements
        * DatasetName.satellite_model
        * DatasetName.satellite_track

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param tle_path:
        A `str` to the directory containing the Two Line Element data.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    acq = acquisition
    century = calculate_julian_century(acq.acquisition_datetime)
    geobox = acq.gridded_geo_box()

    # longitude and latitude datasets
    longitude = lon_lat_group[DatasetName.lon.value]
    latitude = lon_lat_group[DatasetName.lat.value]

    # Min and Max lat extents
    # This method should handle northern and southern hemispheres
    # include a 1 degree buffer
    min_lat = min(min(geobox.ul_lonlat[1], geobox.ur_lonlat[1]),
                  min(geobox.ll_lonlat[1], geobox.lr_lonlat[1])) - 1
    max_lat = max(max(geobox.ul_lonlat[1], geobox.ur_lonlat[1]),
                  max(geobox.ll_lonlat[1], geobox.lr_lonlat[1])) + 1

    # Get the lat/lon of the scene centre
    # check if we have a file with GPS satellite track points
    # which can be used for cases of image granules/tiles, eg Sentinel-2A
    if acq.gps_file:
        points = acq.read_gps_file()
        subs = points[(points.latitude >= min_lat) &
                      (points.latitude <= max_lat)]
        idx = subs.shape[0] // 2 - 1
        centre_xy = (subs.iloc[idx].longitude, subs.iloc[idx].latitude)
    else:
        centre_xy = geobox.centre_lonlat

    # Get the earth spheroidal paramaters
    spheroid = setup_spheroid(geobox.crs.ExportToWkt())

    # Get the satellite orbital elements
    orbital_elements = setup_orbital_elements(acq, tle_path)

    # Get the satellite model paramaters
    smodel = setup_smodel(centre_xy[0], centre_xy[1], spheroid[0],
                          orbital_elements[0])

    # Get the times and satellite track information
    track = setup_times(min_lat, max_lat, spheroid[0], orbital_elements[0],
                        smodel[0])

    # Initialise the output files
    if out_group is None:
        fid = h5py.File('satellite-solar-angles.h5', driver='core',
                        backing_store=False)
    else:
        fid = out_group

    if GroupName.sat_sol_group.value not in fid:
        fid.create_group(GroupName.sat_sol_group.value)

    grp = fid[GroupName.sat_sol_group.value]

    # store the parameter settings used with the satellite and solar angles
    # function
    params = {'dimensions': (acq.lines, acq.samples),
              'lines': acq.lines,
              'samples': acq.samples,
              'century': century,
              'decimal_hour': acq.decimal_hour(),
              'acquisition_datetime': acq.acquisition_datetime,
              'centre_longitude_latitude': centre_xy,
              'minimum_latiude': min_lat,
              'maximum_latiude': max_lat,
              'latitude_buffer': 1.0,
              'max_view_angle': acq.maximum_view_angle}
    _store_parameter_settings(grp, spheroid[1], orbital_elements[1],
                              smodel[1], track[1], params)

    out_dtype = 'float32'
    no_data = -999
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=acq.tile_size)
    kwargs['shape'] = (acquisition.lines, acquisition.samples)
    kwargs['fillvalue'] = no_data
    kwargs['dtype'] = out_dtype

    sat_v_ds = grp.create_dataset(DatasetName.satellite_view.value, **kwargs)
    sat_az_ds = grp.create_dataset(DatasetName.satellite_azimuth.value,
                                   **kwargs)
    sol_z_ds = grp.create_dataset(DatasetName.solar_zenith.value, **kwargs)
    sol_az_ds = grp.create_dataset(DatasetName.solar_azimuth.value, **kwargs)
    rel_az_ds = grp.create_dataset(DatasetName.relative_azimuth.value,
                                   **kwargs)
    time_ds = grp.create_dataset(DatasetName.acquisition_time.value, **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.transform.to_gdal(),
             'no_data_value': no_data}
    desc = "Contains the satellite viewing angle in degrees."
    attrs['description'] = desc
    attach_image_attributes(sat_v_ds, attrs)

    desc = "Contains the satellite azimuth angle in degrees."
    attrs['description'] = desc
    attach_image_attributes(sat_az_ds, attrs)

    desc = "Contains the solar zenith angle in degrees."
    attrs['description'] = desc
    attach_image_attributes(sol_z_ds, attrs)

    desc = "Contains the solar azimuth angle in degrees."
    attrs['description'] = desc
    attach_image_attributes(sol_az_ds, attrs)

    desc = "Contains the relative azimuth angle in degrees."
    attrs['description'] = desc
    attach_image_attributes(rel_az_ds, attrs)

    desc = ("Contains the satellite acquisition time grid in seconds before "
            "and after the scene acquisition datetime.")
    attrs['description'] = desc
    attach_image_attributes(time_ds, attrs)

    # Initialise centre line variables
    x_cent = np.zeros((acq.lines), dtype=out_dtype)
    n_cent = np.zeros((acq.lines), dtype=out_dtype)

    # TODO: fix so that we can process in a way that doesn't require
    #       all columns to parsed through
    for tile in acq.tiles():
        idx = (slice(tile[0][0], tile[0][1]), slice(tile[1][0], tile[1][1]))

        # read the lon and lat tile
        lon_data = longitude[idx]
        lat_data = latitude[idx]

        # may not be processing full row wise (all columns)
        dims = lon_data.shape
        col_offset = idx[1].start

        view = np.full(dims, no_data, dtype=out_dtype)
        azi = np.full(dims, no_data, dtype=out_dtype)
        asol = np.full(dims, no_data, dtype=out_dtype)
        soazi = np.full(dims, no_data, dtype=out_dtype)
        rela_angle = np.full(dims, no_data, dtype=out_dtype)
        time = np.full(dims, no_data, dtype=out_dtype)

        # loop each row within each tile (which itself could be a single row)
        for i in range(lon_data.shape[0]):
            row_id = idx[0].start + i + 1 # FORTRAN 1 based index
            stat = angle(dims[1], acq.lines, row_id, col_offset, lat_data[i],
                         lon_data[i], spheroid[0], orbital_elements[0],
                         acq.decimal_hour(), century, 12, smodel[0], track[0],
                         view[i], azi[i], asol[i], soazi[i], rela_angle[i],
                         time[i], x_cent, n_cent)
                         # x_cent[idx[0]], n_cent[idx[0]])

            if stat != 0:
                msg = ("Error in calculating angles at row: {}.\n"
                       "No interval found in track!")
                raise RuntimeError(msg.format(row_id - 1))

        # output to disk
        sat_v_ds[idx] = view
        sat_az_ds[idx] = azi
        sol_z_ds[idx] = asol
        sol_az_ds[idx] = soazi
        rel_az_ds[idx] = rela_angle
        time_ds[idx] = time

    # outputs
    # TODO: rework create_boxline so that it reads tiled data effectively
    create_centreline_dataset(geobox, x_cent, n_cent, grp)
    create_boxline(geobox, sat_v_ds[:], grp[DatasetName.centreline.value], grp,
                   acq.maximum_view_angle)

    if out_group is None:
        return fid
