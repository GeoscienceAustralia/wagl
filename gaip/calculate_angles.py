#!/usr/bin/env python

"""
Satellite and Solar angle claculations over a 2D grid.
"""

import math
import ephem
import numpy as np
import h5py
import pandas

from osgeo import osr
from eotools import tiling
from gaip import gridded_geo_box
from gaip import load_tle
from gaip import angle
from gaip import set_satmod
from gaip import set_times
from gaip import cstart_cend
from gaip import dataset_compression_kwargs
from gaip import attach_image_attributes
from gaip import attach_table_attributes

CRS = "EPSG:4326"


def _calculate_angles_wrapper(acquisition, lon_fname, lon_dname, lat_fname,
                              lat_dname, out_fname, npoints=12,
                              compression='lzf', max_angle=9.0, tle_path=None):
    """
    A private wrapper for dealing with the internal custom workings of the
    NBAR workflow.
    """
    if lon_fname == lat_fname:
        with h5py.File(lon_fname, 'r') as src:
            lon_ds = src[lon_dname]
            lat_ds = src[lat_dname]
            fid = calculate_angles(acquisition, lon_ds, lat_ds, npoints,
                                   out_fname, compression, max_angle=max_angle,
                                   tle_path=tle_path)
    else:
        with h5py.File(lon_fname, 'r') as lon_src,\
            h5py.File(lat_fname, 'r') as lat_src:

            lon_ds = lon_src[lon_dname]
            lat_ds = lat_src[lat_dname]
            fid = calculate_angles(acquisition, lon_ds, lat_ds, npoints,
                                   out_fname, compression, max_angle=max_angle,
                                   tle_path=tle_path)

    fid.close()

    return


def create_centreline_dataset(geobox, y, x, n):
    """
    Creates the centre line dataset.

    :param geobox:
        An instance of a GriddedGeoBox object.

    :param y:
        A 1D np array of type int with the same shape as x & n.
        Details the row number starting at 1.

    :param x:
        A 1D np array of type int with the same shape as y & n.
        Details the column number starting at 0.

    :param n:
        A 1D np array of type int with the same shape as y & x.
        Details whether or not the track point coordinate is
        averaged.

    :return:
        A `NumPy` dataset with the following datatype:

        * [('row_index', 'int64'), ('col_index', 'int64'),
           ('n_pixels', 'float'), ('latitude', 'float'),
           ('longitude', 'float')]
    """
    rows = y.shape[0]

    # Define the TO_CRS for lon & lat outputs
    sr = osr.SpatialReference()
    sr.SetFromUserInput(CRS)

    dtype = np.dtype([('row_index', 'int64'), ('col_index', 'int64'),
                      ('n_pixels', 'float'), ('latitude', 'float'),
                      ('longitude', 'float')])
    data = np.zeros(rows, dtype=dtype)

    for r in range(rows):
        # We offset by -1 to get the zero based col and row id
        mapXY = geobox.convert_coordinates((x[r] - 1, y[r] - 1))
        lon, lat = geobox.transform_coordinates(mapXY, to_crs=sr)
        data['row_index'][r] = y[r]
        data['col_index'][r] = x[r]
        data['n_pixels'][r] = n[r]
        data['latitude'][r] = lat
        data['longitude'][r] = lon

    return data


def create_boxline_coordinator(view_angle_dataset, line, ncentre, npoints,
                               max_angle=9.0):
    """
    Creates the boxline and coordinator datasets.

    :param view_angle_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the satellite view
        angles when index/sliced.

    :param line:
        A 1D NumPy array containing the values 1->n for n lines.
        1 based index.

    :param ncentre:
        A 1D NumPy array containing the values for the column
        coordinate for the central index; 1 based index.

    :param npoints:
        A 1D NumPy array containing the values for the number of points
        used to determine the pixel index.

    :param max_angle:
        The maximum viewing angle. Default is 9.0 degrees.

    :return:
        2 `NumPy` datasets:

    boxline_dtype = np.dtype([('row_index', 'int64'),
                                 ('bisection_index', 'int64'),
                                 ('start_index', 'int64'),
                                 ('end_index', 'int64')])
    coordinator_dtype = np.dtype([('row_index', 'int64'),
                                  ('col_index', 'int64')])
        * boxline; dtype = [('row_index', 'int64'),
                            ('bisection_index', 'int64'),
                            ('start_index', 'int64'),
                            ('end_index', 'int64')]
        * coordinator; dtype = [('row_index', 'int64'),
                                ('bisection_index', 'int64'),
                                ('start_index', 'int64'),
                                ('end_index', 'int64')]
    """
    rows, cols = view_angle_dataset.shape

    # allocate the output arrays
    istart = np.zeros(rows, dtype='int')
    iend = np.zeros(rows, dtype='int')

    # calculate the column start and end indices
    cstart_cend(cols, rows, max_angle, view_angle_dataset[:].transpose(),
                line, ncentre, istart, iend)

    # TODO: Fuqin to document, and these results are only used for
    # granules that do not contain the satellite track path
    kk = ll = -1
    for i in range(rows):
        if npoints[i] > 0.5:
            if kk == -1:
                kk = i
            ll = i

    coordinator = np.zeros((9, 2), dtype='int64')
    mid_col = cols // 2
    mid_row = rows // 2
    mid_row_idx = mid_row - 1

    # Sentinel-2a doesn't require a start end index, but Landsat does
    # the following should satisfy both use cases, that way we don't
    # require two separate bilinear functions, and two separate boxline
    # (or bisection co-ordinator functions)
    # only used for the case "normal center line"
    start_col = np.zeros(ncentre.shape, dtype='int')
    end_col = np.zeros(ncentre.shape, dtype='int')
    start_col.fill(1)
    end_col.fill(cols)

    np.maximum(istart, start_col, out=start_col)
    np.minimum(iend, end_col, out=end_col)

    df = pandas.DataFrame({'line': line,
                           'bisection': ncentre,
                           'npoints': npoints})

    # do we have a hit with the satellite track falling within the granule
    if npoints[0] > 0.5:
        if npoints[-1] >= 0.5:
            # normal center line
            coordinator[0] = [line[0], start_col[0]]
            coordinator[1] = [line[0], ncentre[0]]
            coordinator[2] = [line[0], end_col[0]]
            coordinator[3] = [line[mid_row_idx], start_col[mid_row_idx]]
            coordinator[4] = [line[mid_row_idx], ncentre[mid_row_idx]]
            coordinator[5] = [line[mid_row_idx], end_col[mid_row_idx]]
            coordinator[6] = [line[-1], start_col[-1]]
            coordinator[7] = [line[-1], ncentre[-1]]
            coordinator[8] = [line[-1], end_col[-1]]

            data = df[['line', 'bisection']].copy()
            data['start'] = start_col
            data['end'] = end_col
        elif npoints[-1] < 0.5:
            # first half centerline
            coordinator[0] = [line[0], 1]
            coordinator[1] = [line[0], ncentre[0]]
            coordinator[2] = [line[0], cols]
            coordinator[3] = [line[ll], 1]
            coordinator[4] = [line[ll], ncentre[0]]
            coordinator[5] = [line[ll], cols]
            coordinator[6] = [line[-1], cols]
            coordinator[7] = [line[-1], cols]
            coordinator[8] = [line[-1], cols]

            data = df[['line']].copy()
            data['bisection'] = ncentre[0]
            data['start'] = start_col
            data['end'] = end_col
        elif npoints[-1] >= 0.5:
            # last half of center line
            coordinator[0] = [line[0], 1]
            coordinator[1] = [line[0], ncentre[ll]]
            coordinator[2] = [line[0], cols]
            coordinator[3] = [line[kk], 1]
            coordinator[4] = [line[kk], ncentre[ll]]
            coordinator[5] = [line[kk], cols]
            coordinator[6] = [line[-1], 1]
            coordinator[7] = [line[-1], ncentre[ll]]
            coordinator[8] = [line[-1], cols]

            data = df[['line']].copy()
            data['bisection'] = ncentre[ll]
            data['start'] = start_col
            data['end'] = end_col
    else:
        # no centre line
        coordinator[0] = [line[0], 1]
        coordinator[1] = [line[0], mid_col]
        coordinator[2] = [line[0], cols]
        coordinator[3] = [line[mid_row_idx], 1]
        coordinator[4] = [line[mid_row_idx], mid_col]
        coordinator[5] = [line[mid_row_idx], cols]
        coordinator[6] = [line[-1], 1]
        coordinator[7] = [line[-1], mid_col]
        coordinator[8] = [line[-1], cols]

        data = df[['line']].copy()
        data['bisection'] = mid_col
        data['start'] = start_col
        data['end'] = end_col

    columns = ['row_index', 'bisection_index', 'start_index', 'end_index']
    boxline_dtype = np.dtype([(col, 'int64') for col in columns])
    coordinator_dtype = np.dtype([('row_index', 'int64'),
                                  ('col_index', 'int64')])

    boxline_dset = np.zeros(rows, dtype=boxline_dtype)
    boxline_dset['row_index'] = data['line'].values
    boxline_dset['bisection_index'] = data['bisection'].values
    boxline_dset['start_index'] = data['start'].values
    boxline_dset['end_index'] = data['end'].values

    coordinator_dset = np.zeros(9, dtype=coordinator_dtype)
    coordinator_dset[0] = coordinator[0]
    coordinator_dset[1] = coordinator[1]
    coordinator_dset[2] = coordinator[2]
    coordinator_dset[3] = coordinator[3]
    coordinator_dset[4] = coordinator[4]
    coordinator_dset[5] = coordinator[5]
    coordinator_dset[6] = coordinator[6]
    coordinator_dset[7] = coordinator[7]
    coordinator_dset[8] = coordinator[8]

    return boxline_dset, coordinator_dset


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
    epoch = ephem.date((2000, 01, 01, 12.00))
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

    # Initialise the spheroid array
    spheroid = np.zeros((4))

    # Define the spatial reference
    sr = osr.SpatialReference()
    sr.ImportFromWkt(proj_wkt)

    # Spheroid major axis
    dset['semi_major_axis'] = sr.GetSemiMajor()

    # Inverse flattening
    dset['inverse_flattening'] = sr.GetInvFlattening()

    # Eccentricity squared
    dset['eccentricity_squared'] = 1.0 - (1.0 - 1.0 / spheroid[1]) ** 2

    # Earth rotational angular velocity rad/sec
    # Other sources such as:
    # http://www.oosa.unvienna.org/pdf/icg/2012/template/WGS_84.pdf
    # state 0.000072921150 as the mean value
    dset['earth_rotational_angular_velocity'] = 0.000072722052

    return np.array(dset.tolist()).squeeze(), dset


def setup_orbital_elements(ephemeral, datetime, acquisition):
    """
    Given an ephemeral object and a datetime object, calculate the
    satellite orbital paramaters used for calculating the angle grids.

    :param ephemeral:
        A pyephem object already instantiated by loading a TLE (Two
        Line Element).

    :param datetime:
        A datetime object containing the date to be used in computing
        the satellite orbital paramaters.

    :param acquisition:
        An `Acquisition` object.

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
        ephemeral.compute(datetime)
        pi = np.pi
        n = ephemeral._n  # number or orbits per day
        s = 24 * 60 * 60  # Seconds in a day
        mu = 398600441800000.0  # Earth Gravitational parameter m^3s^-2

        # orbital inclination (degrees)
        dset['orbital_inclination'] = np.rad2deg(ephemeral._inc)

        # semi_major radius (m)
        # http://smallsats.org/2012/12/06/two-line-element-set-tle/
        dset['semi_major_radius'] = (mu / (2 * pi * n / s) ** 2) ** (1. / 3)

        # angular velocity (rad sec-1)
        dset['angular_velocity'] = (2 * pi * n) / s

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


def setup_times(ymin, ymax, spheroid, orbital_elements, smodel, npoints=12):
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

    :param npoints:
        The number of time sample points to be calculated along the
        satellite track. Default is 12.

    :return:
        A floating point np array of [npoints,8] containing the
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
    track, _ = set_times(ymin, ymax, npoints, spheroid, orbital_elements,
                         smodel)

    columns = ['t', 'rho', 'phi_p', 'lam', 'beta', 'hxy', 'mj', 'skew']
    dtype = np.dtype([(col, 'float64') for col in columns])
    track_dset = np.zeros(npoints, dtype=dtype)
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
    group = fid.create_group('parameters')

    for key in params:
        group.attrs[key] = params[key]

    # sheroid
    desc = "The spheroid used in the satelite and solar angles calculation."
    attrs = {'Description': desc}
    sph_dset = group.create_dataset('spheroid', data=spheriod)
    attach_table_attributes(sph_dset, title='Spheroid', attrs=attrs)

    # orbital elements
    desc = ("The satellite orbital parameters used in the satellite and "
            "solar angles calculation.")
    attrs = {'Description': desc}
    orb_dset = group.create_dataset('orbital-elements', data=orbital_elements)
    attach_table_attributes(orb_dset, title='Orbital Elements', attrs=attrs)

    # satellite model
    desc = ("The satellite model used in the satelite and solar angles "
            "calculation.")
    attrs = {'Description': desc}
    sat_dset = group.create_dataset('satellite-model', data=satellite_model)
    attach_table_attributes(sat_dset, title='Satellite Model', attrs=attrs)

    # satellite track
    desc = ("The satellite track information used in the satelite and solar "
            "angles calculation.")
    attrs = {'Description': desc}
    track_dset = group.create_dataset('satellite-track', data=satellite_track)
    attach_table_attributes(track_dset, title='Satellite Track', attrs=attrs)

    # flush
    fid.flush()


def calculate_angles(acquisition, lon_dataset, lat_dataset, npoints=12,
                     out_fname=None, compression='lzf', max_angle=9.0,
                     tle_path=None):
    """
    Calculate the satellite view, satellite azimuth, solar zenith,
    solar azimuth, and relative aziumth angle grids, as well as the
    time grid. All grids are output as float32 ENVI files.
    A wrapper routine for the ``angle_all`` Fortran module built via
    ``F2Py``.

    :param acquisition:
        An instance of an `Acquisition` object.

    :param lon_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the longitude
        values when index/sliced.
        The dimensions must match that of the `acquisition` objects's
        samples (x) and lines (y) parameters.

    :param lat_dataset:
        A `NumPy` or `NumPy` like dataset that allows indexing
        and returns a `NumPy` dataset containing the latitude
        values when index/sliced.
        The dimensions must match that of the `acquisition` objects's
        samples (x) and lines (y) parameters.

    :param npoints:
        The number of time sample points to be calculated along the
        satellite track. Default is 12

    :param out_fname:
        If set to None (default) then the results will be returned
        as an in-memory hdf5 file, i.e. the `core` driver.
        Otherwise it should be a string containing the full file path
        name to a writeable location on disk in which to save the HDF5
        file.

        The dataset names will be as follows:

        * satellite-view
        * satellite-azimuth
        * solar-zenith
        * solar-azimuth
        * relative-azimuth
        * acquisition-time

    :param compression:
        The compression filter to use. Default is 'lzf'.
        Options include:

        * 'lzf' (Default)
        * 'lz4'
        * 'mafisc'
        * An integer [1-9] (Deflate/gzip)

    :param max_angle:
        The maximum satellite view angle to use within the workflow.
        Default is 9.0 degrees.

    :return:
        An opened `h5py.File` object, that is either in-memory using the
        `core` driver, or on disk.
    """
    # Get the datetime of the acquisition
    dt = acquisition.scene_center_datetime

    # Compute the geobox
    geobox = gridded_geo_box(acquisition)

    # Image projection
    prj = geobox.crs.ExportToWkt()

    # Min and Max lat extents
    # This method should handle northern and southern hemispheres
    min_lat = min(min(geobox.ul_lonlat[1], geobox.ur_lonlat[1]),
                  min(geobox.ll_lonlat[1], geobox.lr_lonlat[1]))
    max_lat = max(max(geobox.ul_lonlat[1], geobox.ur_lonlat[1]),
                  max(geobox.ll_lonlat[1], geobox.lr_lonlat[1]))

    # temporary lat/lon buffer for satellite track calculations
    min_lat -= 1
    max_lat += 1

    # Get the lat/lon of the scene centre
    # check if we have a file with GPS satellite track points
    # which can be used for cases of image granules/tiles, eg Sentinel-2A
    if acquisition.gps_file:
        points = acquisition.read_gps_file()
        subs = points[(points.lat >= min_lat) & (points.lat <= max_lat)]
        idx = subs.shape[0] // 2 - 1
        centre_xy = (subs.iloc[idx].lon, subs.iloc[idx].lat)
    else:
        centre_xy = geobox.centre_lonlat

    # Get the earth spheroidal paramaters
    spheroid, spheroid_dset = setup_spheroid(prj)

    # Get the satellite orbital elements
    sat_ephemeral = load_tle(acquisition, tle_path)

    orbital_elements, orb_dset = setup_orbital_elements(sat_ephemeral, dt,
                                                        acquisition)

    # Scene centre in time stamp in decimal hours
    hours = acquisition.decimal_hour

    # Calculate the julian century past JD2000
    century = calculate_julian_century(dt)

    # Get the satellite model paramaters
    smodel, smodel_dset = setup_smodel(centre_xy[0], centre_xy[1], spheroid,
                                       orbital_elements)

    # Get the times and satellite track information
    track, track_dset = setup_times(min_lat, max_lat, spheroid,
                                    orbital_elements, smodel, npoints)

    # Array dimensions
    cols = acquisition.samples
    rows = acquisition.lines
    dims = (rows, cols)

    # Initialise 1D arrays to hold the angles
    out_dtype = 'float32'
    view = np.zeros((1, cols), dtype=out_dtype)
    azi = np.zeros((1, cols), dtype=out_dtype)
    asol = np.zeros((1, cols), dtype=out_dtype)
    soazi = np.zeros((1, cols), dtype=out_dtype)
    rela_angle = np.zeros((1, cols), dtype=out_dtype)
    time = np.zeros((1, cols), dtype=out_dtype)

    # Initialise the output files
    if out_fname is None:
        fid = h5py.File('satellite-solar-angles.h5', driver='core',
                        backing_store=False)
    else:
        fid = h5py.File(out_fname, 'w')

    # store the parameter settings used with the satellite and solar angles
    # function
    params = {'dimensions': dims,
              'lines': rows,
              'samples': cols,
              'century': century,
              'hours': hours,
              'scene_acquisition_datetime_iso': dt.isoformat(),
              'centre_longitude_latitude': centre_xy,
              'minimum_latiude': min_lat,
              'maximum_latiude': max_lat,
              'latitude_buffer': '1.0 degrees',
              'max_satellite_viewing_angle': max_angle}
    _store_parameter_settings(fid, spheroid_dset, orb_dset,
                              smodel_dset, track_dset, params)

    no_data = -999
    kwargs = dataset_compression_kwargs(compression=compression,
                                        chunks=(1, cols))
    kwargs['shape'] = dims
    kwargs['fillvalue'] = no_data
    kwargs['dtype'] = out_dtype

    sat_v_ds = fid.create_dataset('satellite-view', **kwargs)
    sat_az_ds = fid.create_dataset('satellite-azimuth', **kwargs)
    sol_z_ds = fid.create_dataset('solar-zenith', **kwargs)
    sol_az_ds = fid.create_dataset('solar-azimuth', **kwargs)
    rel_az_ds = fid.create_dataset('relative-azimuth', **kwargs)
    time_ds = fid.create_dataset('acquisition-time', **kwargs)

    # attach some attributes to the image datasets
    attrs = {'crs_wkt': geobox.crs.ExportToWkt(),
             'geotransform': geobox.affine.to_gdal(),
             'no_data_value': no_data}
    desc = "The satellite viewing angle in degrees."
    attrs['Description'] = desc
    attach_image_attributes(sat_v_ds, attrs)

    desc = "The satellite azimuth angle in degrees."
    attrs['Description'] = desc
    attach_image_attributes(sat_az_ds, attrs)

    desc = "The solar zenith angle in degrees."
    attrs['Description'] = desc
    attach_image_attributes(sol_z_ds, attrs)

    desc = "The solar azimuth angle in degrees."
    attrs['Description'] = desc
    attach_image_attributes(sol_az_ds, attrs)

    desc = "The relative azimuth angle in degrees."
    attrs['Description'] = desc
    attach_image_attributes(rel_az_ds, attrs)

    desc = ("The satellite acquisition time grid in seconds before and after "
            "the scene acquisition datetime.")
    attrs['Description'] = desc
    attach_image_attributes(time_ds, attrs)

    # Initialise centre line variables
    y_cent = np.arange(1, rows + 1).astype('float32')
    x_cent = np.zeros((rows), dtype='float32')
    n_cent = np.zeros((rows), dtype='float32')

    # Initialise the tile generator for processing
    # Process 1 row of data at a time
    tiles = tiling.generate_tiles(cols, rows, cols, 1)

    for i, tile in enumerate(tiles):
        idx = (slice(tile[0][0], tile[0][1]), slice(tile[1][0], tile[1][1]))

        # read the lon and lat tile
        lon_data = lon_dataset[idx]
        lat_data = lat_dataset[idx]

        # set to null value
        view[:] = no_data
        azi[:] = no_data
        asol[:] = no_data
        soazi[:] = no_data
        rela_angle[:] = no_data
        time[:] = no_data

        stat = angle(cols, rows, i + 1, lat_data, lon_data, spheroid,
                     orbital_elements, hours, century, npoints, smodel,
                     track, view[0], azi[0], asol[0], soazi[0], rela_angle[0],
                     time[0], x_cent, n_cent)

        if stat != 0:
            msg = ("Error in calculating angles at row: {}.\n"
                   "No interval found in track!")
            raise RuntimeError(msg.format(i))

        # output to disk
        sat_v_ds[idx] = view
        sat_az_ds[idx] = azi
        sol_z_ds[idx] = asol
        sol_az_ds[idx] = soazi
        rel_az_ds[idx] = rela_angle
        time_ds[idx] = time

    # centreline
    # here need code to write the track in the image as an ascii file
    # if more than one pixel in a line was a track point the coordinates
    # are averaged
    wh = n_cent > 1.5
    x_cent[wh] = x_cent[wh] / n_cent[wh]

    # check whether there is no centre pixel in the line. It is assumed that
    # at least the adjacent lines have pixel
    wh = n_cent < 0.5
    temp = x_cent[0:2].copy()
    x_cent[wh] = np.roll(x_cent, 1)[wh]
    # account for first element potentially being changed with the
    # last element
    if wh[0]:
        x_cent[0] = temp[1]

    # convert X & Y centre points to integers (basically array co-ordinates)
    y_cent = np.rint(y_cent)
    x_cent = np.rint(x_cent)

    # create the dataset and save to the HDF5 file
    centreline_dataset = create_centreline_dataset(geobox, y_cent, x_cent,
                                                   n_cent)
    kwargs = dataset_compression_kwargs(compression=compression)
    cent_dset = fid.create_dataset('centreline', data=centreline_dataset,
                                   **kwargs)
    desc = ("Contains the array, latitude and longitude coordinates of the "
            "satellite track path.")
    attrs = {'Description': desc,
             'array_coordinate_offset': -1}
    attach_table_attributes(cent_dset, title='Centreline', attrs=attrs)

    # boxline and coordinator
    boxline, coordinator = create_boxline_coordinator(sat_v_ds, y_cent, x_cent,
                                                      n_cent,
                                                      max_angle=max_angle)
    desc = ("Contains the bi-section, column start and column end array "
            "coordinates.")
    attrs['Description'] = desc
    box_dset = fid.create_dataset('boxline', data=boxline, **kwargs)
    attach_table_attributes(box_dset, title='Boxline', attrs=attrs)

    desc = ("Contains the row and column array coordinates used for the "
            "atmospheric calculations.")
    attrs['Description'] = desc
    coord_dset = fid.create_dataset('coordinator', data=coordinator, **kwargs)
    attach_table_attributes(coord_dset, title='Coordinator', attrs=attrs)

    fid.flush()

    return fid
