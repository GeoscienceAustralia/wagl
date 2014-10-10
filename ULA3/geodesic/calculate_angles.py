#!/usr/bin/env python
import numpy
from osgeo import gdal
from osgeo import osr
import ephem

from EOtools.DatasetDrivers import SceneDataset

from ULA3.set_satmod import set_satmod
from ULA3.set_times import set_times
from ULA3.angle_all import angle

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
        A floating point NumPy array of 4 elements containing
        spheroidal paramaters.
        Index 0 contains the spheroid Major Axis.
        Index 1 contains the spheroid Inverse Flattening.
        Index 2 contains the spheroid Squared Eccentricity.
        Index 3 contains the Earth rotational angular velocity in
        radians/second.
    """

    # Initialise the spheroid array
    spheroid = numpy.zeros((4))

    # Define the spatial reference
    sr = osr.SpatialReference()
    sr.ImportFromWkt(proj_wkt)

    # Spheroid major axis
    spheroid[0] = sr.GetSemiMajor()

    # Inverse flattening
    spheroid[1] = sr.GetInvFlattening()

    # Eccentricity squared
    spheroid[2] = 1.0 - (1.0 - 1.0/spheroid[1])**2

    # Earth rotational angular velocity rad/sec
    # Other sources such as:
    # http://www.oosa.unvienna.org/pdf/icg/2012/template/WGS_84.pdf
    # state 0.000072921150 as the mean value
    spheroid[3] = 0.000072722052

    return spheroid 

def setup_orbital_elements(ephemeral, datetime):
    """
    Given an ephemeral object and a datetime object, calculate the
    satellite orbital paramaters used for calculating the angle grids.

    :param ephemeral:
        A pyephem object already instantiated by loading a TLE (Two
        Line Element).

    :param datetime:
        A datetime object containing the date to be used in computing
        the satellite orbital paramaters.

    :return:
        A floating point NumPy array of 3 elements containing the
        satellite ephemeral bodies orbital paramaters.
        Index 0 contains the obrital inclination in degrees.
        Index 1 contains the semi major raidus in metres.
        Index 2 contains the angular velocity in radians/sec^1.
    """

    ephemeral.compute(datetime)
    pi = numpy.pi
    n = ephemeral._n # number or orbits per day
    s = 24*60*60 # Seconds in a day
    mu = 398600441800000.0 # Earth Gravitational parameter m^3s^-2

    orbital_elements = numpy.zeros((3))

    # orbital inclination (degrees)
    orbital_elements[0] = numpy.rad2deg(ephemeral._inc)

    # semi_major radius (m)
    # http://smallsats.org/2012/12/06/two-line-element-set-tle/
    orbital_elements[1] = (mu/(2*pi*n/s)**2)**(1./3)
    
    # angular velocity (rad sec-1)
    orbital_elements[2] = (2*pi*n)/s


    # For testing we'll use static values
    orbital_elements[0] = 98.200000
    orbital_elements[1] = 7083160.000000
    orbital_elements[2] = 0.001059

    return orbital_elements

def setup_smodel(centre_lon, centre_lat, spheroid, orbital_elements):
    """
    Setup the satellite model.

    :param centre_lon:
        The longitude of the scene centre.

    :param centre_lat:
        The lattitude of the scene centre.

    :param spheroid:
        A 4 element floating point array containing the Earth
        spheroidal paramaters.
        Index 0 contains the spheroid Major Axis.
        Index 1 contains the spheroid Inverse Flattening.
        Index 2 contains the spheroid Squared Eccentricity.
        Index 3 contains the Earth rotational angular velocity in
        radians/second.

    :param orbital_elements:
        A 3 element floating point array containing the satellite
        orbital elements.
        Index 0 contains the obrital inclination in degrees.
        Index 1 contains the semi major raidus in metres.
        Index 2 contains the angular velocity in radians/sec^1.

    :return:
        A floating point NumPy array of 12 elements containing the
        satellite model paramaters.
        Index 0 contains phi0.
        Index 1 contains phi0_p.
        Index 2 contains rho0.
        Index 3 contains t0.
        Index 4 contains lam0.
        Index 5 contains gamm0.
        Index 6 contains beta0.
        Index 7 contains rotn0.
        Index 8 contains hxy0.
        Index 9 contains N0.
        Index 10 contains H0.
        Index 11 contains th_ratio0.
    """

    smodel, istat = set_satmod(centre_lon, centre_lat, spheroid, orbital_elements)

    return smodel

def setup_times(ymin, ymax, spheroid, orbital_elements, smodel, npoints=12):
    """
    Setup the satellite track times.

    :param ymin:
        The minimum lattitude in the array extent.

    :param ymax:
        The maximum lattitude in the array extent.

    :param spheroid:
        A 4 element floating point array containing the Earth
        spheroidal paramaters.
        Index 0 contains the spheroid Major Axis.
        Index 1 contains the spheroid Inverse Flattening.
        Index 2 contains the spheroid Squared Eccentricity.
        Index 3 contains the Earth rotational angular velocity in
        radians/second.

    :param orbital_elements:
        A 3 element floating point array containing the satellite
        orbital elements.
        Index 0 contains the obrital inclination in degrees.
        Index 1 contains the semi major raidus in metres.
        Index 2 contains the angular velocity in radians/sec^1.

    :param smodel:
        A floating point NumPy array of 12 elements containing the
        satellite model paramaters.
        Index 0 contains phi0.
        Index 1 contains phi0_p.
        Index 2 contains rho0.
        Index 3 contains t0.
        Index 4 contains lam0.
        Index 5 contains gamm0.
        Index 6 contains beta0.
        Index 7 contains rotn0.
        Index 8 contains hxy0.
        Index 9 contains N0.
        Index 10 contains H0.
        Index 11 contains th_ratio0.

    :param npoints:
        The number of time sample points to be calculated along the
        satellite track. Default is 12

    :return:
        A floating point NumPy array of [npoints,8] containing the
        satellite track times and other information.
        Index 0 t
        Index 1 rho
        Index 2 phi_p
        Index 3 lam
        Index 4 beta
        Index 5 hxy
        Index 6 mj
        Index 7 skew
    """

    track, istat = set_times(ymin, ymax, npoints, spheroid, orbital_elements, smodel)

    return track

def calculate_angles(scene_dataset, lon_array, lat_array, npoints=12):
    """
    Calcualte the satellite view, satellite azimuth, solar zenith,
    solar azimuth, and relative aziumth angle grids, as well as the
    time grid. All grids are output as float32 ENVI files.
    The CENTRELINE is also calculated and output to disk.

    :param scene_dataset:
        An object based on the SceneDataset Class.

    :param lon_array:
        A 2D float64 NumPy array containing the longitude values.

    :param lat_array:
        A 2D float64 NumPy array containing the lattitude values.

    :param npoints:
        The number of time sample points to be calculated along the
        satellite track. Default is 12

    :return:
        6 float32 NumPy arrays of the same shape as lon_array:
        Satellite zenith angle.
        Satellite azimuth angle.
        Solar zenith angle.
        Solar azimuth angle.
        Relative azimuth angle.
        Time.
        3 float32 NumPy arrays with the same shape as
        lon_array.shape[0] containing the array coodinates of the
        satellite track line and N_Cent(??):
        Y_coordinates (Starting at 1)
        X_coordinates
        N_cent (Value 2 if centre x coordinate was averaged and
                1 if centre x coordinate was not averaged)
    """

    # Image projection, geotransform and scene centre time stamp
    prj = scene_dataset.GetProjection()
    geoT = scene_dataset.GetGeoTransform()
    centre_datetime = scene_dataset.scene_centre_datetime

    # Get the lat/lon of the scene centre
    centre_xy = numpy.array(scene_dataset.lonlats['CENTRE'])

    # Get the earth spheroidal paramaters
    spheroid = setup_spheroid(scene_dataset.GetProjection())

    # Get the satellite orbital elements
    tle_dir = '/g/data1/v10/eoancillarydata/sensor-specific'
    sat_ephemeral = scene_dataset.satellite.load_tle(centre_datetime, tle_dir)
    orbital_elements = setup_orbital_elements(sat_ephemeral, centre_datetime)

    # Min and Max lat extents
    lat_min_max = scene_dataset.get_bounds()[1]

    # Scene centre in time stamp in decimal hours
    hours = scene_dataset.decimal_hour

    # Calculate the julian century past JD2000
    century = calculate_julian_century(centre_datetime)

    # Need something to determine max satellite view angle
    # Currently not even used in Fuqin's code
    view_max = 9.0

    # Get the satellite model paramaters
    smodel = setup_smodel(centre_xy[0], centre_xy[1], spheroid, orbital_elements)

    # Get the times and satellite track information
    track = setup_times(lat_min_max[0], lat_min_max[1], spheroid, orbital_elements, smodel, npoints)

    # Array dimensions
    dims = lon_array.shape
    cols = dims[1]
    rows = dims[0]

    # Initialise the arrays to hold the angles
    view = numpy.zeros(dims, dtype='float32')
    azi = numpy.zeros(dims, dtype='float32')
    asol = numpy.zeros(dims, dtype='float32')
    soazi = numpy.zeros(dims, dtype='float32')
    rela_angle = numpy.zeros(dims, dtype='float32')
    time = numpy.zeros(dims, dtype='float32')

    # Set to null value
    view[:] = -999
    azi[:] = -999
    asol[:] = -999
    soazi[:] = -999
    rela_angle[:] = -999
    time[:] = -999

    # Initialise centre line variables
    Y_cent = numpy.arange(1,rows+1).astype('float32')
    X_cent = numpy.zeros((rows), dtype='float32')
    N_cent = numpy.zeros((rows), dtype='float32')

    # Loop over each row
    for i in range(rows):
        istat = angle(cols, rows, i+1, lat_array[i], lon_array[i], spheroid,
                      orbital_elements, hours, century, npoints, smodel, track,
                      view[i], azi[i], asol[i], soazi[i], rela_angle[i],
                      time[i], X_cent, N_cent)

    # Centreline
    # here need code to write the track in the image as an ascii file
    # if more than one pixel in a line was a track point the coordinates
    # are averaged
    wh = N_cent > 1.5
    X_cent[wh] = X_cent[wh] / N_cent[wh]

    # check whether there is no centre pixel in the line. It is assumed that
    # at least the adjacent lines have pixel
    wh = N_cent < 0.5
    temp = X_cent[0:2].copy()
    X_cent[wh] = numpy.roll(X_cent, 1)[wh]
    # Account for first element potentially being changed with the last element
    if wh[0]:
        X_cent[0] = temp[1]

    # Convert X & Y centre points to integers (basically array co-ordinates)
    Y_cent = numpy.rint(Y_cent)
    X_cent = numpy.rint(X_cent)

    return (view, azi, asol, soazi, rela_angle, time, Y_cent, X_cent, N_cent)
