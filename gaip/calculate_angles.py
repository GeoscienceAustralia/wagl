"""
Grid Calculations.
"""

import math
import rasterio
import ephem
import numpy as np
import os

from osgeo import gdal
from osgeo import osr
from EOtools import tiling
from gaip import acquisitions
from gaip import find_file
from gaip import gridded_geo_box
from gaip import load_tle
from gaip import angle
from gaip import set_satmod
from gaip import set_times

CRS = "EPSG:4326"
TLE_DIR = '/g/data1/v10/eoancillarydata/sensor-specific'

# To be used as a template while gaip is restructured


def sat_sol_grid_workflow(l1t_path, work_path, lonlat_path):
    """
    Workflow to generate the satellite and solar grids.

    :note:
        This routine is only utilised by the unittesting framework.
    """
    # Retrieve an acquisitions object
    acqs = acquisitions(l1t_path)

    # Create the geobox
    geobox = gridded_geo_box(acqs[0])

    # Find and open the longitude and lattitude files
    # Avoiding DataManger here. find_file will be used sparingly until a
    # proper workflow is written.
    lon_fname = find_file(lonlat_path, 'LON.bin')
    lat_fname = find_file(lonlat_path, 'LAT.bin')

    # Get the array dimensions from the first acquisistion
    # The dimensions should match for all bands except the panchromatic
    cols = acqs[0].samples

    # Define the output file names
    sat_view_zenith_fname = os.path.join(work_path, 'SATELLITE_VIEW.bin')
    sat_azimuth_fname = os.path.join(work_path, 'SATELLITE_AZIMUTH.bin')
    solar_zenith_fname = os.path.join(work_path, 'SOLAR_ZENITH.bin')
    solar_azimuth_fname = os.path.join(work_path, 'SOLAR_AZIMUTH.bin')
    relative_azimuth_fname = os.path.join(work_path, 'RELATIVE_AZIMUTH.bin')
    time_fname = os.path.join(work_path, 'TIME.bin')

    out_fnames = [sat_view_zenith_fname, sat_azimuth_fname,
                  solar_zenith_fname, solar_azimuth_fname,
                  relative_azimuth_fname, time_fname]

    # Get the angles, time, & satellite track coordinates
    (satellite_zenith, satellite_azimuth, solar_zenith,
     solar_azimuth, relative_azimuth, time,
     y_cent, x_cent, n_cent) = calculate_angles(acqs[0], lon_fname,
                                                lat_fname, npoints=12,
                                                to_disk=out_fnames)

    # Write out the CENTRELINE file
    create_centreline_file(geobox, y_cent, x_cent, n_cent, cols,
                           view_max=9.0)


def create_centreline_file(geobox, y, x, n, cols, view_max,
                           outfname='CENTRELINE'):
    """
    Creates the centre line text file.

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

    :param cols:
        An integer indicating the number of columns in the original
        array.

    :param view_max:
        An float indicating the maximum view angle of the satellite
        used for determining the centreline.

    :param outfname:
        A string containing the full file system pathname of the
        output file.
        Default is CENTRELINE produced in the current working
        directory.
    """

    rows = y.shape[0]

    # Define the TO_CRS for lon & lat outputs
    sr = osr.SpatialReference()
    sr.SetFromUserInput(CRS)

    with open(outfname, 'w') as outfile:

        # Right justified with length of 14 per item
        outfile.write('{view_max:>14}\n'.format(view_max=view_max))
        outfile.write('{rows:>14}{cols:>14}\n'.format(rows=rows, cols=cols))

        for r in range(rows):
            # We offset by -1 to get the zero based col and row id
            mapXY = geobox.convert_coordinates((x[r] - 1, y[r] - 1))
            lon, lat = geobox.transform_coordinates(mapXY, to_crs=sr)
            # Right justified at various lengths
            msg = '{row:>14}{col:>14}{n:>14}{lat:>21}{lon:>21}\n'
            msg = msg.format(row=int(y[r]), col=int(x[r]), n=n[r], lat=lat,
                             lon=lon)
            outfile.write(msg)


def create_header_angle_file(acquisition, view_max, outfname='HEADERANGLE'):
    """
    Creates the header angle text file.

    :param acquisition:
        An instance of an acquisitions object.

    :param view_max:
        An float indicating the maximum view angle of the satellite
        used for determining the centreline.

    :param outfname:
        A string containing the full file system pathname of the
        output file.
        Default is HEADERANGLE produced in the current working
        directory..
    """
    # Get the satellite orbital elements
    sat_ephemeral = load_tle(acquisition, TLE_DIR)

    # If we have None, then no suitable TLE was found, so use values gathered
    # by the acquisition object
    if sat_ephemeral is None:
        # orbital inclination (degrees)
        orb_incl = math.degrees(acquisition.inclination)
        # semi_major radius (m)
        orb_radius = acquisition.semi_major_axis
        # angular velocity (rad sec-1)
        omega = acquisition.omega
        orbital_elements = np.array([orb_incl, orb_radius, omega],
                                    dtype='float')
    else:
        dt = acquisition.scene_center_datetime
        orbital_elements = setup_orbital_elements(sat_ephemeral, dt)
        orb_incl = orbital_elements[0]
        orb_radius = orbital_elements[1]
        omega = orbital_elements[2]

    # Scene centre datetime info
    year = acquisition.scene_center_datetime.year
    month = acquisition.scene_center_datetime.month
    day = acquisition.scene_center_datetime.day
    hours = acquisition.decimal_hour

    # Spatial info
    geobox = acquisition.gridded_geo_box()
    cx, cy = geobox.centre_lonlat
    samples, lines = geobox.getShapeXY()

    # Output to disk
    with open(outfname, 'w') as outfile:
        header_file = ("{year} {month} {day} {hours}\n"
                       "{lines} {samples}\n"
                       "{centre_lat} {centre_lon}\n"
                       "{radius} {inclination} {velocity}\n"
                       "{max_angle}")

        outfile.write(header_file.format(year=year, month=month, day=day,
                                         hours=hours, lines=lines,
                                         samples=samples, centre_lat=cy,
                                         centre_lon=cx, radius=orb_radius,
                                         inclination=orb_incl,
                                         velocity=omega, max_angle=view_max))


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
        Index 0 contains the spheroid Major Axis.
        Index 1 contains the spheroid Inverse Flattening.
        Index 2 contains the spheroid Squared Eccentricity.
        Index 3 contains the Earth rotational angular velocity in
        radians/second.
    """

    # Initialise the spheroid array
    spheroid = np.zeros((4))

    # Define the spatial reference
    sr = osr.SpatialReference()
    sr.ImportFromWkt(proj_wkt)

    # Spheroid major axis
    spheroid[0] = sr.GetSemiMajor()

    # Inverse flattening
    spheroid[1] = sr.GetInvFlattening()

    # Eccentricity squared
    spheroid[2] = 1.0 - (1.0 - 1.0 / spheroid[1]) ** 2

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
        A floating point np array of 3 elements containing the
        satellite ephemeral bodies orbital paramaters.
        Index 0 contains the obrital inclination in degrees.
        Index 1 contains the semi major raidus in metres.
        Index 2 contains the angular velocity in radians/sec^1.
    """

    ephemeral.compute(datetime)
    pi = np.pi
    n = ephemeral._n  # number or orbits per day
    s = 24 * 60 * 60  # Seconds in a day
    mu = 398600441800000.0  # Earth Gravitational parameter m^3s^-2

    orbital_elements = np.zeros((3))

    # orbital inclination (degrees)
    orbital_elements[0] = np.rad2deg(ephemeral._inc)

    # semi_major radius (m)
    # http://smallsats.org/2012/12/06/two-line-element-set-tle/
    orbital_elements[1] = (mu / (2 * pi * n / s) ** 2) ** (1. / 3)

    # angular velocity (rad sec-1)
    orbital_elements[2] = (2 * pi * n) / s

    return orbital_elements


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
        A floating point np array of 12 elements containing the
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

    smodel, istat = set_satmod(centre_lon, centre_lat, spheroid,
                               orbital_elements)

    return smodel


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
        A floating point np array of 12 elements containing the
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
        A floating point np array of [npoints,8] containing the
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

    track, istat = set_times(ymin, ymax, npoints, spheroid, orbital_elements,
                             smodel)

    return track


def calculate_angles(acquisition, lon_fname, lat_fname, npoints=12,
                     to_disk=None):
    """
    Calculate the satellite view, satellite azimuth, solar zenith,
    solar azimuth, and relative aziumth angle grids, as well as the
    time grid. All grids are output as float32 ENVI files.
    A wrapper routine for the ``angle_all`` Fortran module built via
    ``F2Py``.

    :param acquisition:
        An instance of an acquisitions object.

    :param lon_fname:
        A string containing the full file path name to the location
        of the image containing the containing the longitude values.

    :param lat_fname:
        A string containing the full file path name to the location
        of the image containing the containing the latitude values.

    :param npoints:
        The number of time sample points to be calculated along the
        satellite track. Default is 12

    :param to_disk:
        If set to None (default) then the results will be returned
        in memory and no disk space will be used. Otherwise to_disk
        should be a list of length 6 containing file path names for
        the computed arrays. These arrays will be written directly
        to disk in a tiled fashion. Setting this keyword reduces
        memory consumption. When set, then 6 filepathnames will
        be returned instead of np arrays.
        The order is important, and is given as follows:
        Satellite zenith angle.
        Satellite azimuth angle.
        Solar zenith angle.
        Solar azimuth angle.
        Relative azimuth angle.
        Time.

    :return:
        6 float32 np arrays of the same shape as lon_array unless
        the outfilenames is set in which case 6 filepath
        names will be returned:

        1. Satellite zenith angle.
        2. Satellite azimuth angle.
        3. Solar zenith angle.
        4. Solar azimuth angle.
        5. Relative azimuth angle.
        6. Time.

        3 float32 np arrays with the same shape as
        lon_array.shape[0] containing the array coodinates of the
        satellite track line and N_Cent:

        4. Y_coordinates (Starting at 1)
        5. X_coordinates
        6. n_cent (Value 2 if centre x coordinate was averaged and
           1 if centre x coordinate was not averaged)
    """
    # Get the datetime of the acquisition
    dt = acquisition.scene_center_datetime

    # Compute the geobox
    geobox = gridded_geo_box(acquisition)

    # Image projection, geotransform and scene centre time stamp
    prj = geobox.crs.ExportToWkt()
    geoT = geobox.affine.to_gdal()

    # Get the lat/lon of the scene centre
    centre_xy = geobox.centre_lonlat

    # Get the earth spheroidal paramaters
    spheroid = setup_spheroid(prj)

    # Get the satellite orbital elements
    sat_ephemeral = load_tle(acquisition, TLE_DIR)

    # If we have None, then no suitable TLE was found, so use values gathered
    # by the acquisition object
    if sat_ephemeral is None:
        # orbital inclination (degrees)
        orb_incl = math.degrees(acquisition.inclination)
        # semi_major radius (m)
        orb_radius = acquisition.semi_major_axis
        # angular velocity (rad sec-1)
        omega = acquisition.omega
        orbital_elements = np.array([orb_incl, orb_radius, omega],
                                    dtype='float')
    else:
        orbital_elements = setup_orbital_elements(sat_ephemeral, dt)

    # Min and Max lat extents
    # This method should handle northern and southern hemispheres
    min_lat = min(min(geobox.ul_lonlat[1], geobox.ur_lonlat[1]),
                  min(geobox.ll_lonlat[1], geobox.lr_lonlat[1]))
    max_lat = max(max(geobox.ul_lonlat[1], geobox.ur_lonlat[1]),
                  max(geobox.ll_lonlat[1], geobox.lr_lonlat[1]))

    # Scene centre in time stamp in decimal hours
    hours = acquisition.decimal_hour

    # Calculate the julian century past JD2000
    century = calculate_julian_century(dt)

    # Need something to determine max satellite view angle
    # Currently not even used in Fuqin's code
    view_max = 9.0

    # Get the satellite model paramaters
    smodel = setup_smodel(
        centre_xy[0], centre_xy[1], spheroid, orbital_elements)

    # Get the times and satellite track information
    track = setup_times(
        min_lat, max_lat, spheroid, orbital_elements, smodel, npoints)

    # Array dimensions
    cols = acquisition.samples
    rows = acquisition.lines
    dims = (rows, cols)

    if to_disk is None:
        # Initialise 2D arrays to hold the angles
        view = np.zeros(dims, dtype='float32')
        azi = np.zeros(dims, dtype='float32')
        asol = np.zeros(dims, dtype='float32')
        soazi = np.zeros(dims, dtype='float32')
        rela_angle = np.zeros(dims, dtype='float32')
        time = np.zeros(dims, dtype='float32')
    else:
        # Initialise 1D arrays to hold the angles
        view = np.zeros((1, cols), dtype='float32')
        azi = np.zeros((1, cols), dtype='float32')
        asol = np.zeros((1, cols), dtype='float32')
        soazi = np.zeros((1, cols), dtype='float32')
        rela_angle = np.zeros((1, cols), dtype='float32')
        time = np.zeros((1, cols), dtype='float32')

        if len(to_disk) != 6:
            print "Incorrect number of filenames!"
            print "Results will be returned as np arrays"
            to_disk = None

        # Initialise the output files
        output_files = []
        drv = gdal.GetDriverByName("ENVI")
        output_files.append(drv.Create(to_disk[0], cols, rows, 1, 6))
        output_files.append(drv.Create(to_disk[1], cols, rows, 1, 6))
        output_files.append(drv.Create(to_disk[2], cols, rows, 1, 6))
        output_files.append(drv.Create(to_disk[3], cols, rows, 1, 6))
        output_files.append(drv.Create(to_disk[4], cols, rows, 1, 6))
        output_files.append(drv.Create(to_disk[5], cols, rows, 1, 6))

        # Set the projection and geotransfrom
        for outds in output_files:
            outds.SetProjection(prj)
            outds.SetGeoTransform(geoT)

        # Get the band level write access
        out_SAT_V_bnd = output_files[0].GetRasterBand(1)
        out_SAT_AZ_bnd = output_files[1].GetRasterBand(1)
        out_SOL_Z_bnd = output_files[2].GetRasterBand(1)
        out_SOL_AZ_bnd = output_files[3].GetRasterBand(1)
        out_REL_AZ_bnd = output_files[4].GetRasterBand(1)
        out_TIME_bnd = output_files[5].GetRasterBand(1)

    # Set to null value
    view[:] = -999
    azi[:] = -999
    asol[:] = -999
    soazi[:] = -999
    rela_angle[:] = -999
    time[:] = -999

    # Initialise centre line variables
    y_cent = np.arange(1, rows + 1).astype('float32')
    x_cent = np.zeros((rows), dtype='float32')
    n_cent = np.zeros((rows), dtype='float32')

    # Initialise the tile generator for processing
    # Process 1 row of data at a time
    tiles = tiling.generate_tiles(cols, rows, cols, 1, Generator=True)

    # Rather than do 8000+ if checks within the loop, we'll construct two
    # seperate loops
    if to_disk is None:
        with rasterio.open(lon_fname) as lon, rasterio.open(lat_fname) as lat:
            # Loop over each row
            for i in range(rows):
                # Get current tile
                tile = tiles.next()
                ystart = tile[0][0]
                xstart = tile[1][0]

                # Read the lon and lat tile
                lon_array = lon.read_band(1, window=tile)
                lat_array = lat.read_band(1, window=tile)

                istat = angle(cols, rows, i + 1, lat_array, lon_array,
                              spheroid, orbital_elements, hours, century,
                              npoints, smodel, track, view[i], azi[i],
                              asol[i], soazi[i], rela_angle[i], time[i],
                              x_cent, n_cent)
    else:
        with rasterio.open(lon_fname) as lon, rasterio.open(lat_fname) as lat:
            # Loop over each row
            for i in range(rows):
                # Get tile info
                tile = tiles.next()
                ystart = tile[0][0]
                xstart = tile[1][0]

                # Read the lon and lat tile
                lon_array = lon.read_band(1, window=tile)
                lat_array = lat.read_band(1, window=tile)

                # Set to null value
                view[:] = -999
                azi[:] = -999
                asol[:] = -999
                soazi[:] = -999
                rela_angle[:] = -999
                time[:] = -999

                istat = angle(cols, rows, i + 1, lat_array, lon_array,
                              spheroid, orbital_elements, hours, century,
                              npoints, smodel, track, view[0], azi[0],
                              asol[0], soazi[0], rela_angle[0], time[0],
                              x_cent, n_cent)

                # Output to disk
                out_SAT_V_bnd.WriteArray(view, xstart, ystart)
                out_SAT_V_bnd.FlushCache()
                out_SAT_AZ_bnd.WriteArray(azi, xstart, ystart)
                out_SAT_AZ_bnd.FlushCache()
                out_SOL_Z_bnd.WriteArray(asol, xstart, ystart)
                out_SOL_Z_bnd.FlushCache()
                out_SOL_AZ_bnd.WriteArray(soazi, xstart, ystart)
                out_SOL_AZ_bnd.FlushCache()
                out_REL_AZ_bnd.WriteArray(rela_angle, xstart, ystart)
                out_REL_AZ_bnd.FlushCache()
                out_TIME_bnd.WriteArray(time, xstart, ystart)
                out_TIME_bnd.FlushCache()

    if to_disk is not None:
        # Close all image files opened for writing
        output_files = None
        outds = None
        out_SAT_V_bnd = None
        out_SAT_AZ_bnd = None
        out_SOL_Z_bnd = None
        out_SOL_AZ_bnd = None
        out_REL_AZ_bnd = None
        out_TIME_bnd = None

    # Centreline
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
    # Account for first element potentially being changed with the
    # last element
    if wh[0]:
        x_cent[0] = temp[1]

    # Convert X & Y centre points to integers (basically array co-ordinates)
    y_cent = np.rint(y_cent)
    x_cent = np.rint(x_cent)

    # If we didn't write to disk return np arrays otherwise
    # return the filepath names
    if to_disk is None:
        result = (view, azi, asol, soazi, rela_angle, time, y_cent, x_cent,
                  n_cent)
        return result
    else:
        result = (to_disk[0], to_disk[1], to_disk[2], to_disk[3], to_disk[4],
                  to_disk[5], y_cent, x_cent, n_cent)
        return result
