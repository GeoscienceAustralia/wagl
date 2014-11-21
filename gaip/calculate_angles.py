#!/usr/bin/env python
import numpy
from osgeo import gdal
from osgeo import osr
import ephem

from EOtools.DatasetDrivers import SceneDataset

from ULA3.set_satmod import set_satmod
from ULA3.set_times import set_times
from ULA3.angle_all import angle

from gaip import find_file
from gaip import read_img


# To be used as a template while gaip is restructured
def sat_sol_grid_workflow(DATA, CONFIG):
    """
    Workflow to generate the satellite and solar grids.
    """

    # Track our working directory so we can save files
    work_dir = CONFIG.work_path

    # Find and open the longitude and lattitude files
    # Avoiding DataManger here. find_file will be used sparingly until a proper
    # workflow is written.
    lon_fname = find_file(work_dir, 'LON.tif')
    lat_fname = find_file(work_dir, 'LAT.tif')

    # We should be able to change the workflow and pass these as filename
    # strings.  Internally, the calculate_angles could read a row at a time
    # rather than pass the entire array. Internally angles are calculated one
    # row at a time.
    lon_arr = read_img(lon_fname)
    lat_arr = read_img(lat_fname)

    # Get the array dimensions
    dims = lon_arr.shape
    cols = dims[1]
    rows = dims[0]

    # Get the L1T data. (Code borrowed from above)
    L1T_dataset = DATA.get_item(CONFIG.input['l1t']['path'], SceneDataset)
    assert L1T_dataset, 'Unable to retrieve SceneDataset object for L1T input scene dataset'
    logger.debug( 'SceneDataset object for %s retrieved', L1T_dataset.pathname)

    # Define the output file names
    sat_view_zenith_fname  = os.path.join(work_dir, 'SAT_V.bin')
    sat_azimuth_fname      = os.path.join(work_dir, 'SAT_AZ.bin')
    solar_zenith_fname     = os.path.join(work_dir, 'SOL_Z.bin')
    solar_azimuth_fname    = os.path.join(work_dir, 'SOL_AZ.bin')
    relative_azimuth_fname = os.path.join(work_dir, 'REL_AZ.bin')
    time_fname             = os.path.join(work_dir, 'TIME.bin')

    out_fnames = [sat_view_zenith_fname, sat_azimuth_fname, solar_zenith_fname
        solar_azimuth_fname, relative_azimuth_fname, time_fname]

    # Get the angles, time, & satellite track coordinates
    (satellite_zenith, satellite_azimuth, solar_zenith,
     solar_azimuth, relative_azimuth, time,
     Y_cent, X_cent, N_cent) = ca.calculate_angles(L1T_dataset, lon_arr,
        lat_arr, npoints=12, to_disk=out_fnames)

    # Write out the CENTRELINE file
    create_centreline_file(Y_cent, X_cent, N_cent, cols, view_max=9.0,
                           outdir=work_dir)


def create_centreline_file(y, x, n, cols, view_max, outdir, outfname='CENTRELINE'):
    """
    :param y:
        A 1D NumPy array of type int with the same shape as x & n.
        Details the row number starting at 1.

    :param x:
        A 1D NumPy array of type int with the same shape as y & n.
        Details the column number starting at 0.

    :param n:
        A 1D NumPy array of type int with the same shaoe as y & x.
        Details whether or not the track point coordinate is
        averaged.

    :param cols:
        An integer indicating the number of columns in the original
        array.

    :param view_max:
        An float indicating the maximum view angle of the satellite
        used for determining the centreline.

    :param outdir:
        A string containing the output directory.

    :param outfname:
        A string containing the output filename.
        Default is CENTRELINE.
    """

    rows = y.shape[0]

    # I'm not sure of the formatting used by FORTRAN, nothing was specified
    # but the outputs had spaces.
    # It might be more ideal to created it as a csv???
    fname = os.path.join(outdir, outfname)
    outf  = open(fname, 'w')

    outf.write('%f\n' %view_max)
    outf.write('%i      %i\n' %(rows, cols))

    for r in range(rows):
        outf.write('%i      %i       %f     %f     %f\n'%(y[r], x[r], n[r], n[r], n[r]))

    outf.close()


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


def calculate_angles(scene_dataset, lon_array, lat_array, npoints=12,
        to_disk=None):
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

    :param to_disk:
        If set to None (default) then the results will be returned
        in memory and no disk space will be used. Otherwise to_disk
        should be a list of length 6 containing file path names for
        the computed arrays. These arrays will be written directly
        to disk in a tiled fashion. Setting this keyword reduces
        memory consumption. When set, then 6 filepathnames will
        be returned instead of NumPy arrays.
        The order is important, and is given as follows:
        Satellite zenith angle.
        Satellite azimuth angle.
        Solar zenith angle.
        Solar azimuth angle.
        Relative azimuth angle.
        Time.

    :return:
        6 float32 NumPy arrays of the same shape as lon_array unless
        the outfilenames is set in which case 6 filepath
        names will be returned:
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

    if to_disk is None:
        # Initialise 2D arrays to hold the angles
        view = numpy.zeros(dims, dtype='float32')
        azi = numpy.zeros(dims, dtype='float32')
        asol = numpy.zeros(dims, dtype='float32')
        soazi = numpy.zeros(dims, dtype='float32')
        rela_angle = numpy.zeros(dims, dtype='float32')
        time = numpy.zeros(dims, dtype='float32')
    else:
        # Initialise 1D arrays to hold the angles
        view = numpy.zeros(dims, dtype='float32')
        azi = numpy.zeros(dims, dtype='float32')
        asol = numpy.zeros(dims, dtype='float32')
        soazi = numpy.zeros(dims, dtype='float32')
        rela_angle = numpy.zeros(dims, dtype='float32')
        time = numpy.zeros(dims, dtype='float32')

        # Generate a list of tiles for processing
        # Process 1 row of data at a time
        tiles = tiling.generate_tiles(cols, rows, cols, 1)

        if len(to_disk) != 6:
            print "Incorrect number of filenames!"
            print "Results will be returned as NumPy arrays"
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
        out_SAT_V_bnd  = output_files[0].GetRasterBand(1)
        out_SAT_AZ_bnd = output_files[1].GetRasterBand(1)
        out_SOL_Z_bnd  = output_files[2].GetRasterBand(1)
        out_SOL_AZ_bnd = output_files[3].GetRasterBand(1)
        out_REL_AZ_bnd = output_files[4].GetRasterBand(1)
        out_TIME_bnd   = output_files[5].GetRasterBand(1)

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

    # Rather than do 8000+ if checks within the loop, we'll construct
    # different loops
    if to_disk is None:
        # Loop over each row
        for i in range(rows):
            istat = angle(cols, rows, i+1, lat_array[i], lon_array[i],
                spheroid, orbital_elements, hours, century, npoints, smodel,
                track, view[i], azi[i], asol[i], soazi[i], rela_angle[i],
                time[i], X_cent, N_cent)
    else:
        # Loop over each row
        for i in range(rows):
            # Set to null value
            view[:] = -999
            azi[:] = -999
            asol[:] = -999
            soazi[:] = -999
            rela_angle[:] = -999
            time[:] = -999

            # Get current tile
            tile = tiles[i]
            xstart = tile[2]
            ystart = tile[0]

            istat = angle(cols, rows, i+1, lat_array[i], lon_array[i],
                spheroid, orbital_elements, hours, century, npoints,
                smodel, track, view, azi, asol, soazi, rela_angle,
                time, X_cent, N_cent)

            # Output to disk
            out_SAT_V_bnd.WriteArray(view, xstart, ystart).FlushCache()
            out_SAT_AZ_bnd.WriteArray(azi, xstart, ystart).FlushCache()
            out_SOL_Z_bnd.WriteArray(asol, xstart, ystart).FlushCache()
            out_SOL_AZ_bnd.WriteArray(soazi, xstart, ystart).FlushCache()
            out_REL_AZ_bnd.WriteArray(rela_angle, xstart, ystart).FlushCache()
            out_TIME_bnd.WriteArray(time, xstart, ystart).FlushCache()

    if to_disk is not None:
        # Close all image files opened for writing
        output_files = None
        outds = None
        out_SAT_V_bnd  = None
        out_SAT_AZ_bnd = None
        out_SOL_Z_bnd  = None
        out_SOL_AZ_bnd = None
        out_REL_AZ_bnd = None
        out_TIME_bnd   = None

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

    # If we didn't write to disk return NumPy arrays otherwise
    # return the filepath names
    if to_disk is None:
        result = (view, azi, asol, soazi, rela_angle, time, Y_cent, X_cent,
            N_cent)
        return result
    else:
        result = (to_disk[0], to_disk[1], to_disk[2], to_disk[3], to_disk[4],
            to_disk[5], Y_cent, X_cent, N_cent)
        return result

