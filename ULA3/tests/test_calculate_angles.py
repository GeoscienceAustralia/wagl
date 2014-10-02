#!/usr/bin/env python
import sys
import os
import argparse
import gc
import unittest

import numpy
import numpy.testing as npt
from osgeo import gdal
from osgeo import osr
import ephem

from EOtools.DatasetDrivers import SceneDataset

from set_satmod import set_satmod
from set_times import set_times
from angle_all import angle

import pdb

def write_img(array, name='', format='ENVI', projection=None, geotransform=None):
    """
    Just a small function to write float32 2D arrays.
    """

    dims = array.shape
    dtype = 6
    bands = 1

    drv = gdal.GetDriverByName(format)
    outds = drv.Create(name, dims[1], dims[0], bands, dtype)

    if prj:
        outds.SetProjection(prj)

    if geot:
        outds.SetGeoTransform(geot)

    band = outds.GetRasterBand(1)
    band.WriteArray(array)
    band.FlushCache()

    outds = None

def calculate_julian_century(datetime):
    """
    Given a datetime object return the julian century from the 2000 epoch.
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

    return orbital_elements

def setup_smodel(centre_lon, centre_lat, spheroid, orbital_elements):
    """
    Setup the satellite model.
    """

    smodel, istat = set_satmod(centre_lon, centre_lat, spheroid, orbital_elements)

    return smodel

def setup_times(ymin, ymax, spheroid, orbital_elements, smodel, npoints=12):
    """
    Setup the satellite track times.
    """

    track, istat = set_times(ymin, ymax, npoints, spheroid, orbital_elements, smodel)

    return track

def calculate_angles(scene_dataset, lon_array, lat_array, npoints=12):
    """
    Calcualte the satellite view, satellite azimuth, solar zenith,
    solar azimuth, and relative aziumth angle grids, as well as the
    time grid. All grids are output as float32 ENVI files.
    The CENTRELINE is also calculated and output to disk.
    """

    prj = scene_dataset.GetProjection()
    geoT = scene_dataset.GetGeoTransform()
    centre_datetime = scene_dataset.scene_centre_datetime

    centre_xy = numpy.array(scene_dataset.lonlats['CENTRE'])

    print "Setup the spheroid"
    spheroid = setup_spheroid(scene_dataset.GetProjection())

    print "Setup the orbital_elements"
    tle_dir = '/g/data1/v10/eoancillarydata/sensor-specific'
    sat_ephemeral = scene_dataset.satellite.load_tle(centre_datetime, tle_dir)
    orbital_elements = setup_orbital_elements(sat_ephemeral, centre_datetime)

    lat_min_max = scene_dataset.get_bounds()[1]
    hours = scene_dataset.decimal_hour

    print "Calculate the julian century past JD2000"
    century = calculate_julian_century(centre_datetime)

    # Need something to determine max satellite view angle
    # Currently not even used in Fuqin's code
    view_max = 9.0

    print "Setup the smodel"
    smodel = setup_smodel(centre_xy[0], centre_xy[1], spheroid, orbital_elements)

    print "Setup the times and track info"
    track = setup_times(lat_min_max[0], lat_min_max[1], spheroid, orbital_elements, smodel, npoints)

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

    print "Calculting the angles by looping over the rows"
    # Loop over each row
    for i in range(rows):
        istat = angle(cols, rows, i+1, lat_array[i], lon_array[i], spheroid,
                      orbital_elements, hours, century, npoints, smodel, track,
                      view[i], azi[i], asol[i], soazi[i], rela_angle[i],
                      time[i], X_cent, N_cent)

    # Angle file names
    sat_view_zenith_fname = 'SAT_V.bin'
    sat_azimuth_fname = 'SAT_AZ.bin'
    solar_zenith_fname = 'SOL_Z.bin'
    solar_azimuth_fname = 'SOL_AZ.bin'
    relative_azimuth_fname = 'REL_AZ.bin'
    time_fname = 'TIME.bin'

    print "Writing the angle images"

    print "Writing satelite view zenith angle: %s" %sat_view_zenith_fname
    write_img(view, sat_view_zenith_fname, projection=prj, geotransform=geoT)

    print "Writing satellite azimuth angle: %s" %sat_azimuth_fname
    write_img(azi, sat_azimuth_fname, projection=prj, geotransform=geoT)

    print "Writing solar zenith angle: %s" %solar_zenith_fname
    write_img(asol, solar_zenith_fname, projection=prj, geotransform=geoT)

    print "Writing solar azimith angle: %s" %solar_azimuth_fname
    write_img(soazi, solar_azimuth_fname, projection=prj, geotransform=geoT)

    print "Writing relative azimuth angle: %s" %relative_azimuth_fname
    write_img(rela_angle, relative_azimuth_fname, projection=prj, geotransform=geoT)

    print "Writing time array: %s" %time_fname
    write_img(time, time_fname, projection=prj, geotransform=geoT)

    del view, azi, asol, soazi, rela_angle, time
    gc.collect()

    print "Computing the centreline"
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

    print "Writing out the centreline file"
    # Write the centreline to disk
    outf = open('CENTRELINE', 'w')

    outf.write('%f\n'%view_max)
    outf.write('%i, %i\n'%(rows, cols))
    for r in range(rows):
        outf.write('%i, %i, %f\n'%(Y_cent[r], X_cent[r], N_cent[r]))

    outf.close()

#class AnglesOutputsTester(unittest.TestCase):
class AnglesOutputsTester(npt.TestCase):
    """
    Unittesting for the satellite and solar angle calculations.
    """

    def __init__(self, reference_dir, test_dir, decimal_precision=4,
                 integer_precision=1):
        """
        Unittests will occur for the following files:
        SAT_V.bin
        SAT_AZ.bin
        SOL_Z.bin
        SOL_AZ.bin
        REL_AZ.bin
        TIME.bin
        CENTRELINE

        :param reference_dir:
            A full file pathname to the directory containing the
            reference data.

        :param test_dir:
            A full file pathname to the directory containing the
            test data.

        :param decimal_precision:
            The decimal precision to be used during array comparison.
            Default is 4, i.e. values must be correct up to 4 d.p. in
            order to Pass.

        :param integer_precision:
            The intger precision to be used during array comparison.
            Default is 1, i.e. values must be correct within 1 integer
            in order to Pass.
        """

        self.reference_dir = reference_dir
        self.test_dir = test_dir
        self.dec_precision = decimal_precision
        self.int_precision = integer_precision

        self.fname_sat_v = 'SAT_V.bin'
        self.fname_sat_az = 'SAT_AZ.bin'
        self.fname_sol_z = 'SOL_Z.bin'
        self.fname_sol_az = 'SOL_AZ.bin'
        self.fname_rel_az = 'REL_AZ.bin'
        self.fname_time = 'TIME.bin'
        self.fname_centreline = 'CENTRELINE'

        # Read and store the centreline data in memory
        self.centreline_ref = None
        self.centreline_test = None
        self._read_centreline_files()

    def find_file(self, dir, file):
        """
        
        """
        fname = os.path.join(dir, file)
        if os.path.isfile(fname):
            return fname
        else:
            print "Error! file not found: %s" %fname
            sys.exit(2)

    def _read_centreline_files(self):
        """
        Read the centreline reference and test text files.
        """

        # Get the filenames for both the reference and test files
        ref_fname = self.find_file(self.reference_dir, self.fname_centreline)
        test_fname = self.find_file(self.test_dir, self.fname_centreline)

        # Open and read the reference data
        f = open(ref_fname)
        ref_data = f.readlines()
        f.close()

        # Open and read the test data
        f = open(test_fname)
        test_data = f.readlines()
        f.close()

        self.centreline_ref = ref_data
        self.centreline_test = test_data

    def test_satellite_view(self):
        """
        Test the satellite view angle array.
        """

        # Get the filenames for both the reference and test files
        ref_fname = self.find_file(self.reference_dir, self.fname_sat_v)
        test_fname = self.find_file(self.test_dir, self.fname_sat_v)

        # Open the image datasets
        ref_ds = gdal.Open(ref_fname)
        test_ds = gdal.Open(test_fname)

        # Get the image data
        ref_img = ref_ds.ReadAsArray()
        test_img = test_ds.ReadAsArray()

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=self.dec_precision))

    def test_satellite_azimuth(self):
        """
        Test the satellite azimuth angle array.
        """

        # Get the filenames for both the reference and test files
        ref_fname = self.find_file(self.reference_dir, self.fname_sat_az)
        test_fname = self.find_file(self.test_dir, self.fname_sat_az)

        # Open the image datasets
        ref_ds = gdal.Open(ref_fname)
        test_ds = gdal.Open(test_fname)

        # Get the image data
        ref_img = ref_ds.ReadAsArray()
        test_img = test_ds.ReadAsArray()

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=self.dec_precision))

    def test_solar_zenith(self):
        """
        Test the solar zenith angle array.
        """

        # Get the filenames for both the reference and test files
        ref_fname = self.find_file(self.reference_dir, self.fname_sol_z)
        test_fname = self.find_file(self.test_dir, self.fname_sol_z)

        # Open the image datasets
        ref_ds = gdal.Open(ref_fname)
        test_ds = gdal.Open(test_fname)

        # Get the image data
        ref_img = ref_ds.ReadAsArray()
        test_img = test_ds.ReadAsArray()

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=self.dec_precision))

    def test_solar_azimuth(self):
        """
        Test the solar azimuth angle array.
        """

        # Get the filenames for both the reference and test files
        ref_fname = self.find_file(self.reference_dir, self.fname_sol_az)
        test_fname = self.find_file(self.test_dir, self.fname_sol_az)

        # Open the image datasets
        ref_ds = gdal.Open(ref_fname)
        test_ds = gdal.Open(test_fname)

        # Get the image data
        ref_img = ref_ds.ReadAsArray()
        test_img = test_ds.ReadAsArray()

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=self.dec_precision))

    def test_relative_azimuth(self):
        """
        Test the relative azimuth angle array.
        """

        # Get the filenames for both the reference and test files
        ref_fname = self.find_file(self.reference_dir, self.fname_rel_az)
        test_fname = self.find_file(self.test_dir, self.fname_rel_az)

        # Open the image datasets
        ref_ds = gdal.Open(ref_fname)
        test_ds = gdal.Open(test_fname)

        # Get the image data
        ref_img = ref_ds.ReadAsArray()
        test_img = test_ds.ReadAsArray()

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=self.dec_precision))

    def test_time_array(self):
        """
        Test the time array.
        """

        # Get the filenames for both the reference and test files
        ref_fname = self.find_file(self.reference_dir, self.fname_time)
        test_fname = self.find_file(self.test_dir, self.fname_time)

        # Open the image datasets
        ref_ds = gdal.Open(ref_fname)
        test_ds = gdal.Open(test_fname)

        # Get the image data
        ref_img = ref_ds.ReadAsArray()
        test_img = test_ds.ReadAsArray()

        self.assertIsNone(npt.assert_almost_equal(test_img, ref_img,
                                                  decimal=self.dec_precision))

    def test_centreline_max_view_angle(self):
        """
        Test the maximum satellite view angle.
        This is the first line of the centreline file.
        """

        ref_data = float(self.centreline_ref[0].split())
        test_data = float(self.centreline_test[0].split())

        self.assertEqual(ref_data, test_data)

    def test_centreline_rows(self):
        """
        Test the number of rows.
        This is the first element of the second line in the
        centreline file.
        """

        ref_rows = int(self.centreline_ref[1].split[0])
        test_rows = int(self.centreline_test[1].split[0])

        self.assertEqual(ref_rows, test_rows)

    def test_centreline_columns(self):
        """
        Test the number of columns.
        This is the second element of the second line in the
        centreline file.
        """

        ref_cols = int(self.centreline_ref[1].split[1])
        test_cols = int(self.centreline_test[1].split[1])

        self.assertEqual(ref_cols, test_cols)

    def test_number_centreline_points(self):
        """
        Test that the number of centreline points are same.
        These start at the third line of the centreline file and
        contain 3 elements.
        """

        ref_data = self.centreline_ref[2:]
        test_data = self.centreline_test[2:]

        self.assertEqual(len(ref_data), len(test_data))

    def test_centreline_points(self):
        """
        Test that the centreline points are roughly the same.
        These start at the third line of the centreline file and
        contain 3 elements. We'll allow for variation by 1 integer.

        eg ['1742','4624','1.00000']
        """

        ref_data = self.centreline_ref[2:]
        test_data = self.centreline_test[2:]

        int_prec = self.integer_precision

        ref_points = numpy.zeros((len(ref_data), 3), dtype='int')
        test_points = numpy.zeros((len(ref_data), 3))

        for i in range(len(ref_data)):
            rx, ry, rz = ref_data[i].split()
            tx, ty, tz = test_data[i].split()
            ref_points[i,0], ref_points[i,1], ref_points[i,2] = rx, ry, float(rz)
            test_points[i,0], test_points[i,1], test_points[i,2] = tx, ty, float(tz)

        self.assertIsNone(npt.assert_allclose(test_points, ref_points, atol=int_prec))

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Calculates satellite and solar angles.')

    parser.add_argument('--L1T_dir', required=True, help='A directory path of a L1T scene.')
    parser.add_argument('--nbar_work_dir', required=True, help='A directory path to the associated NBAR working directory.')
    parser.add_argument('--outdir', required=True, help='A directory path that will contain the output files.')
    parser.add_argument('--dec_precision', default=4, help='The decimal precision used for array comparison')
    parser.add_argument('--int_precision', default=1, help='The integer precision used for array comparison')

    parsed_args = parser.parse_args()

    #f = '/g/data1/v10/NBAR_validation_reference/Nov2013/L1T_Input/LS7_90-84_2013-10-03/EQR/LS7_ETM_OTH_P51_GALPGS04-002_090_084_20131003'
    #lon_f = '/short/v10/jps547/nbar/test_skew/work/LON.tif'
    #lat_f = '/short/v10/jps547/nbar/test_skew/work/LAT.tif'

    L1T_dir = parsed_args.L1T_dir
    nbar_work_dir = parsed_args.nbar_work_dir
    outdir = parsed_args.outdir
    dec_precision = parsed_args.dec_precision
    int_precision = parsed_args.int_precision

    # Check the output directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Change to the output directory that will contain the results
    os.chdir(outdir)

    # Open the L1T dataset
    ds = SceneDataset(L1T_dir)

    # Find and open the longitude and lattitude files
    lon_f = os.path.join(nbar_work_dir, 'LON.tif')
    lat_f = os.path.join(nbar_work_dir, 'LAT.tif')

    if os.path.isfile(lon_f):
        lon_ds = gdal.Open(lon_f)
    else:
        print "Error! No longitude file could be found:\n %s"%lon_f
        sys.exit(2)
    
    if os.path.isfile(lat_f):
        lat_ds = gdal.Open(lat_f)
    else:
        print "Error! No lattitude file could be found:\n %s"%lat_f
        sys.exit(2)

    #lon_ds = gdal.Open(lon_f)
    #lat_ds = gdal.Open(lat_f)

    lon_arr = lon_ds.ReadAsArray()
    lat_arr = lat_ds.ReadAsArray()

    # Close the lat/long datasets
    lon_ds = None
    lat_ds = None

    print "Starting main routine"
    calculate_angles(ds, lon_arr, lat_arr)

    print "Running unittests on the angle grids and the CENTRELINE file."
    suite = AnglesOutputsTester(reference_dir=nbar_work_dir, test_dir=outdir,
                                decimal_precision=dec_precision,
                                integer_precision=int_precision)

    # Run each test
    #suite.test_satellite_view()
    #suite.test_satellite_azimuth()
    #suite.test_solar_zenith()
    #suite.test_solar_azimuth()
    #suite.test_relative_azimuth()
    ##suite.test_time_array() #Currently no time array output from Fuqins code
    #suite.test_centreline_max_view_angle()
    #suite.test_centreline_rows()
    #suite.test_centreline_columns()
    #suite.test_number_centreline_points()
    #suite.test_centreline_points()

    Suite = unittest.TestSuite()
    Suite.addTest(suite)
    unittest.TextTestRunner(verbosity=2).run(Suite)

    print "Finished"
