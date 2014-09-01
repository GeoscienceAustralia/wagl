#!/usr/bin/env python
import numpy
from osgeo import gdal
from osgeo import osr
import gc
import ephem

from EOtools.DatasetDrivers import SceneDataset
from image_tools import write_img

from set_satmod import set_satmod
from set_times import set_times
from angle_all import angle

import pdb

def calculate_julian_century(datetime):
    """
    
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

def get_lattitude_min_max():
    """
    
    """

def setup_spheroid(proj_wkt):
    """
    
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

def setup_orbital_elements():
    """
    
    """

    # For the moment hard code the values
    # Need to be able to get these from a TLE
    orbital_elements = numpy.zeros((3))

    # orbital inclination (degrees)
    orbital_elements[0] = 98.200000

    # semi_major radius (m)
    orbital_elements[1] = 7083160.000000

    # angular velocity (rad sec-1)
    orbital_elements[2] = 0.001059

    return orbital_elements

def setup_smodel(centre_lon, centre_lat, spheroid, orbital_elements):
    """
    
    """

    smodel, istat = set_satmod(centre_lon, centre_lat, spheroid, orbital_elements)

    return smodel

def setup_times(ymin, ymax, spheroid, orbital_elements, smodel, npoints=12):
    """
    
    """

    track, istat = set_times(ymin, ymax, npoints, spheroid, orbital_elements, smodel)

    return track

def calculate_angles(scene_dataset, lon_array, lat_array, npoints=12):
    """
    
    """

    prj = scene_dataset.GetProjection()
    geoT = scene_dataset.GetGeoTransform()

    centre_xy = numpy.array(scene_dataset.lonlats['CENTRE'])

    print "Setup the spheroid"
    spheroid = setup_spheroid(scene_dataset.GetProjection())

    print "Setup the orbital_elements"
    orbital_elements = setup_orbital_elements()

    lat_min_max = scene_dataset.get_bounds()[1]
    hours = scene_dataset.decimal_hour

    print "Calculate the julian century past JD2000"
    century = calculate_julian_century(scene_dataset.scene_centre_datetime)

    # Need something to determine max satellite view angle
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

    print "Writing the angle images"
    print "Writing view"
    write_img(view, 'view', projection=prj, geotransform=geoT)
    print "Writing azi"
    write_img(azi, 'azi', projection=prj, geotransform=geoT)
    print "Writing asol"
    write_img(asol, 'asol', projection=prj, geotransform=geoT)
    print "Writing soazi"
    write_img(soazi, 'soazi', projection=prj, geotransform=geoT)
    print "Writing rela_angle"
    write_img(rela_angle, 'rela_angle', projection=prj, geotransform=geoT)
    print "Writing time"
    write_img(time, 'time', projection=prj, geotransform=geoT)

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
    outf = open('centreline.txt', 'w')
    for r in range(rows):
        outf.write('%i, %i, %f\n'%(Y_cent[r], X_cent[r], N_cent[r]))

    outf.close()


if __name__ == '__main__':

    f = '/g/data1/v10/NBAR_validation_reference/Nov2013/L1T_Input/LS7_90-84_2013-10-03/EQR/LS7_ETM_OTH_P51_GALPGS04-002_090_084_20131003'
    lon_f = '/short/v10/jps547/nbar/test_skew/work/LON.tif'
    lat_f = '/short/v10/jps547/nbar/test_skew/work/LAT.tif'

    ds = SceneDataset(f)

    lon_ds = gdal.Open(lon_f)
    lat_ds = gdal.Open(lat_f)

    lon_arr = lon_ds.ReadAsArray()
    lat_arr = lat_ds.ReadAsArray()

    # Close the lat/long datasets
    lon_ds = None
    lat_ds = None

    print "Starting main routine"
    calculate_angles(ds, lon_arr, lat_arr)
    print "Finished"
