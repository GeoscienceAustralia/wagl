#!/usr/bin/env python

from gaip import Buffers
from gaip import calculate_angles as ca
from gaip import read_img
from gaip import run_slope

def calculate_self_shadow(acquisition, DSM_fname, buffer,
        solar_zenith_fname, solar_azimuth_fname,
        satellite_view_fname, satellite_azimuth_fname,
        out_fnames=None):
    """
    Computes the self shadow mask, slope, aspect, incident, exiting,
    azimuth incident, azimuth exiting and relative slope angles.

    :param acquisition:
        An instance of an acquisition object.

    :param DSM_fname:
        A string containing the full file path name to the Digital
        Surface Model to be used in deriving the surface angles.

    :param buffer:
        An object with members top, bottom, left and right giving the
        size of the buffer (in pixels) which have been added to the
        corresponding sides of DSM.

    :param solar_zenith_fname:
        A string containing the full file path name to the solar
        zenith angle image.

    :param solar_azimuth_fname:
        A string containing the full file path name to the solar
        azimuth angle image.

    :param satellite_view_fname:
        A string containing the full file path name to the satellite
        view angle image.

    :param satellite_azimuth_fname:
        A string containing the full file path name to the satellite
        azimuth angle image.

    :param out_fnames:
        A list of length 8 containing full file path names for each
        of the outputs in the following order:
            * Self shadow mask
            * Slope
            * Aspect
            * Incident angle
            * Exiting angle
            * Azimuth incident angle
            * Azimuth exiting angle
            * Relative slope
        Default is None, in which case default filenames will be
        used and the results written to the current working
        directory.

    :return:
        None. Outputs are written to disk.
    """
    # Setup the geobox
    geobox = acquisition.gridded_geo_box()

    # Retrive the spheroid parameters
    # (used in calculating pixel size in metres per lat/lon)
    spheroid = ca.setup_spheroid(geobox.crs.ExportToWkt())

    # Are we in projected or geographic space
    is_utm = not geobox.crs.IsGeographic()

    # Read the DSM and angle arrays into memory
    DSM = read_img(DSM_fname)
    solar_zenith = read_img(solar_zenith_fname)
    solar_azimuth = read_img(solar_azimuth_fname)
    satellite_view = read_img(satellite_view_fname)
    satellite_azimuth = read_img(satellite_azimuth_fname)

    # Define Top, Bottom, Left, Right pixel buffers
    pixel_buf = Buffers(buffer)

    # Compute self shadow, slope and various other angles
    slope_results = run_slope(acquisition, DSM, solar_zenith, satellite_view,
        solar_azimuth, satellite_azimuth, pixel_buf, is_utm, spheroid)

    # Output the results
    slope_results.write_arrays(out_fnames=out_fnames, geo_box=geobox)
