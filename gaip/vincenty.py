#!/usr/bin/env python

# ===============================================================================
# Copyright 2015 Geoscience Australia
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ===============================================================================

'''
Algorithms from Geocentric Datum of Australia Technical Manual

http://www.anzlic.org.au/icsm/gdatum/chapter4.html

This page last updated 11 May 1999

Computations on the Ellipsoid

There are a number of formulae that are available
to calculate accurate geodetic positions,
azimuths and distances on the ellipsoid.

Vincenty's formulae (Vincenty, 1975) may be used
for lines ranging from a few cm to nearly 20,000 km,
with millimetre accuracy.
The formulae have been extensively tested
for the Australian region, by comparison with results
from other formulae (Rainsford, 1955 & Sodano, 1965).

* Inverse problem: azimuth and distance from known
                    latitudes and longitudes
* Direct problem: Latitude and longitude from known
                    position, azimuth and distance.
* Sample data
* Excel spreadsheet

Vincenty's Inverse formulae
Given: latitude and longitude of two points
                    (phi1, lembda1 and phi2, lembda2),
Calculate: the ellipsoidal distance (s) and
forward and reverse azimuths between the points (alpha12, alpha21).
'''

from __future__ import absolute_import, print_function
import math

import numpy

__version__ = '1.0.1'


class GreatCircle(object):

    """
    formula for perfect sphere from Ed Williams' 'Aviation Formulary'
    (http://williams.best.vwh.net/avform.htm)

    code for ellipsoid posted to GMT mailing list by Jim Leven in Dec 1999

    Contact: Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
    """

    def __init__(self, rmajor, rminor, lon1, lat1, lon2, lat2):
        """
        Define a great circle by specifying:
        rmajor - radius of major axis of ellipsoid
        rminor - radius of minor axis of ellipsoid.
        lon1 - starting longitude of great circle
        lat1 - starting latitude
        lon2 - ending longitude
        lat2 - ending latitude
        All must be given in degrees.

        Instance variables:
        distance - distance along great circle in radians.
        lon1,lat1,lon2,lat2 - start and end points (in radians).
        """
        # convert to radians from degrees.
        lat1 = math.radians(lat1)
        lon1 = math.radians(lon1)
        lat2 = math.radians(lat2)
        lon2 = math.radians(lon2)
        self.a = rmajor
        self.f = (rmajor - rminor) / rmajor
        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2
        # distance along geodesic in meters.
        d, a12, a21 = vinc_dist(self.f, self.a, lat1, lon1, lat2, lon2)
        self.distance = d
        self.azimuth12 = a12
        self.azimuth21 = a21
        # great circle arc-length distance (in radians).
        self.gcarclen = 2. * math.asin(math.sqrt((math.sin((lat1 - lat2) / 2))**2 +
                                                 math.cos(lat1) * math.cos(lat2) * (math.sin((lon1 - lon2) / 2))**2))
        # check to see if points are antipodal (if so, route is undefined).
        if self.gcarclen == math.pi:
            self.antipodal = True
        else:
            self.antipodal = False

    def points(self, npoints):
        """
        compute arrays of npoints equally spaced
        intermediate points along the great circle.

        input parameter npoints is the number of points
        to compute.

        Returns lons, lats (lists with longitudes and latitudes
        of intermediate points in degrees).

        For example npoints=10 will return arrays lons,lats of 10
        equally spaced points along the great circle.
        """
        # must ask for at least 2 points.
        if npoints <= 1:
            raise ValueError('npoints must be greater than 1')
        elif npoints == 2:
            return [
                math.degrees(
                    self.lon1), math.degrees(
                    self.lon2)], [
                math.degrees(
                    self.lat1), math.degrees(
                        self.lat2)]
        # can't do it if endpoints are antipodal, since
        # route is undefined.
        if self.antipodal:
            raise ValueError('cannot compute intermediate points on a great circle whose endpoints are antipodal')
        d = self.gcarclen
        delta = 1.0 / (npoints - 1)
        f = delta * numpy.arange(npoints)  # f=0 is point 1, f=1 is point 2.
        incdist = self.distance / (npoints - 1)
        lat1 = self.lat1
        lat2 = self.lat2
        lon1 = self.lon1
        lon2 = self.lon2
        # perfect sphere, use great circle formula
        if self.f == 0.:
            A = numpy.sin((1 - f) * d) / math.sin(d)
            B = numpy.sin(f * d) / math.sin(d)
            x = A * math.cos(lat1) * math.cos(lon1) + B * math.cos(lat2) * math.cos(lon2)
            y = A * math.cos(lat1) * math.sin(lon1) + B * math.cos(lat2) * math.sin(lon2)
            z = A * math.sin(lat1) + B * math.sin(lat2)
            lats = numpy.arctan2(z, numpy.sqrt(x**2 + y**2))
            lons = numpy.arctan2(y, x)
            lons = map(math.degrees, lons.tolist())
            lats = map(math.degrees, lats.tolist())
        # use ellipsoid formulas
        else:
            latpt = self.lat1
            lonpt = self.lon1
            azimuth = self.azimuth12
            lons = [math.degrees(lonpt)]
            lats = [math.degrees(latpt)]
            for n in range(npoints - 2):
                latptnew, lonptnew, alpha21 = vinc_pt(self.f, self.a, latpt, lonpt, azimuth, incdist)
                d, azimuth, a21 = vinc_dist(self.f, self.a, latptnew, lonptnew, lat2, lon2)
                lats.append(math.degrees(latptnew))
                lons.append(math.degrees(lonptnew))
                latpt = latptnew
                lonpt = lonptnew
            lons.append(math.degrees(self.lon2))
            lats.append(math.degrees(self.lat2))
        return lons, lats

# ---------------------------------------------------------------------
# |                                                                    |
# |     geodetic.py -  a collection of geodetic functions              |
# |                                                                    |
# ---------------------------------------------------------------------
#
#
# ----------------------------------------------------------------------
# | Algrothims from Geocentric Datum of Australia Technical Manual      |
# |                                                                     |
# | http://www.anzlic.org.au/icsm/gdatum/chapter4.html                  |
# |                                                                     |
# | This page last updated 11 May 1999                                  |
# |                                                                     |
# | Computations on the Ellipsoid                                       |
# |                                                                     |
# | There are a number of formulae that are available                   |
# | to calculate accurate geodetic positions,                           |
# | azimuths and distances on the ellipsoid.                            |
# |                                                                     |
# | Vincenty's formulae (Vincenty, 1975) may be used                    |
# | for lines ranging from a few cm to nearly 20,000 km,                |
# | with millimetre accuracy.                                           |
# | The formulae have been extensively tested                           |
# | for the Australian region, by comparison with results               |
# | from other formulae (Rainsford, 1955 & Sodano, 1965).               |
# |                                                                     |
# | * Inverse problem: azimuth and distance from known                  |
# |                     latitudes and longitudes                        |
# | * Direct problem: Latitude and longitude from known                 |
# |                     position, azimuth and distance.                 |
# | * Sample data                                                       |
# | * Excel spreadsheet                                                 |
# |                                                                     |
# | Vincenty's Inverse formulae                                         |
# | Given: latitude and longitude of two points                         |
# |                     (phi1, lembda1 and phi2, lembda2),              |
# | Calculate: the ellipsoidal distance (s) and                         |
# | forward and reverse azimuths between the points (alpha12, alpha21). |
# |                                                                     |
# ----------------------------------------------------------------------


def vinc_dist(f, a, phi1, lembda1, phi2, lembda2):
    """
    Returns the distance between two geographic points on the ellipsoid
    and the forward and reverse azimuths between these points.
    lats, longs and azimuths are in radians, distance in metres

    Arguments:
        f: flattening
        a: equatorial radius (metres)
        phi1: latitude of first point
        lembda1: longitude of first point
        phi2: latitude of second point
        lembda2: longitude of second point

    Returns ( s, alpha12,  alpha21 ) as a tuple
    """

    if (abs(phi2 - phi1) < 1e-8) and (abs(lembda2 - lembda1) < 1e-8):
        return 0.0, 0.0, 0.0

    two_pi = 2.0 * math.pi

    b = a * (1.0 - f)

    TanU1 = (1 - f) * math.tan(phi1)
    TanU2 = (1 - f) * math.tan(phi2)

    U1 = math.atan(TanU1)
    U2 = math.atan(TanU2)

    lembda = lembda2 - lembda1
    last_lembda = -4000000.0                # an impossibe value
    omega = lembda

    # Iterate the following equations,
    #  until there is no significant change in lembda

    while (last_lembda < -3000000.0 or lembda != 0 and abs((last_lembda - lembda) / lembda) > 1.0e-9):

        sqr_sin_sigma = pow(math.cos(U2) * math.sin(lembda), 2) + \
            pow((math.cos(U1) * math.sin(U2) -
                 math.sin(U1) * math.cos(U2) * math.cos(lembda)), 2)

        Sin_sigma = math.sqrt(sqr_sin_sigma)

        Cos_sigma = math.sin(U1) * math.sin(U2) + math.cos(U1) * math.cos(U2) * math.cos(lembda)

        sigma = math.atan2(Sin_sigma, Cos_sigma)

        Sin_alpha = math.cos(U1) * math.cos(U2) * math.sin(lembda) / math.sin(sigma)
        alpha = math.asin(Sin_alpha)

        Cos2sigma_m = math.cos(sigma) - (2 * math.sin(U1) * math.sin(U2) / pow(math.cos(alpha), 2))

        C = (f / 16) * pow(math.cos(alpha), 2) * (4 + f * (4 - 3 * pow(math.cos(alpha), 2)))

        last_lembda = lembda

        lembda = omega + (1 - C) * f * math.sin(alpha) * \
                         (
                             sigma +
                             C * math.sin(sigma) *
                             (Cos2sigma_m + C * math.cos(sigma) * (-1 + 2 * pow(Cos2sigma_m, 2)))
                         )

    u2 = pow(math.cos(alpha), 2) * (a * a - b * b) / (b * b)

    A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))

    B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))

    delta_sigma = B * Sin_sigma * (Cos2sigma_m + (B / 4) *
                                   (Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2)) -
                                    (B / 6) * Cos2sigma_m * (-3 + 4 * sqr_sin_sigma) *
                                    (-3 + 4 * pow(Cos2sigma_m, 2))))

    s = b * A * (sigma - delta_sigma)

    alpha12 = math.atan2((math.cos(U2) * math.sin(lembda)),
                         (math.cos(U1) * math.sin(U2) - math.sin(U1) * math.cos(U2) * math.cos(lembda)))

    alpha21 = math.atan2((math.cos(U1) * math.sin(lembda)),
                         (-math.sin(U1) * math.cos(U2) + math.cos(U1) * math.sin(U2) * math.cos(lembda)))

    if (alpha12 < 0.0):
        alpha12 = alpha12 + two_pi
    if (alpha12 > two_pi):
        alpha12 = alpha12 - two_pi

    alpha21 = alpha21 + two_pi / 2.0
    if (alpha21 < 0.0):
        alpha21 = alpha21 + two_pi
    if (alpha21 > two_pi):
        alpha21 = alpha21 - two_pi

    return s, alpha12, alpha21

    # END of Vincenty's Inverse formulae


# ----------------------------------------------------------------------------
# Vincenty's Direct formulae                                                |
# Given: latitude and longitude of a point (phi1, lembda1) and              |
# the geodetic azimuth (alpha12)                                            |
# and ellipsoidal distance in metres (s) to a second point,                 |
#                                                                           |
# Calculate: the latitude and longitude of the second point (phi2, lembda2) |
# and the reverse azimuth (alpha21).                                        |
#                                                                           |
# ----------------------------------------------------------------------------

def vinc_pt(f, a, phi1, lembda1, alpha12, s):
    """
    Returns the lat and long of projected point and reverse azimuth
    given a reference point and a distance and azimuth to project.
    lats, longs and azimuths are passed in decimal degrees

    Returns ( phi2,  lambda2,  alpha21 ) as a tuple
    """

    two_pi = 2.0 * math.pi

    if (alpha12 < 0.0):
        alpha12 = alpha12 + two_pi
    if (alpha12 > two_pi):
        alpha12 = alpha12 - two_pi

    b = a * (1.0 - f)

    TanU1 = (1 - f) * math.tan(phi1)
    U1 = math.atan(TanU1)
    sigma1 = math.atan2(TanU1, math.cos(alpha12))
    Sinalpha = math.cos(U1) * math.sin(alpha12)
    cosalpha_sq = 1.0 - Sinalpha * Sinalpha

    u2 = cosalpha_sq * (a * a - b * b) / (b * b)
    A = 1.0 + (u2 / 16384) * (4096 + u2 * (-768 + u2 *
                                           (320 - 175 * u2)))
    B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))

    # Starting with the approximation
    sigma = (s / (b * A))

    last_sigma = 2.0 * sigma + 2.0  # something impossible

    # Iterate the following three equations
    # until there is no significant change in sigma
    # two_sigma_m , delta_sigma
    while (abs((last_sigma - sigma) / sigma) > 1.0e-9):
        two_sigma_m = 2 * sigma1 + sigma

        delta_sigma = B * math.sin(sigma) * (
            math.cos(two_sigma_m) +
            (B / 4) * (
                math.cos(sigma) *
                (
                    -1 + 2 * math.pow(math.cos(two_sigma_m), 2) -
                    (B / 6) * math.cos(two_sigma_m) *
                    (-3 + 4 * math.pow(math.sin(sigma), 2)) *
                    (-3 + 4 * math.pow(math.cos(two_sigma_m), 2))
                )
            )
        )

        last_sigma = sigma
        sigma = (s / (b * A)) + delta_sigma

    phi2 = math.atan2((math.sin(U1) *
                       math.cos(sigma) +
                       math.cos(U1) *
                       math.sin(sigma) *
                       math.cos(alpha12)), ((1 -
                                             f) *
                                            math.sqrt(math.pow(Sinalpha, 2) +
                                                      pow(math.sin(U1) *
                                                          math.sin(sigma) -
                                                          math.cos(U1) *
                                                          math.cos(sigma) *
                                                          math.cos(alpha12), 2))))

    lembda = math.atan2((math.sin(sigma) * math.sin(alpha12)), (math.cos(U1) * math.cos(sigma) -
                                                                math.sin(U1) * math.sin(sigma) * math.cos(alpha12)))

    C = (f / 16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq))

    omega = lembda - (1 - C) * f * Sinalpha *  \
        (sigma + C * math.sin(sigma) * (math.cos(two_sigma_m) +
                                        C * math.cos(sigma) * (-1 + 2 * math.pow(math.cos(two_sigma_m), 2))))

    lembda2 = lembda1 + omega

    alpha21 = math.atan2(Sinalpha, (-math.sin(U1) * math.sin(sigma) +
                                    math.cos(U1) * math.cos(sigma) * math.cos(alpha12)))

    alpha21 = alpha21 + two_pi / 2.0
    if (alpha21 < 0.0):
        alpha21 = alpha21 + two_pi
    if (alpha21 > two_pi):
        alpha21 = alpha21 - two_pi

    return phi2, lembda2, alpha21

    # END of Vincenty's Direct formulae

# #---------------------------------------------------------------------------
# Notes:
#
# * "The inverse formulae may give no solution over a line
#       between two nearly antipodal points. This will occur when
#       lembda ... is greater than pi in absolute value". (Vincenty, 1975)
#
# * In Vincenty (1975) L is used for the difference in longitude,
#       however for consistency with other formulae in this Manual,
#       omega is used here.
#
# * Variables specific to Vincenty's formulae are shown below,
#       others common throughout the manual are shown in the Glossary.
#
#
# alpha = Azimuth of the geodesic at the equator
# U = Reduced latitude
# lembda = Difference in longitude on an auxiliary sphere (lembda1 & lembda2
#               are the geodetic longitudes of points 1 & 2)
# sigma = Angular distance on a sphere, from point 1 to point 2
# sigma1 = Angular distance on a sphere, from the equator to point 1
# sigma2 = Angular distance on a sphere, from the equator to point 2
# sigma_m = Angular distance on a sphere, from the equator to the
#               midpoint of the line from point 1 to point 2
# u, A, B, C = Internal variables
#
#
# Sample Data
#
# Flinders Peak
# -37o57'03.72030"
# 144o25'29.52440"
# Buninyong
# -37o39'10.15610"
# 143o55'35.38390"
# Ellipsoidal Distance
# 54,972.271 m
#
# Forward Azimuth
# 306o52'05.37"
#
# Reverse Azimuth
# 127o10'25.07"
#
#
# #*******************************************************************
