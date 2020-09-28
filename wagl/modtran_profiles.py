#!/usr/bin/env python
# pragma pylint: disable=line-too-long
# flake8: noqa
"""
These are the templates for the inputs into MODTRAN.
Parameters to be inserted are:
    * albedo
    * water vapour; units: g/cm^2
    * ozone; units: ATM-cm
    * filter function (spectral response of sensor)
    * aerosol visibility; optical depth at 550 nm (negative); units: km
    * elevation (sourced from a DSM); units: km above sea level
    * satellite height; units: km
    * satellite view angle (zenith); corrected (180 - angle); units: degrees
    * julian day of year
    * latitude; units: degrees
    * longitude; corrected (360 - angle); units: degrees
    * time; units: decimal hours in UTC
    * satellite azimuth angle; corrected (angle + 180); units: degrees
    * satellite view offset; (180 - angle); units: degrees
"""

from __future__ import absolute_import, print_function


MIDLAT_SUMMER_ALBEDO = """\
TM{binary} 2    2    2    1    2    2    2    2    2    2    1    1    0  10.000{albedo:7.2f}
TFF  8   0   375.000  g{water:7.5f}    a{ozone:5.3f}     T f f          t      0.3         0      0.70         0         0
{filter_function:<75}
    1    0    0    0    0    0{visibility:10.3f}     0.000     0.000     0.000{elevation:10.3f}
{sat_height:10.3f}{elevation:10.3f}{sat_view:10.3f}     0.000     0.000     0.000    0          0.000
    1    0{doy:5d}    0
{lat:10.3f}{lon:10.3f}     0.000     0.000{time:10.3f}{sat_azimuth:10.3f}     0.000     0.667
     350.0    2600.0       1.0       1.0RN#       NT    T                                                     
    0                                                                                                         
"""

TROPICAL_ALBEDO = """\
TM{binary} 1    2    2    1    1    1    1    1    1    1    1    1    0  10.000{albedo:7.2f}
TFF  8   0   375.000  g{water:7.5f}    a{ozone:5.3f}     T f f          t      0.3         0      0.70         0         0
{filter_function:<75}
    1    0    0    0    0    0{visibility:10.3f}     0.000     0.000     0.000{elevation:10.3f}
{sat_height:10.3f}{elevation:10.3f}{sat_view:10.3f}     0.000     0.000     0.000    0          0.000
    1    0{doy:5d}    0
{lat:10.3f}{lon:10.3f}     0.000     0.000{time:10.3f}{sat_azimuth:10.3f}     0.000     0.667
     350.0    2600.0       1.0       1.0RN#       NT    T                                                     
    0                                                                                                         
"""

MIDLAT_SUMMER_TRANSMITTANCE = """\
TM{binary} 2    2    2    1    2    2    2    2    2    2    1    1    0  10.000{albedo:7.2f}
TFF  8   0   375.000  g{water:7.5f}    a{ozone:5.3f}     T f f          t      0.3         0      0.70         0         0
{filter_function:<75}
    1    0    0    0    0    0{visibility:10.5f}     0.000     0.000     0.000{elevation:10.3f}
{sat_height:10.3f}{elevation:10.3f}{sat_view:10.3f}     0.000     0.000     0.000    0          0.000
    2    0{doy:5d}    0
     0.000{sat_view_offset:10.3f}                                                       0.667
     350.0    2600.0       1.0       1.0RN#       NT    T                                                     
    0                                                                                                         
"""

TROPICAL_TRANSMITTANCE = """\
TM{binary} 1    2    2    1    1    1    1    1    1    1    1    1    0  10.000{albedo:7.2f}
TFF  8   0   375.000  g{water:7.5f}    a{ozone:5.3f}     T f f          t      0.3         0      0.70         0         0
{filter_function:<75}
    1    0    0    0    0    0{visibility:10.5f}     0.000     0.000     0.000{elevation:10.3f}
{sat_height:10.3f}{elevation:10.3f}{sat_view:10.3f}     0.000     0.000     0.000    0          0.000
    2    0{doy:5d}    0
     0.000{sat_view_offset:10.3f}                                                       0.667
     350.0    2600.0       1.0       1.0RN#       NT    T                                                     
    0                                                                                                         
"""

THERMAL_TRANSMITTANCE = """\
T {binary} 7    2    1    0    2    2    2    2    2    2    1    1    0  10.000   0.00
F   0F   0   375.000       1.0    a{ozone:5.3f} 4   T                       
{filter_function:<75}
    1    0    0    0    0    0{visibility:10.3f}     0.000     0.000     0.000{gpheight:10.3f}
{n:5d}    0    0                               {atmospheric_profile}
    50.000 9.510E-01 2.550E+00 7.127E-05 0.000E+00 0.000E+00ABH
    55.000 5.150E-01-3.850E+00 5.991E-05 0.000E+00 0.000E+00ABH
    60.000 2.720E-01-1.605E+01 7.753E-05 0.000E+00 0.000E+00ABH
    70.000 6.700E-02-5.505E+01 7.083E-04 0.000E+00 0.000E+00ABH
    80.000 1.200E-02-9.905E+01 7.237E-02 0.000E+00 0.000E+00ABH
   100.000 0.100E-03-8.265E+01 1.461E-05 0.000E+00 0.000E+00ABH
{sat_height:10.3f}{gpheight:10.3f}{sat_view:10.3f}
    7000.0   14100.0      10.0      20.0TN        NT                           
    3                                                                          
{gpheight:10.3f}{sat_height:10.3f}    55.770                                   
    0                                                                          
"""

SBT_FORMAT = "\n{gpheight:10.3f}{pressure:10.3E}{airtemp:10.3E}{humidity:10.3E}{zero:10.3E}{zero:10.3E}ABH"
