#!/usr/bin/env python
# pragma pylint: disable=line-too-long
"""
These are the templates for the inputs into MODTRAN.
Parameters to be inserted are:
    * albedo
    * water vapour; units: g/cm^2
    * ozone; units: ATM-cm
    * filter function (spectral response of sensor)
    * aerosol_type
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


MIDLAT_SUMMER_ALBEDO = ("""\
TM  2    2    2   -1    2    2    2    2    2    2    1    1    0  10.000{albedo:7.2f}     1
TTT  8   0   375.000  g{water:7.5f}    a{ozone:5.3f} 1   T f f -0.988000t      1.0         0      0.00         0         0
 {filter_function:<75}
    {aerosol_type}    0    0    0    0    0{visibility:10.3f}     0.000     0.000     0.000{elevation:10.3f}
{sat_height:10.3f}{elevation:10.3f}{sat_view:10.3f}     0.000     0.000     0.000    0          0.000
    1    2{doy:5d}    0
{lat:10.5f}{lon:10.5f}   0.00000   0.00000{time:10.5f}{sat_azimuth:10.5f}   0.00000   0.66700
     350.0    2600.0       1.0       1.0RN#       NT       1
    0
""")

TROPICAL_ALBEDO = ("""\
TM  1    2    2   -1    1    1    1    1    1    1    1    1    0  10.000{albedo:7.2f}     1
TTT  8   0   375.000  g{water:7.5f}    a{ozone:5.3f} 1   T f f -0.988000t      1.0         0      0.00         0         0
 {filter_function:<75}
    {aerosol_type}    0    0    0    0    0{visibility:10.3f}     0.000     0.000     0.000{elevation:10.3f}
{sat_height:10.3f}{elevation:10.3f}{sat_view:10.3f}     0.000     0.000     0.000    0          0.000
    1    2{doy:5d}    0
{lat:10.3f}{lon:10.3f}     0.000     0.000{time:10.3f}{sat_azimuth:10.3f}     0.000     0.667
     350.0    2600.0       1.0       1.0RN#       NT       1
    0
""")
