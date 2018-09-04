#!/usr/bin/env python
# pragma pylint: disable=line-too-long

from __future__ import absolute_import, print_function
import pandas as pd

"""
These are the templates for the inputs into MODTRAN.
Parameters to be inserted are:
    * name
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

def midlat_summer_albedo(name, water, ozone, visibility, doy, lats, lons, time, sat_azimuth,\
                         elevation, sat_height, sat_view, albedo, filter_function, binary):
    lat = lats
    lon = lons
    MIDLAT_SUMMER_ALBEDO = {
  "MODTRAN": [
  {
  "MODTRANINPUT":{
    "NAME":name,
    "CASE":0, 
    "RTOPTIONS":{
      "MODTRN":"RT_MODTRAN", 
      "LYMOLC":False, 
      "T_BEST":False, 
      "IEMSCT":"RT_SOLAR_AND_THERMAL", 
      "IMULT":"RT_DISORT",
      "DISALB":True,
      "NSTR":8, 
      "SOLCON":0.0000000000e+00 
    },
    "ATMOSPHERE":{
      "MODEL":"ATM_MIDLAT_SUMMER", 
      "M1":"ATM_MIDLAT_SUMMER", 
      "M2":"ATM_MIDLAT_SUMMER", 
      "M3":"ATM_MIDLAT_SUMMER", 
      "M4":"ATM_MIDLAT_SUMMER", 
      "M5":"ATM_MIDLAT_SUMMER", 
      "M6":"ATM_MIDLAT_SUMMER", 
      "MDEF":1, 
      "CO2MX":3.7500000000e+02, 
      "H2OSTR":water, 
      "H2OUNIT":"G", 
      "O3STR":ozone,
      "O3UNIT":"A", 
      "C_PROF":0, 
      "AERRH":0.0000000000e+00 
    },
    "AEROSOLS":{
      "H2OAER":False, 
      "CDASTM":"t", 
      "ASTMC":3.0000001192e-01, 
      "ASTMX":0.0000000000e+00, 
      "ASTMO":6.9999998808e-01, 
      "APLUS":"  ", 
      "IHAZE":"AER_RURAL", 
      "CNOVAM":False, 
      "ISEASN":"SEASN_AUTO", 
      "ARUSS":"   ", 
      "IVULCN":"STRATO_BACKGROUND", 
      "ICSTL":0, 
      "ICLD":"CLOUD_NONE", 
      "IVSA":False, 
      "VIS":visibility,
      "WSS":0.0000000000e+00, 
      "WHH":0.0000000000e+00, 
      "RAINRT":0.0000000000e+00, 
      "IPH":2, 
      "HGPF":6.6699999571e-01 
    },
    "GEOMETRY":{
      "ITYPE":2, 
      "IPARM":1, 
      "IDAY":doy,
      "PARM1":lat, 
      "PARM2":lon, 
      "PARM3":0.0000000000e+00, 
      "PARM4":0.0000000000e+00, 
      "GMTIME":time, 
      "TRUEAZ":sat_azimuth,
      "ANGLEM":0.0000000000e+00, 
      "H1ALT":sat_height,  
      "H2ALT":elevation,  
      "OBSZEN":sat_view,
      "HRANGE":0.0000000000e+00, 
      "BETA":0.0000000000e+00, 
      "RAD_E":0.0000000000e+00, 
      "LENN":0, 
      "BCKZEN":0.0000000000e+00, 
      "CKRANG":0.0000000000e+00 
    },
    "SURFACE":{
      "SURFTYPE":"REFL_CONSTANT", 
      "TPTEMP":1.0000000000e+01, 
      "SURREF":albedo, 
      "GNDALT":elevation
    },
    "SPECTRAL":{
      "V1":3.5000000000e+02, 
      "V2":2.6000000000e+03, 
      "DV":1.0000000000e+00, 
      "FWHM":1.0000000000e+00, 
      "YFLAG":"R", 
      "XFLAG":"N", 
      "DLIMIT":"#       ", 
      "FLAGS":"NT     ", 
      "MLFLX":1, 
      "VRFRAC":0.0000000000e+00, 
      "SFWHM":0.0000000000e+00, 
      "LSUNFL":"4",
      "LBMNAM":" ", 
      "FILTNM":filter_function
    },
    "FILEOPTIONS":{
      "BINARY":binary,
      "CKPRNT":False,
      "NOPRNT":0 
    }
  }
  }
 ]
}
    
    
    return MIDLAT_SUMMER_ALBEDO



def tropical_albedo(name, water, ozone, visibility, doy, lats, lons, time, sat_azimuth, sat_height, elevation,\
                    sat_view, albedo, filter_function, binary):
    lat = lats
    lon = lons
    TROPICAL_ALBEDO ={
  "MODTRAN": [
  {
  "MODTRANINPUT":{
    "NAME":name, 
    "CASE":0, 
    "RTOPTIONS":{
      "MODTRN":"RT_MODTRAN", 
      "LYMOLC":False, 
      "T_BEST":False, 
      "IEMSCT":"RT_SOLAR_AND_THERMAL", 
      "IMULT":"RT_DISORT",
      "DISALB":True,
      "NSTR":8, 
      "SOLCON":0.0000000000e+00 
    },
    "ATMOSPHERE":{
      "MODEL":"ATM_TROPICAL", 
      "M1":"ATM_TROPICAL", 
      "M2":"ATM_TROPICAL", 
      "M3":"ATM_TROPICAL", 
      "M4":"ATM_TROPICAL", 
      "M5":"ATM_TROPICAL", 
      "M6":"ATM_TROPICAL", 
      "MDEF":1, 
      "CO2MX":3.7500000000e+02, 
      "H2OSTR":water, 
      "H2OUNIT":"G", 
      "O3STR":ozone, 
      "O3UNIT":"A", 
      "C_PROF":0, 
      "AERRH":0.0000000000e+00 
    },
    "AEROSOLS":{
      "H2OAER":False, 
      "CDASTM":"t", 
      "ASTMC":3.0000001192e-01, 
      "ASTMX":0.0000000000e+00, 
      "ASTMO":6.9999998808e-01, 
      "APLUS":"  ", 
      "IHAZE":"AER_RURAL", 
      "CNOVAM":False, 
      "ISEASN":"SEASN_AUTO", 
      "ARUSS":"   ", 
      "IVULCN":"STRATO_BACKGROUND", 
      "ICSTL":0, 
      "ICLD":"CLOUD_NONE", 
      "IVSA":False, 
      "VIS":visibility, 
      "WSS":0.0000000000e+00, 
      "WHH":0.0000000000e+00, 
      "RAINRT":0.0000000000e+00, 
      "IPH":2, 
      "HGPF":6.6699999571e-01 
    },
    "GEOMETRY":{
      "ITYPE":2, 
      "IPARM":1, 
      "IDAY":doy, 
      "PARM1":lat, 
      "PARM2":lon, 
      "PARM3":0.0000000000e+00, 
      "PARM4":0.0000000000e+00, 
      "GMTIME":time, 
      "TRUEAZ":sat_azimuth, 
      "ANGLEM":0.0000000000e+00, 
      "H1ALT":sat_height, 
      "H2ALT":elevation, 
      "OBSZEN":sat_view, 
      "HRANGE":0.0000000000e+00, 
      "BETA":0.0000000000e+00, 
      "RAD_E":0.0000000000e+00, 
      "LENN":0, 
      "BCKZEN":0.0000000000e+00, 
      "CKRANG":0.0000000000e+00 
    },
    "SURFACE":{
      "SURFTYPE":"REFL_CONSTANT", 
      "TPTEMP":1.0000000000e+01, 
      "SURREF":albedo, 
      "GNDALT":elevation 
    },
    "SPECTRAL":{
      "V1":3.5000000000e+02, 
      "V2":2.6000000000e+03, 
      "DV":1.0000000000e+00, 
      "FWHM":1.0000000000e+00, 
      "YFLAG":"R", 
      "XFLAG":"N", 
      "DLIMIT":"#       ", 
      "FLAGS":"NT     ", 
      "MLFLX":1, 
      "VRFRAC":0.0000000000e+00, 
      "SFWHM":0.0000000000e+00, 
      "LSUNFL":"4",
      "LBMNAM":" ", 
      "FILTNM":filter_function
    },
    "FILEOPTIONS":{
      "BINARY":binary, 
      "CKPRNT":False,
      "NOPRNT":0 
    }
  }
  }
 ]
} 

    return TROPICAL_ALBEDO


def thermal_transmittance(name, ozone, n, prof_alt, prof_pres, prof_temp, prof_water, visibility, sat_height, gpheight,\
                          sat_view, filter_function, binary):
    
    prof_altitude = prof_alt  + [5.0000000000e+01, 5.5000000000e+01, 6.0000000000e+01, 7.0000000000e+01, 8.0000000000e+01, 1.0000000000e+02 ]
    prof_pressure = prof_pres + [9.5099997520e-01, 5.1499998569e-01, 2.7200001478e-01, 6.7000001669e-02, 1.2000000104e-02, 9.9999997474e-05 ]
    prof_temperature = prof_temp + [2.5499999523e+00, -3.8499999046e+00, -1.6049999237e+01, -5.5049999237e+01, -9.9050003052e+01 -8.2650001526e+01 ]
    prof_h20 = prof_water + [7.1269998443e-05, 5.9909998527e-05, 7.7529999544e-05, 7.0829998003e-04, 7.2370000184e-02, 1.4610000107e-05 ]
    
    
    THERMAL_TRANSMITTANCE = {
  "MODTRAN": [
  {
  "MODTRANINPUT":{
    "NAME":name, 
    "CASE":0, 
    "RTOPTIONS":{
      "MODTRN":"RT_MODTRAN", 
      "LYMOLC":False, 
      "T_BEST":False, 
      "IEMSCT":"RT_THERMAL_ONLY", 
      "IMULT":"RT_NO_MULTIPLE_SCATTER", 
      "DISALB":False, 
      "NSTR":0, 
      "SOLCON":0.0000000000e+00 
    },
    "ATMOSPHERE":{
      "MODEL":"ATM_USER_ALT_PROFILE", 
      "M1":"ATM_MIDLAT_SUMMER", 
      "M2":"ATM_MIDLAT_SUMMER", 
      "M3":"ATM_MIDLAT_SUMMER", 
      "M4":"ATM_MIDLAT_SUMMER", 
      "M5":"ATM_MIDLAT_SUMMER", 
      "M6":"ATM_MIDLAT_SUMMER", 
      "MDEF":1, 
      "CO2MX":3.7500000000e+02, 
      "H2OSTR":1.0000000000e+00, 
      "O3STR":ozone, 
      "O3UNIT":"A", 
      "C_PROF":0, 
      "AERRH":0.0000000000e+00, 
      "AYRANG":False, 
      "E_MASS":0.0000000000e+00, 
      "AIRMWT":0.0000000000e+00, 
      "NLAYERS":n, 
      "NPROF":4, 
      "PROFILES":[
           {
                "TYPE":"PROF_ALTITUDE",
                "UNITS":"UNT_KILOMETERS",
                "PROFILE": prof_altitude
           },
           {
                "TYPE":"PROF_PRESSURE",
                "UNITS":"UNT_PMILLIBAR",
                "PROFILE":prof_pressure
           },
           {
                "TYPE":"PROF_TEMPERATURE",
                "UNITS":"UNT_TCELSIUS",
                "PROFILE":prof_temperature
           },
           {
                "TYPE":"PROF_H2O",
                "UNITS":"UNT_REL_HUMIDITY",
                "PROFILE":prof_h20
           }
           ]
    },
    "AEROSOLS":{
      "H2OAER":False, 
      "CDASTM":" ", 
      "ASTMC":0.0000000000e+00, 
      "ASTMX":0.0000000000e+00, 
      "ASTMO":0.0000000000e+00, 
      "APLUS":"  ", 
      "IHAZE":"AER_RURAL", 
      "CNOVAM":False, 
      "ISEASN":"SEASN_AUTO", 
      "ARUSS":"   ", 
      "IVULCN":"STRATO_BACKGROUND", 
      "ICSTL":0, 
      "ICLD":"CLOUD_NONE", 
      "IVSA":False, 
      "VIS":visibility, 
      "WSS":0.0000000000e+00, 
      "WHH":0.0000000000e+00, 
      "RAINRT":0.0000000000e+00 
    },
    "GEOMETRY":{
      "ITYPE":2, 
      "H1ALT":sat_height, 
      "H2ALT":gpheight, 
      "OBSZEN":sat_view, 
      "HRANGE":0.0000000000e+00, 
      "BETA":0.0000000000e+00, 
      "RAD_E":0.0000000000e+00, 
      "LENN":0, 
      "BCKZEN":0.0000000000e+00, 
      "CKRANG":0.0000000000e+00 
    },
    "SURFACE":{
      "SURFTYPE":"REFL_CONSTANT", 
      "TPTEMP":1.0000000000e+01, 
      "SURREF":0.0000000000e+00, 
      "GNDALT":gpheight 
    },
    "SPECTRAL":{
      "V1":7.0000000000e+03, 
      "V2":1.4100000000e+04, 
      "DV":1.0000000000e+01, 
      "FWHM":2.0000000000e+01, 
      "YFLAG":"T", 
      "XFLAG":"N", 
      "FLAGS":"NT     ", 
      "MLFLX":0, 
      "VRFRAC":0.0000000000e+00, 
      "SFWHM":0.0000000000e+00, 
      "LSUNFL":"4", 
      "LBMNAM":" ", 
      "FILTNM":filter_function 
    },
    "FILEOPTIONS":{

      "BINARY":binary,
      "CKPRNT":False,
      "NOPRNT":0 
    }
  }
  },
  {
  "MODTRANINPUT":{
    "NAME":name, 
    "CASE":1, 
    "RTOPTIONS":{
      "MODTRN":"RT_MODTRAN", 
      "LYMOLC":False, 
      "T_BEST":False, 
      "IEMSCT":"RT_THERMAL_ONLY", 
      "IMULT":"RT_NO_MULTIPLE_SCATTER", 
      "DISALB":False, 
      "NSTR":0, 
      "SOLCON":0.0000000000e+00 
    },
    "ATMOSPHERE":{
      "MODEL":"ATM_USER_ALT_PROFILE", 
      "M1":"ATM_MIDLAT_SUMMER", 
      "M2":"ATM_MIDLAT_SUMMER", 
      "M3":"ATM_MIDLAT_SUMMER", 
      "M4":"ATM_MIDLAT_SUMMER", 
      "M5":"ATM_MIDLAT_SUMMER", 
      "M6":"ATM_MIDLAT_SUMMER", 
      "MDEF":1, 
      "CO2MX":3.7500000000e+02, 
      "H2OSTR":1.0000000000e+00, 
      "O3STR":ozone, 
      "O3UNIT":"A", 
      "C_PROF":0, 
      "AERRH":0.0000000000e+00, 
      "AYRANG":False, 
      "E_MASS":0.0000000000e+00, 
      "AIRMWT":0.0000000000e+00, 
      "NLAYERS":n, 
      "NPROF":4, 
      "PROFILES":[
           {
                "TYPE":"PROF_ALTITUDE",
                "UNITS":"UNT_KILOMETERS",
                "PROFILE":prof_altitude
           },
           {
                "TYPE":"PROF_PRESSURE",
                "UNITS":"UNT_PMILLIBAR",
                "PROFILE":prof_pressure
           },
           {
                "TYPE":"PROF_TEMPERATURE",
                "UNITS":"UNT_TCELSIUS",
                "PROFILE":prof_temperature
           },
           {
                "TYPE":"PROF_H2O",
                "UNITS":"UNT_REL_HUMIDITY",
                "PROFILE":prof_h20
           }
           ]
    },
    "AEROSOLS":{
      "H2OAER":False, 
      "CDASTM":" ", 
      "ASTMC":0.0000000000e+00, 
      "ASTMX":0.0000000000e+00, 
      "ASTMO":0.0000000000e+00, 
      "APLUS":"  ", 
      "IHAZE":"AER_RURAL", 
      "CNOVAM":False, 
      "ISEASN":"SEASN_AUTO", 
      "ARUSS":"   ", 
      "IVULCN":"STRATO_BACKGROUND", 
      "ICSTL":0, 
      "ICLD":"CLOUD_NONE", 
      "IVSA":False, 
      "VIS":visibility, 
      "WSS":0.0000000000e+00, 
      "WHH":0.0000000000e+00, 
      "RAINRT":0.0000000000e+00 
    },
    "GEOMETRY":{
      "ITYPE":2, 
      "H1ALT":gpheight, 
      "H2ALT":sat_height, 
      "OBSZEN":5.5770000000e+01, 
      "HRANGE":0.0000000000e+00, 
      "BETA":0.0000000000e+00, 
      "RAD_E":0.0000000000e+00, 
      "LENN":0, 
      "BCKZEN":0.0000000000e+00, 
      "CKRANG":0.0000000000e+00 
    },
    "SURFACE":{
      "SURFTYPE":"REFL_CONSTANT", 
      "TPTEMP":1.0000000000e+01, 
      "SURREF":0.0000000000e+00, 
      "GNDALT":gpheight 
    },
    "SPECTRAL":{
      "V1":7.0000000000e+03, 
      "V2":1.4100000000e+04, 
      "DV":1.0000000000e+01, 
      "FWHM":2.0000000000e+01, 
      "YFLAG":"T", 
      "XFLAG":"N", 
      "FLAGS":"NT     ", 
      "MLFLX":0, 
      "VRFRAC":0.0000000000e+00, 
      "SFWHM":0.0000000000e+00, 
      "LSUNFL":"4", 
      "LBMNAM":" ", 
      "FILTNM":filter_function
    },
    "FILEOPTIONS":{
      "BINARY":binary, 
      "CKPRNT":False,
      "NOPRNT":0 
    }
  }
  }
 ]
}
    return THERMAL_TRANSMITTANCE





























    
