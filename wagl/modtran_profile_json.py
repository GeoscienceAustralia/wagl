#!/usr/bin/env python
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
# pragma pylint: disable=line-too-long

from __future__ import absolute_import, print_function


def midlat_summer_albedo(
    name,
    water,
    ozone,
    visibility,
    doy,
    lat,
    lon,
    time,
    sat_azimuth,
    elevation,
    sat_height,
    sat_view,
    albedo,
    filter_function,
    binary,
):
    """
    MODTRAN 6.0.1 input: 'json' format template for mid latitude summer albedo
    """

    _midlat_summer_albedo = {
        "MODTRAN": [
            {
                "MODTRANINPUT": {
                    "NAME": name,
                    "CASE": 0,
                    "RTOPTIONS": {
                        "MODTRN": "RT_MODTRAN",
                        "LYMOLC": False,
                        "T_BEST": False,
                        "IEMSCT": "RT_SOLAR_AND_THERMAL",
                        "IMULT": "RT_DISORT",
                        "DISALB": True,
                        "NSTR": 8,
                        "SOLCON": -0.98799997568,
                    },
                    "ATMOSPHERE": {
                        "MODEL": "ATM_MIDLAT_SUMMER",
                        "M1": "ATM_MIDLAT_SUMMER",
                        "M2": "ATM_MIDLAT_SUMMER",
                        "M3": "ATM_MIDLAT_SUMMER",
                        "M4": "ATM_MIDLAT_SUMMER",
                        "M5": "ATM_MIDLAT_SUMMER",
                        "M6": "ATM_MIDLAT_SUMMER",
                        "MDEF": 1,
                        "CO2MX": 375.0,
                        "H2OSTR": water,
                        "H2OUNIT": "G",
                        "O3STR": ozone,
                        "O3UNIT": "A",
                        "C_PROF": 0,
                        "AERRH": 0.0,
                    },
                    "AEROSOLS": {
                        "H2OAER": False,
                        "CDASTM": "t",
                        "ASTMC": 0.30000001192,
                        "ASTMX": 0.0,
                        "ASTMO": 0.69999998808,
                        "APLUS": "  ",
                        "IHAZE": "AER_RURAL",
                        "CNOVAM": False,
                        "ISEASN": "SEASN_AUTO",
                        "ARUSS": "   ",
                        "IVULCN": "STRATO_BACKGROUND",
                        "ICSTL": 0,
                        "ICLD": "CLOUD_NONE",
                        "IVSA": False,
                        "VIS": visibility,
                        "WSS": 0.0,
                        "WHH": 0.0,
                        "RAINRT": 0.0,
                        "IPH": 2,
                        "HGPF": 0.66699999571,
                    },
                    "GEOMETRY": {
                        "ITYPE": 2,
                        "IPARM": 1,
                        "IDAY": doy,
                        "PARM1": lat,
                        "PARM2": lon,
                        "PARM3": 0.0,
                        "PARM4": 0.0,
                        "GMTIME": time,
                        "TRUEAZ": sat_azimuth,
                        "ANGLEM": 0.0,
                        "H1ALT": sat_height,
                        "H2ALT": elevation,
                        "OBSZEN": sat_view,
                        "HRANGE": 0.0,
                        "BETA": 0.0,
                        "RAD_E": 0.0,
                        "LENN": 0,
                        "BCKZEN": 0.0,
                        "CKRANG": 0.0,
                    },
                    "SURFACE": {
                        "SURFTYPE": "REFL_CONSTANT",
                        "TPTEMP": 10.0,
                        "SURREF": albedo,
                        "GNDALT": elevation,
                    },
                    "SPECTRAL": {
                        "V1": 350.0,
                        "V2": 2600.0,
                        "DV": 1.0,
                        "FWHM": 1.0,
                        "YFLAG": "R",
                        "XFLAG": "N",
                        "DLIMIT": "#       ",
                        "FLAGS": "NT     ",
                        "MLFLX": 1,
                        "VRFRAC": 0.0,
                        "SFWHM": 0.0,
                        "LSUNFL": "1",
                        "LBMNAM": " ",
                        "FILTNM": filter_function,
                    },
                    "FILEOPTIONS": {"BINARY": binary, "CKPRNT": False, "NOPRNT": 0},
                }
            }
        ]
    }
    return _midlat_summer_albedo


def tropical_albedo(
    name,
    water,
    ozone,
    visibility,
    doy,
    lat,
    lon,
    time,
    sat_azimuth,
    sat_height,
    elevation,
    sat_view,
    albedo,
    filter_function,
    binary,
):
    """
    MODTRAN 6.0.1 input: 'json' format template for tropical albedo
    """

    _tropical_albedo = {
        "MODTRAN": [
            {
                "MODTRANINPUT": {
                    "NAME": name,
                    "CASE": 0,
                    "RTOPTIONS": {
                        "MODTRN": "RT_MODTRAN",
                        "LYMOLC": False,
                        "T_BEST": False,
                        "IEMSCT": "RT_SOLAR_AND_THERMAL",
                        "IMULT": "RT_DISORT",
                        "DISALB": True,
                        "NSTR": 8,
                        "SOLCON": -0.98799997568,
                    },
                    "ATMOSPHERE": {
                        "MODEL": "ATM_TROPICAL",
                        "M1": "ATM_TROPICAL",
                        "M2": "ATM_TROPICAL",
                        "M3": "ATM_TROPICAL",
                        "M4": "ATM_TROPICAL",
                        "M5": "ATM_TROPICAL",
                        "M6": "ATM_TROPICAL",
                        "MDEF": 1,
                        "CO2MX": 375.0,
                        "H2OSTR": water,
                        "H2OUNIT": "G",
                        "O3STR": ozone,
                        "O3UNIT": "A",
                        "C_PROF": 0,
                        "AERRH": 0.0,
                    },
                    "AEROSOLS": {
                        "H2OAER": False,
                        "CDASTM": "t",
                        "ASTMC": 0.30000001192,
                        "ASTMX": 0.0,
                        "ASTMO": 0.69999998808,
                        "APLUS": "  ",
                        "IHAZE": "AER_RURAL",
                        "CNOVAM": False,
                        "ISEASN": "SEASN_AUTO",
                        "ARUSS": "   ",
                        "IVULCN": "STRATO_BACKGROUND",
                        "ICSTL": 0,
                        "ICLD": "CLOUD_NONE",
                        "IVSA": False,
                        "VIS": visibility,
                        "WSS": 0.0,
                        "WHH": 0.0,
                        "RAINRT": 0.0,
                        "IPH": 2,
                        "HGPF": 0.66699999571,
                    },
                    "GEOMETRY": {
                        "ITYPE": 2,
                        "IPARM": 1,
                        "IDAY": doy,
                        "PARM1": lat,
                        "PARM2": lon,
                        "PARM3": 0.0,
                        "PARM4": 0.0,
                        "GMTIME": time,
                        "TRUEAZ": sat_azimuth,
                        "ANGLEM": 0.0,
                        "H1ALT": sat_height,
                        "H2ALT": elevation,
                        "OBSZEN": sat_view,
                        "HRANGE": 0.0,
                        "BETA": 0.0,
                        "RAD_E": 0.0,
                        "LENN": 0,
                        "BCKZEN": 0.0,
                        "CKRANG": 0.0,
                    },
                    "SURFACE": {
                        "SURFTYPE": "REFL_CONSTANT",
                        "TPTEMP": 10.0,
                        "SURREF": albedo,
                        "GNDALT": elevation,
                    },
                    "SPECTRAL": {
                        "V1": 350.0,
                        "V2": 2600.0,
                        "DV": 1.0,
                        "FWHM": 1.0,
                        "YFLAG": "R",
                        "XFLAG": "N",
                        "DLIMIT": "#       ",
                        "FLAGS": "NT     ",
                        "MLFLX": 1,
                        "VRFRAC": 0.0,
                        "SFWHM": 0.0,
                        "LSUNFL": "1",
                        "LBMNAM": " ",
                        "FILTNM": filter_function,
                    },
                    "FILEOPTIONS": {"BINARY": binary, "CKPRNT": False, "NOPRNT": 0},
                }
            }
        ]
    }

    return _tropical_albedo


def thermal_transmittance(
    name,
    ozone,
    n,
    prof_alt,
    prof_pres,
    prof_temp,
    prof_water,
    visibility,
    sat_height,
    gpheight,
    sat_view,
    filter_function,
    binary,
):
    """
    MODTRAN 6.0.1 input: 'json' format template for thermal transmittance
    """

    prof_altitude = prof_alt + [
        5.0000000000e01,
        5.5000000000e01,
        6.0000000000e01,
        7.0000000000e01,
        8.0000000000e01,
        1.0000000000e02,
    ]

    prof_pressure = prof_pres + [
        9.5099997520e-01,
        5.1499998569e-01,
        2.7200001478e-01,
        6.7000001669e-02,
        1.2000000104e-02,
        9.9999997474e-05,
    ]

    prof_temperature = prof_temp + [
        2.5499999523e00,
        -3.8499999046e00,
        -1.6049999237e01,
        -5.5049999237e01,
        -9.9050003052e01,
        -8.2650001526e01,
    ]

    prof_h20 = prof_water + [
        7.1269998443e-05,
        5.9909998527e-05,
        7.7529999544e-05,
        7.0829998003e-04,
        7.2370000184e-02,
        1.4610000107e-05,
    ]

    _thermal_transmittance = {
        "MODTRAN": [
            {
                "MODTRANINPUT": {
                    "NAME": name,
                    "CASE": 0,
                    "RTOPTIONS": {
                        "MODTRN": "RT_MODTRAN",
                        "LYMOLC": False,
                        "T_BEST": False,
                        "IEMSCT": "RT_THERMAL_ONLY",
                        "IMULT": "RT_NO_MULTIPLE_SCATTER",
                        "DISALB": False,
                        "NSTR": 0,
                        "SOLCON": -0.98799997568,
                    },
                    "ATMOSPHERE": {
                        "MODEL": "ATM_USER_ALT_PROFILE",
                        "M1": "ATM_MIDLAT_SUMMER",
                        "M2": "ATM_MIDLAT_SUMMER",
                        "M3": "ATM_MIDLAT_SUMMER",
                        "M4": "ATM_MIDLAT_SUMMER",
                        "M5": "ATM_MIDLAT_SUMMER",
                        "M6": "ATM_MIDLAT_SUMMER",
                        "MDEF": 1,
                        "CO2MX": 375.0,
                        "H2OSTR": 1.0,
                        "O3STR": ozone,
                        "O3UNIT": "A",
                        "C_PROF": 0,
                        "AERRH": 0.0,
                        "AYRANG": False,
                        "E_MASS": 0.0,
                        "AIRMWT": 0.0,
                        "NLAYERS": n,
                        "NPROF": 4,
                        "PROFILES": [
                            {
                                "TYPE": "PROF_ALTITUDE",
                                "UNITS": "UNT_KILOMETERS",
                                "PROFILE": prof_altitude,
                            },
                            {
                                "TYPE": "PROF_PRESSURE",
                                "UNITS": "UNT_PMILLIBAR",
                                "PROFILE": prof_pressure,
                            },
                            {
                                "TYPE": "PROF_TEMPERATURE",
                                "UNITS": "UNT_TCELSIUS",
                                "PROFILE": prof_temperature,
                            },
                            {
                                "TYPE": "PROF_H2O",
                                "UNITS": "UNT_REL_HUMIDITY",
                                "PROFILE": prof_h20,
                            },
                        ],
                    },
                    "AEROSOLS": {
                        "H2OAER": False,
                        "CDASTM": " ",
                        "ASTMC": 0.0,
                        "ASTMX": 0.0,
                        "ASTMO": 0.0,
                        "APLUS": "  ",
                        "IHAZE": "AER_RURAL",
                        "CNOVAM": False,
                        "ISEASN": "SEASN_AUTO",
                        "ARUSS": "   ",
                        "IVULCN": "STRATO_BACKGROUND",
                        "ICSTL": 0,
                        "ICLD": "CLOUD_NONE",
                        "IVSA": False,
                        "VIS": visibility,
                        "WSS": 0.0,
                        "WHH": 0.0,
                        "RAINRT": 0.0,
                    },
                    "GEOMETRY": {
                        "ITYPE": 2,
                        "H1ALT": sat_height,
                        "H2ALT": gpheight,
                        "OBSZEN": sat_view,
                        "HRANGE": 0.0,
                        "BETA": 0.0,
                        "RAD_E": 0.0,
                        "LENN": 0,
                        "BCKZEN": 0.0,
                        "CKRANG": 0.0,
                    },
                    "SURFACE": {
                        "SURFTYPE": "REFL_CONSTANT",
                        "TPTEMP": 10.0,
                        "SURREF": 0.0,
                        "GNDALT": gpheight,
                    },
                    "SPECTRAL": {
                        "V1": 7000.0,
                        "V2": 14100.0,
                        "DV": 10.0,
                        "FWHM": 20.0,
                        "YFLAG": "T",
                        "XFLAG": "N",
                        "FLAGS": "NT     ",
                        "MLFLX": 0,
                        "VRFRAC": 0.0,
                        "SFWHM": 0.0,
                        "LSUNFL": "1",
                        "LBMNAM": " ",
                        "FILTNM": filter_function,
                    },
                    "FILEOPTIONS": {"BINARY": binary, "CKPRNT": False, "NOPRNT": 0},
                }
            },
            {
                "MODTRANINPUT": {
                    "NAME": name,
                    "CASE": 1,
                    "RTOPTIONS": {
                        "MODTRN": "RT_MODTRAN",
                        "LYMOLC": False,
                        "T_BEST": False,
                        "IEMSCT": "RT_THERMAL_ONLY",
                        "IMULT": "RT_NO_MULTIPLE_SCATTER",
                        "DISALB": False,
                        "NSTR": 0,
                        "SOLCON": -0.98799997568,
                    },
                    "ATMOSPHERE": {
                        "MODEL": "ATM_USER_ALT_PROFILE",
                        "M1": "ATM_MIDLAT_SUMMER",
                        "M2": "ATM_MIDLAT_SUMMER",
                        "M3": "ATM_MIDLAT_SUMMER",
                        "M4": "ATM_MIDLAT_SUMMER",
                        "M5": "ATM_MIDLAT_SUMMER",
                        "M6": "ATM_MIDLAT_SUMMER",
                        "MDEF": 1,
                        "CO2MX": 375.0,
                        "H2OSTR": 1.0,
                        "O3STR": ozone,
                        "O3UNIT": "A",
                        "C_PROF": 0,
                        "AERRH": 0.0,
                        "AYRANG": False,
                        "E_MASS": 0.0,
                        "AIRMWT": 0.0,
                        "NLAYERS": n,
                        "NPROF": 4,
                        "PROFILES": [
                            {
                                "TYPE": "PROF_ALTITUDE",
                                "UNITS": "UNT_KILOMETERS",
                                "PROFILE": prof_altitude,
                            },
                            {
                                "TYPE": "PROF_PRESSURE",
                                "UNITS": "UNT_PMILLIBAR",
                                "PROFILE": prof_pressure,
                            },
                            {
                                "TYPE": "PROF_TEMPERATURE",
                                "UNITS": "UNT_TCELSIUS",
                                "PROFILE": prof_temperature,
                            },
                            {
                                "TYPE": "PROF_H2O",
                                "UNITS": "UNT_REL_HUMIDITY",
                                "PROFILE": prof_h20,
                            },
                        ],
                    },
                    "AEROSOLS": {
                        "H2OAER": False,
                        "CDASTM": " ",
                        "ASTMC": 0.0,
                        "ASTMX": 0.0,
                        "ASTMO": 0.0,
                        "APLUS": "  ",
                        "IHAZE": "AER_RURAL",
                        "CNOVAM": False,
                        "ISEASN": "SEASN_AUTO",
                        "ARUSS": "   ",
                        "IVULCN": "STRATO_BACKGROUND",
                        "ICSTL": 0,
                        "ICLD": "CLOUD_NONE",
                        "IVSA": False,
                        "VIS": visibility,
                        "WSS": 0.0,
                        "WHH": 0.0,
                        "RAINRT": 0.0,
                    },
                    "GEOMETRY": {
                        "ITYPE": 2,
                        "H1ALT": gpheight,
                        "H2ALT": sat_height,
                        "OBSZEN": 55.77,
                        "HRANGE": 0.0,
                        "BETA": 0.0,
                        "RAD_E": 0.0,
                        "LENN": 0,
                        "BCKZEN": 0.0,
                        "CKRANG": 0.0,
                    },
                    "SURFACE": {
                        "SURFTYPE": "REFL_CONSTANT",
                        "TPTEMP": 10.0,
                        "SURREF": 0.0,
                        "GNDALT": gpheight,
                    },
                    "SPECTRAL": {
                        "V1": 7000.0,
                        "V2": 14100.0,
                        "DV": 10.0,
                        "FWHM": 20.0,
                        "YFLAG": "T",
                        "XFLAG": "N",
                        "FLAGS": "NT     ",
                        "MLFLX": 0,
                        "VRFRAC": 0.0,
                        "SFWHM": 0.0,
                        "LSUNFL": "1",
                        "LBMNAM": " ",
                        "FILTNM": filter_function,
                    },
                    "FILEOPTIONS": {"BINARY": binary, "CKPRNT": False, "NOPRNT": 0},
                }
            },
        ]
    }

    return _thermal_transmittance
