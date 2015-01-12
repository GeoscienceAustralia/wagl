from enum import Enum, unique 

class BandType(Enum):
    REF = 'Reflective'
    THM = 'Thermal'
    PAN = 'Panchromatic'
    ATM = 'Atmosphere'
    BQA = 'Quality'

class Band(object):
    """
    A Band represents a region of the spectrum in which a specific
    Sensor is sensitive.

    A Sensor may have several Bands or regions of spectral sensitivity.
    """

    def __init__(self, id="1", desc="Band 1", resolution_M=25.0, \
        band_type=BandType.REF, wavelength_nm=[0.5,0.6], \
        other_props={}):

        assert(isinstance(desc, str))
        assert(isinstance(resolution_M, float))
        assert(isinstance(band_type, BandType))
        assert(isinstance(wavelength_nm, list))
#        assert(len(wavelength_nm)==2)

        self.id = id
        self.desc = desc
        self.resolution_M = resolution_M
        self.band_type = band_type
        self.wavelength_nm = wavelength_nm
        self.__dict__.update(other_props)

    def __str__(self):
        return "%s, %s (%s) %s" % (self.id, self.desc, \
            self.band_type.name, str(self.wavelength_nm))

    def __repr__(self):
        return "Band: %s" % (str(self), )


class Sensor(object):
    """
    A Sensor is an instrument package, carried by a SensorPlatform
    (e.g. a Satellite or Spacecraft), capable of acquiring irradiance
    data in one or Bands.

    """ 
    def __init__(self, id="ID", desc="A sensor", bands=[], \
        other_props={}):

        self.id = id
        self.desc = desc
        self.bands = bands
        self.__dict__.update(other_props)


#        self.K1 = K1
#        self.K2 = K2
#        self.rgb_bands = rgb_bands
#        self.root_band = root_band
#        self.acquisition_seconds = acquisition_seconds

    def __str__(self):
        return "%s (%s)" % (self.id, self.desc)

    def __repr__(self):
        return "Sensor: %s" % (str(self))

@unique
class Spacecraft(Enum):
    """
    A spacecraft is a SensorPlatform which flies in space (typically in 
    earth orbit).

    Within GA's NEMO division we are concerned with a limited number of 
    Spacecraft. We have chosen to define these as a python Enumeration.
    As new Spacecraft are launched and commissioned, this Enumeration 
    will be updated with the new spacecraft's information.

    TODO: Make enumeration members immutable
    """

    __order__ = 'LS5 LS7 LS8'

    LS5 = ( [  # sensors \
            Sensor("MSS", "Multi-Spectral Scanner", bands=[ \
                Band("1", "Visible Green", 60.0, BandType.REF, [0.6,0.7]), \
                Band("2", "Visible Red", 60.0, BandType.REF, [0.6, 0.7]), \
                Band("3", "Near Infrared 1", 60.0, BandType.REF, [0.7, 0.8]), \
                Band("4", "Near Infrared 2", 60.0, BandType.REF, [0.8, 1.1]) \
            ], other_props={ \
                "K1": 0,
                "K2": 0,
                "rgb_bands": "",
                "root_band": "10",
                "acquisition_seconds": 0}),
            Sensor("TM", "Thematic Mapper", bands=[ \
                Band("1", "Visible Blue", 25.0, BandType.REF, [0.45, 0.52]), \
                Band("2", "Visible Green", 25.0, BandType.REF, [0.52, 0.60]), \
                Band("3", "Visible Red", 25.0, BandType.REF, [0.63, 0.69]), \
                Band("4", "Near Infrared", 25.0, BandType.REF, [0.76, 0.90]), \
                Band("5", "Middle Infrared 1", 25.0, BandType.REF, [1.55, 1.75]), \
                Band("6", "Thermal Infrared", 100.0, BandType.THM, [10.4, 12.5]), \
                Band("7", "Middle Infrared 2", 25.0, BandType.REF, [2.08, 2.35]) \
            ], other_props={ \
                "K1": 607.76,
                "K2": 1260.56,
                "acquisition_seconds": 24}) \
        ], {  # other spacecraft properties \
            "altitude": 705000.0, \
            "classification": "U", \
            "format": "EOSAT Fast Format", \
            "inclination": 1.7139133254584316445390643346558, \
            "intl_designator": "84021A", \
            "name_pattern": "^L.*5$", \
            "satellite_name": "Landsat-5", \
            "nominal_pixel_degrees": 0.000250, \
            "number": 14780, \
            "omega": 0.001059, \
            "projection": "EQUIRECT_0_0", \
            "radius": 7285600.0, \
            "semi_major_axis": 7083160.0, \
            "solar_irrad_file": "solar_irrad_landsat5.txt", \
            "spectral_filter_file": "landsat5_vsir.flt", \
            "sweep_period": 0.07342143906020558002936857562408, \
            "tag": "LS5", \
            "tle_format": "l5_%4d%s_norad.txt" } \
        ) 

    LS7 = ( [  # sensors \
            Sensor("ETM+", "Enhanced Thematic Mapper Plus", bands=[ \
                Band("1", "Visible Blue", 25.0, BandType.REF, [0.45, 0.52]), \
                Band("2", "Visible Green", 25.0, BandType.REF, [0.52, 0.60]), \
                Band("3", "Visible Red", 25.0, BandType.REF, [0.63, 0.69]), \
                Band("4", "Near Infrared", 25.0, BandType.REF, [0.76, 0.90]), \
                Band("5", "Middle Infrared 1", 25.0, BandType.REF, [1.55, 1.75]), \
                Band("61", "Thermal Infrared (Low Gain)", 50.0, BandType.THM, [10.4, 12.5]), \
                Band("62", "Thermal Infrared (High Gaine)", 50.0, BandType.THM, [10.4, 12.5]), \
                Band("7", "Middle Infrared 2", 25.0, BandType.REF, [2.08, 2.35]), \
                Band("8", "Panchromatic", 12.5, BandType.PAN, [0.52, 0.90]), \
            ], other_props={ \
                "K1": 666.09,
                "K2": 1282.71,
                "acquisition_seconds": 24}) \
        ], {  # other spacecraft properties \
            "altitude": 705000.0, \
            "classification": "U", \
            "format": "FastL7A", \
            "inclination": 1.7139133254584316445390643346558, \
            "intl_designator": "99020A", \
            "satellite_name": "Landsat-7", \
            "name_pattern": "^L.*7$", \
            "nominal_pixel_degrees": 0.000250, \
            "number": 25682, \
            "omega": 0.001059, \
            "projection": "EQUIRECT_0_0", \
            "radius": 7285600.0, \
            "semi_major_axis": 7083160.0, \
            "solar_irrad_file": "solar_irrad_landsat7.txt", \
            "spectral_filter_file": "landsat7_vsir.flt", \
            "sweep_period": 0.07342143906020558002936857562408, \
            "tag": "LS7", \
            "tle_format": "L7%4d%sASNNOR.S00" } \
        ) 

    LS8 = ( [  # sensors \
            Sensor("OLI", "Landsat Data Continuity Missions", bands=[ \
                Band("1", "Coastal Aerosol", 25.0, BandType.REF, [0.433,0.453]), \
                Band("2", "Visible Blue", 25.0, BandType.REF, [0.45, 0.515]), \
                Band("3", "Visible Green", 25.0, BandType.REF, [0.525, 0.60]), \
                Band("4", "Visible Red", 25.0, BandType.REF, [0.63, 0.68]), \
                Band("5", "Near Infrared", 25.0, BandType.REF, [0.845, 0.885]), \
                Band("6", "Short-wave Infrared 1", 25.0, BandType.REF, [1.56, 1.66]), \
                Band("7", "Short-wave Infrared 2", 25.0, BandType.REF, [2.1, 2.3]), \
                Band("8", "Panchromatic", 12.5, BandType.PAN, [0.5, 0.68]), \
                Band("9", "Cirrus", 25.0, BandType.ATM, [1.360, 1.390]), \
                Band("quality", "Quality", 25.0, BandType.BQA, []) \
            ], other_props={ \
                "K1": 0,
                "K2": 0,
                "acquisition_seconds": 24 }),
            Sensor("TIRS", "Thermal InfraRed Sensor", bands=[ \
                Band("10", "Thermal Infrared 1", 25.0, BandType.THM, [10.3, 11.3]), \
                Band("11", "Thermal Infrared 2", 25.0, BandType.THM, [11.5, 12.5]) \
            ], other_props={ \
                "K1": 774.89,
                "K2": 1321.08,
                "acquisition_seconds": 24}), \

            # the following is a hack to support a previous (perhaps questionable)
            # decision to invent the "OLI_TIRS" sensor.
            # "OLI_TIRS" is not a sensor, it is a combination 
            # of two sensors. TODO: Fix this hack througout Image Processor
            # and remove the following OLI_TIRS sensor definition

            Sensor("OLI_TIRS", "Thermal InfraRed Sensor", bands=[ \
                Band("1", "Coastal Aerosol", 25.0, BandType.REF, [0.433,0.453]), \
                Band("2", "Visible Blue", 25.0, BandType.REF, [0.45, 0.515]), \
                Band("3", "Visible Green", 25.0, BandType.REF, [0.525, 0.60]), \
                Band("4", "Visible Red", 25.0, BandType.REF, [0.63, 0.68]), \
                Band("5", "Near Infrared", 25.0, BandType.REF, [0.845, 0.885]), \
                Band("6", "Short-wave Infrared 1", 25.0, BandType.REF, [1.56, 1.66]), \
                Band("7", "Short-wave Infrared 2", 25.0, BandType.REF, [2.1, 2.3]), \
                Band("8", "Panchromatic", 12.5, BandType.PAN, [0.5, 0.68]), \
                Band("9", "Cirrus", 25.0, BandType.ATM, [1.360, 1.390]), \
                Band("quality", "Quality", 25.0, BandType.BQA, []), \
                Band("10", "Thermal Infrared 1", 25.0, BandType.THM, [10.3, 11.3]), \
                Band("11", "Thermal Infrared 2", 25.0, BandType.THM, [11.5, 12.5]) \
            ], other_props={ \
                "K1": 774.89,
                "K2": 1321.08,
                "acquisition_seconds": 24}) \
        ], {  # other spacecraft properties \
            "altitude": 705000.0, \
            "classification": "U", \
            "format": "FastL7A", \
            "inclination": 1.7139133254584316445390643346558, \
            "intl_designator": "13008A", \
            "satellite_name": "Landsat-8", \
            "name_pattern": "^L.*8$", \
            "nominal_pixel_degrees": 0.000250, \
            "number": 39084, \
            "omega": 0.001059, \
            "projection": "EQUIRECT_0_0", \
            "radius": 7285600.0, \
            "semi_major_axis": 7083160.0, \
            "solar_irrad_file": "solar_irrad_landsat8.txt", \
            "spectral_filter_file": "landsat8_vsir.flt", \
            "sweep_period": 0.07342143906020558002936857562408, \
            "tag": "LS8", \
            "tle_format": "L8%4d%sASNNOR.S00" } \
        )
     
    
    def __init__(self, sensors=[], other_props={}):
        """
        Construct a new member of the Spacecraft Enumeration.

        Parameters
        ----------
        sensors : list
            A list of Sensor object representing the Sensors carried by
            this Spacecraft

        other_props : dict
            A dictionary of other properties related to the Spacecraft
        """
        self.sensors = {}
        for sensor in sensors:
            self.sensors[sensor.id] = sensor
        self.__dict__.update(other_props)

    def __repr__(self):
        return "%s (Sensors: %s)" % (self, self.sensors.keys())

