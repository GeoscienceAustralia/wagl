import logging, numpy, math
import ULA3.geodesic as geodesic

logger = logging.getLogger('root.' + __name__)

# Solar parameters
DEG2RAD = math.pi / 180.0
RAD2DEG = 180.0 / math.pi
# Coefficients for time polynomial (possol.f)
AT = (0.000075, 0.001868, 0.032077, 0.014615, 0.040849,)
# Coefficients for Solar declination polynomial (possol.f)
AD = (0.006918, 0.399912, 0.070257, 0.006758, 0.000907, 0.002697, 0.001480,)

def eval_sol_grids(lon_array, lat_array, time_array, day_of_year):
    """Calculate solar zenith and azimuth angle grids.

    Code follows possol.f (QDERM)

    Arguments:
        lon_array: longitude (numpy array, radians)
        lat_array: latitude (numpy array, radians)
        time_array: time (numpy array, float)
        day_of_year: day of year (int)

    Returns:
        2-tuple (solar_zenith, solar_azimuth)
    """

    tet = float(day_of_year) * ((math.pi * 2) / 365)
    tet2 = tet * 2
    tet3 = tet * 3
    ct = math.cos(tet)
    ct2 = math.cos(tet2)
    ct3 = math.cos(tet3)
    st = math.sin(tet)
    st2 = math.sin(tet2)
    st3 = math.sin(tet3)

    # possol.f
    # et=a1+a2*cos(tet)-a3*sin(tet)-a4*cos(2.*tet)-a5*sin(2.*tet)

    et = ((AT[0] + AT[1]*ct - AT[2]*st - AT[3]*ct2 - AT[4]*st2) * (float(12 * 60) / math.pi))

    # Hour angle

    # Time values: DOY.fractional_day
    #ftime, fday = math.modf(time_array[0, 0])
    #h = ((time_array - fday) + lon*RAD2DEG/15 + (et/60) - 12) * 15 * DEG2RAD

    # Time values: tu
    h = (time_array + lon_array*RAD2DEG/15 + (et/60) - 12) * 15 * DEG2RAD

    # DEBUG
    #grid.dump(h * RAD2DEG, 'test_H.tif')

    logger.debug('HOUR ANGLE = %s', h)

    # Solar declination
    # possol.f
    # delta=b1-b2*cos(tet)+b3*sin(tet)-b4*cos(2.*tet)+b5*sin(2.*tet)-
    # &b6*cos(3.*tet)+b7*sin(3.*tet)

    delta = (AD[0] - AD[1]*ct  + AD[2]*st - AD[3]*ct2 + AD[4]*st2 - AD[5]*ct3 + AD[6]*st3)
    delta_cos = math.cos(delta)
    delta_sin = math.sin(delta)

    logger.debug('SOLAR DECLINATION = %srad, %sdeg', delta, math.degrees(delta))

    # Geocentric latitude

    latgc = numpy.arctan(numpy.tan(lat_array) * (1.0 - geodesic.earth.ECC2))

    logger.debug('GEOCENTRIC LAT = %s', latgc)

    # Solar zenith/elevation angle

    cos_sol_z = numpy.sin(latgc) * delta_sin + numpy.cos(latgc) * delta_cos * numpy.cos(h)
    sol_z = numpy.arccos(cos_sol_z)

    logger.debug('SOLAR ZENITH = %s', sol_z)

    # Solar azimuth angle

    #azimuth = numpy.arcsin( numpy.sin(h) * delta_cos / numpy.cos(elevation) )

    sol_az = numpy.arccos( (delta_sin - cos_sol_z * numpy.sin(latgc)) /
                           (numpy.sin(sol_z) * numpy.cos(latgc)) )
    #sol_az = numpy.where(h < 0, -sol_az, sol_az)
    numpy.sqrt(sol_az*sol_az, sol_az)

    logger.debug('SOLAR AZIMUTH = %s', sol_az)

#    # Beta (Soler-Eisemann eqn 6)
#
#    beta = numpy.arccos(numpy.sin(latgc) / orbit.INCL_SIN)
#    #beta = numpy.arctan(-1.0 / (numpy.sin(beta) * orbit.INCL_TAN))
#    numpy.arctan(-1.0 / (numpy.sin(beta) * orbit.INCL_TAN), beta)
#
#    print
#    print 'BETA (Soler-Eisemann eqn 6)'
#    print beta
#    grid.dump(beta, 'test_BETA_SE.tif')

    return (sol_z, sol_az)

"""
possol.f (for reference)
========================

      subroutine possol (month,jday,tu,xlon,xlat,
     a                   asol,phi0)

      real    tu,xlon,xlat,asol,phi0
      integer month,jday,ia,nojour

c     solar position (zenithal angle asol,azimuthal angle phi0
c                     in degrees)
c     jday is the number of the day in the month

      ia = 0
      call day_number(jday,month,ia,nojour)

      call  pos_fft (nojour, tu, xlon, xlat, asol, phi0)

c      if(asol.gt.90) call print_error(
c     s 'The sun is not raised')
c      return
      end

      subroutine day_number(jday,month,ia,j)
      integer jday, month, ia, j

      if (month.le.2) then
                      j=31*(month-1)+jday
              return
              endif
      if (month.gt.8) then
                      j=31*(month-1)-((month-2)/2)-2+jday
              else
                      j=31*(month-1)-((month-1)/2)-2+jday
              endif
      if(ia.ne.0 .and. mod(ia,4).eq.0) j=j+1
      return
      end

      subroutine pos_fft (j,tu,xlon,xlat,asol,phi0)
      real    tu, xlat, asol,phi0, tsm, xlon,xla, xj, tet,
     a          a1, a2, a3, a4, a5, et, tsv, ah, b1, b2, b3, b4,
     a          b5, b6, b7, delta, amuzero, elev, az, caz, azim, pi2
      integer j
      parameter (pi=3.14159265,fac=pi/180.)
c     solar position (zenithal angle asol,azimuthal angle phi0
c                     in degrees)
c     j is the day number in the year
c
c    mean solar time (heure decimale)

      tsm=tu+xlon/15.
      xla=xlat*fac
      xj=float(j)
      tet=2.*pi*xj/365.

c    time equation (in mn.dec)
      a1=.000075
      a2=.001868
      a3=.032077
      a4=.014615
      a5=.040849
      et=a1+a2*cos(tet)-a3*sin(tet)-a4*cos(2.*tet)-a5*sin(2.*tet)
      et=et*12.*60./pi

c     true solar time

      tsv=tsm+et/60.
      tsv=(tsv-12.)

c     hour angle

      ah=tsv*15.*fac

c     solar declination   (in radian)

      b1=.006918
      b2=.399912
      b3=.070257
      b4=.006758
      b5=.000907
      b6=.002697
      b7=.001480
      delta=b1-b2*cos(tet)+b3*sin(tet)-b4*cos(2.*tet)+b5*sin(2.*tet)-
     &b6*cos(3.*tet)+b7*sin(3.*tet)

c     elevation,azimuth

      amuzero=sin(xla)*sin(delta)+cos(xla)*cos(delta)*cos(ah)
      elev=asin(amuzero)
      az=cos(delta)*sin(ah)/cos(elev)
      if ( (abs(az)-1.000).gt.0.00000) az = sign(1.,az)
      caz=(-cos(xla)*sin(delta)+sin(xla)*cos(delta)*cos(ah))/cos(elev)
      azim=asin(az)
      if(caz.le.0.) azim=pi-azim
      if(caz.gt.0.and.az.le.0) azim=2*pi+azim
      azim=azim+pi
      pi2=2*pi
      if(azim.gt.pi2) azim=azim-pi2
      elev=elev*180./pi

c     conversion in degrees

      asol=90.-elev
      phi0=azim/fac
      return
      end
"""







#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
# these appear to be unused
#--------------------------------------------------------------------
# def extract_solar_distance(solar_dist_file, DOY):
#     """Extract Earth-Sun distance for this day of the year (varies during orbit)"""
#
#     logger.debug("Reading Solar Distance File: %s" % solar_dist_file)
#
#     fpIn = open(solar_dist_file)
#
#     line = fpIn.readline()
#     if not line:
#         logger.error("ERROR: Solar Distance file %s empty" % solar_dist_file)
#         return 0
#
#     solar_dist = 0
#     while line:
#         temp = line.split()
#         doy = int(temp[0])
#
#         if doy == DOY:
#             solar_dist = float(temp[1])
#         line = fpIn.readline()
#
#     fpIn.close()
#
#     assert solar_dist, 'Solar distance not found for DOY %s in file %s' % (DOY, solar_dist_file)
#
#     logger.debug("Solar Distance: for DOY:%d %f" % (DOY, solar_dist))
#
#     return solar_dist
#
#
#
#
#
# def extract_solar_irrad(solar_irrad_file, refective_bands):
#     """Extract solar irradiance values from the specified file. One for each band
#     Arguments:
#         solar_irrad_file - Lookup file containing solar irradiance values
#         refective_bands - List of reflective bands
#     Returns:
#         solar_irrad - Dict containing value for each reflective band
#     """
#
#     logger.debug("Reading Solar Irradiance File: %s" % solar_irrad_file)
#     fpIn = open(solar_irrad_file)
#
#     line = fpIn.readline()
#     if not line:
#         logger.error("ERROR: Satellite Solar Irrad file %s empty" % solar_irrad_file )
#         return None
#
#     if line.find("band solar irradiance") != 0:
#         logger.error('ERROR: not a "band solar irradiance" file')
#         return None
#
#     solar_irrad = {}
#
#     band_index = 0
#     line = fpIn.readline()
#     while line and band_index < len(refective_bands):
#         temp = line.split()
#         band_number = refective_bands[band_index] # TODO: Need to enhance file for LDCM
#         solar_irrad_value = float(temp[1])
#         solar_irrad[band_number] = float(temp[1])
#         # Skip 0.0 value for thermal band 6 in file
#         if solar_irrad_value:
#             logger.debug('solar_irrad[%s] = %s', band_number, solar_irrad[band_number])
#             band_index += 1
#         line = fpIn.readline()
#
#     fpIn.close()
#
#     return solar_irrad
#
#
#
#
#
#
#
#
# def hour_angle(lon, lat, time, day_of_year):
#     """Calculate hour angle grid.
#
#     Arguments:
#         lon: longitude (radians)
#         lat: latitude (radians)
#         time: time (float, fractional day)
#         day_of_year: day of year (int)
#
#     Returns:
#         numpy array (float)
#     """
#
#     tet = float(day_of_year) * ((math.pi * 2) / 365)
#     tet2 = tet * 2
#     ct = math.cos(tet)
#     ct2 = math.cos(tet2)
#     st = math.sin(tet)
#     st2 = math.sin(tet2)
#
#     et = ((AT[0] + AT[1]*ct - AT[2]*st - AT[3]*ct2 - AT[4]*st2)
#           * (float(12 * 60) / math.pi))
#
#     return (time + lon*RAD2DEG/15.0 + (et/60) - 12) * 15.0 * DEG2RAD
#--------------------------------------------------------------------
#--------------------------------------------------------------------
