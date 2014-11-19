! subroutine solar
!subroutine solar(alat,along,iyear,imonth,iday,tt,solar_zen,sazi)
subroutine solar(alat,along,centry,tt,solar_zen,sazi)
!   program used to calculate solar angle using NOAA method
!   it is more accurate and using julian day
!   the program is written by Fuqin Li in 2010

!   * Re-written as an indepentent subroutine by JS, Aug 2014

!   Inputs:
!       alat
!       along
!       centry
!       tt
!       solar_zen
!       sazi

!   Outputs:
!       solar_zen
!       sazi

    use sys_variables, only : pi, d2r, r2d

    implicit none

    double precision alat, along
    double precision centry, tt
    real solar_zen, sazi
    double precision geom_ss
    double precision sun_long, sun_anom, eccent, sun_eqc
    double precision sun_truelong, sun_trueanom
    double precision sun_rad, sun_app, ecliptic_mean
    double precision obliq_corr, sun_asc, sun_declin
    double precision vary, eq_time, suntime_true, hangle
    double precision solar_zen1, solar_ele1, atmo
    double precision solar_ele, sss, ss_time

    geom_ss = 280.46646d0+centry*(36000.76983d0+centry*0.0003032d0)
    sun_long = dmod(geom_ss,360d0)
    sun_anom = 357.52911d0+centry*(35999.05029d0-0.0001537d0*centry)
    eccent = 0.016708634d0-centry*(0.000042037d0+0.0001537d0*centry)

    sun_eqc = sin(sun_anom*d2r)*(1.914602d0-centry*(0.004817d0+ &
      0.000014d0*centry))+sin(2*sun_anom*d2r)*(0.019993d0-0.000101d0* &
      centry)+sin(3*sun_anom*d2r)*0.000289d0

    sun_truelong = sun_long+sun_eqc
    sun_trueanom = sun_anom+sun_eqc

    sun_rad = (1.000001018d0*(1-eccent**2))/(1+eccent* &
      cos(sun_trueanom*d2r))

    sun_app = sun_truelong-0.00569d0-0.00478d0*sin(d2r*(125.04- &
      1937.136d0*centry))

    ecliptic_mean = 23.0+(26.0+((21.448-centry*(46.815+centry* &
      (0.00059d0-centry*0.001813d0))))/60.0)/60.0

    obliq_corr = ecliptic_mean+0.00256*cos((125.04-1934.136*centry)* &
      d2r)

    sun_asc = r2d*(atan2(cos(sun_app*d2r),cos(obliq_corr*d2r)* &
      sin(sun_app*d2r)))

    sun_declin = r2d*(asin(sin(obliq_corr*d2r)*sin(sun_app*d2r)))
    vary = tan(obliq_corr/2.0*d2r)**2

    eq_time = 4*r2d*(vary*sin(2.0*sun_long*d2r)-2*eccent*sin(d2r* &
      sun_anom)+ 4*eccent*vary*sin(d2r*sun_anom)*cos(2*d2r*sun_long) &
      -0.5*vary**2* sin(4*d2r*sun_long)-1.25*eccent**2*sin(2*d2r* &
      sun_anom))

    ss_time = tt*60+eq_time+4*along*r2d
    suntime_true = dmod(ss_time,1440d0)

    if (suntime_true/4.0 .ge. 0) then
        hangle = suntime_true/4.0-180
    else
        hangle = suntime_true/4.0+180
    endif

    solar_zen1 = r2d*(acos(sin(alat)*sin(d2r*sun_declin)+ &
      cos(alat)*cos(d2r*sun_declin)*cos(d2r*hangle)))

    solar_ele1 = 90.0-solar_zen1

    if (solar_ele1 .gt. 85) atmo = 0

    if (solar_ele1 .gt. 5 .and. solar_ele1 .le.85) atmo = 58.1/tan( &
      solar_ele1*d2r)-0.07/(tan(solar_ele1*d2r)**3)+0.000086d0/ &
      (tan(solar_ele1*d2r)**5)

    if (solar_ele1 .gt. -0.575 .and. solar_ele1 .le.5) atmo = &
      1735.0+solar_ele1*(-518.2+solar_ele1*(103.4+solar_ele1*(-12.79+ &
      solar_ele1*0.711)))

    if (solar_ele1 .le. -0.575 ) atmo = -20.772/tan(solar_ele1*d2r)

    atmo = atmo/3600.0
    solar_ele = solar_ele1+atmo
    solar_zen = 90-solar_ele

    if (hangle .gt. 0) then
        sss = r2d*(acos(((sin(alat)*cos(d2r*solar_zen1))-sin(d2r* &
          sun_declin))/(cos(alat)*sin(d2r*solar_zen1))))+180.0
    else
        sss = 540.0-r2d*(acos(((sin(alat)*cos(d2r*solar_zen1))- &
          sin(d2r*sun_declin))/(cos(alat)*sin(d2r*solar_zen1))))
    endif

    sazi = dmod(sss,360d0)

    return

end subroutine solar
