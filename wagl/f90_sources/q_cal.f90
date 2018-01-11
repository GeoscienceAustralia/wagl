! subroutine q_cal
subroutine q_cal(phip,orb_elements,spheroid,smodel,rhocal,tcal, &
             lamcal,betacal,istat)

!   base subroutine to calculate base track information
!   relative to a given geocentric latitude

!   * Re-written as an indepentent subroutine by JS, Aug 2014

!   Inputs:
!       phip
!       orb_elements
!           1. Orbital inclination (degrees)
!           2. Semi_major radius (m)
!           3. Angular velocity (rad sec-1)
!       spheroid
!           1. Spheroid major axis
!           2. Inverse flattening
!           3. Eccentricity squared
!           4. Earth rotational angular velocity rad/sec
!       smodel
!           1. phi0
!           2. phi0_p
!           3. rho0
!           4. t0
!           5. lam0
!           6. gamm0
!           7. beta0
!           8. rotn0
!           9. hxy0
!           10. N0
!           11. H0
!           12. th_ratio0
!
!   Outputs:
!       rhocal is angle from apogee (rad)
!       tcal is time from apogee (sec)
!       lamcal is longitude at track point (rad)
!       betacal is heading at track point (rad)
!       istat

    use sys_variables, only : pi, d2r, r2d

    implicit none

    double precision phip

    double precision orb_elements(3), spheroid(4)
    double precision oi, ws
    double precision we

    double precision smodel(12)
    double precision t0, gamm0

    double precision rhocal, tcal, lamcal, betacal
    integer istat

!   Satellite orbital parameters
!   orbital inclination (degrees)
!   angular velocity (rad sec-1)
    oi = orb_elements(1)*d2r
    ws = orb_elements(3)

!   Spheroid paramaters
!   Earth rotational angular velocity rad/sec
    we = spheroid(4)

!   smodel parameters
    t0 = smodel(4)
    gamm0 = smodel(6)

!   Initialise the return status
    istat = 0

    rhocal = acos(sin(phip)/sin(oi))
    tcal = (rhocal+pi/2.0d0)/ws
    lamcal = gamm0-pi/2.0d0+atan(tan(rhocal)/ &
      cos(oi))-we*(tcal-t0)
    betacal = atan(-1.0d0/(tan(oi)*sin(rhocal)))

    return

end subroutine q_cal
