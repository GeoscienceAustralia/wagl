! subroutine set_times
SUBROUTINE set_times(ymin,ymax,ntpoints,spheroid,orb_elements, &
             smodel,psx,psy,track,istat)

!   Calculate satellite track times and other info

!   Re-written an independent subroutine by JS, Aug 2014

!   Inputs:
!       ymin
!       ymax
!       ntpoints
!       spheroid
!           1. Spheroid major axis
!           2. Inverse flattening
!           3. Eccentricity squared
!           4. Earth rotational angular velocity rad/sec
!       orb_elements
!           1. Orbital inclination (degrees)
!           2. Semi_major radius (m)
!           3. Angular velocity (rad sec-1)
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
!       track
!           1. t
!           2. rho
!           3. phi_p
!           4. lam
!           5. beta
!           6. hxy
!           7. mj
!           8. skew
!       psx
!           Approximate pixel size (in degrees longitude)
!       psy
!           Approximate pixel size (in degrees latitude)
!
!   Outputs:
!       track
!       istat

    use sys_variables, only : pi, d2r, r2d

    implicit none

    double precision, intent(in) :: ymin, ymax
    integer, intent(in) :: ntpoints
    integer, intent(out) :: istat

    double precision, dimension(4), intent(in) :: spheroid
    double precision, dimension(3), intent(in) :: orb_elements
    double precision, dimension(ntpoints,8), intent(out) :: track
    double precision, dimension(ntpoints) :: tin
    double precision, intent(in) :: psx, psy

!f2py depend(ntpoints), tin, track

    double precision, dimension(12), intent(in) :: smodel
    double precision t_min, t_max, phip_max, phip_min, rhocal, tcal
    double precision lamcal, betacal, delta_t
    integer j

    istat = 0

!   set up the time range to use for the satellite track
    call geod2geo(ymax,orb_elements,spheroid,phip_max,istat)
    call q_cal(phip_max,orb_elements,spheroid,smodel,rhocal,tcal, &
                  lamcal,betacal,istat)
    t_min = tcal-5.0d0

    call geod2geo(ymin,orb_elements,spheroid,phip_min,istat)
    call q_cal(phip_min,orb_elements,spheroid,smodel,rhocal,tcal, &
                  lamcal,betacal,istat)
    t_max = tcal+5.0d0

!   The track has 12 points about 4 sec apart in this case
!   More could be used if accuracy questioned - but it is
!   very accurate

    delta_t = (t_max-t_min)/dble(ntpoints-1)
    do j=1,ntpoints
        tin(j) = dble(j-1)*delta_t+t_min
    enddo

!    print*,'Track starting time(sec)=',t_min
!    print*,'Track ending time(sec)=',t_max
!    print*,'Number of track points=',ntpoints
!    print*,'Track time step(sec)=',delta_t

    call cal_track(ntpoints,tin,orb_elements,spheroid,smodel, &
           psx,psy,track,istat)

    return

END SUBROUTINE set_times
