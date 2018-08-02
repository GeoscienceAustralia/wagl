! subroutine cal_track
SUBROUTINE cal_track(num,tin,orb_elements,spheroid,smodel, &
             track,istat)

!   Cal_Track calculates information at num points in time
!   along the satellite track. The times are in tin.
!   Results are put into the track common block

!   Re-written as an independent subroutine by JS, Aug 2014

!   Inputs:
!       num
!       tin
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
!       track
!           1. t
!           2. rho
!           3. phi_p
!           4. lam
!           5. beta
!           6. hxy
!           7. mj
!           8. skew
!
!   Outputs:
!       track
!       istat

    use sys_variables, only : pi, d2r, r2d

    implicit none

    integer, intent(in) :: num
    integer, intent(inout) :: istat
    double precision, dimension(num), intent(in) :: tin
    double precision, dimension(4), intent(in) :: spheroid
    double precision, dimension(3), intent(in) :: orb_elements
!   track(t,rho,phi_p,lam,beta,hxy,mj,skew)
    double precision, dimension(num,8), intent(out) :: track
    double precision, dimension(num) :: t, rho, phi_p, lam, beta
    double precision, dimension(num) :: hxy, mj, skew

!f2py depend(num), tin, track
!f2py depend(num), t, rho, phi_p, lam, beta, hxy, mj, skew

    double precision oi, ws
    double precision we

!   smodel(phi0,phi0_p,rho0,t0,lam0,gamm0,beta0,rotn0,hxy0,N0,H0,th_ratio0)
    double precision, dimension(12) :: smodel
    double precision rho0, t0, gamm0

    double precision psx, psy, psx_out, psy_out
    integer j

!   Initialise the return status
!   It will be overwritten by other routine calls
    istat = 0

    psx = 1.0d0/3600.0d0
    psy = 1.0d0/3600.0d0

!   smodel parameters
    rho0 = smodel(3)
    t0 = smodel(4)
    gamm0 = smodel(6)

!   Satellite orbital parameters
!   orbital inclination (degrees)
!   angular velocity (rad sec-1)
    oi = orb_elements(1)*d2r
    ws = orb_elements(3)

!   Spheroid parameters
!   Earth rotational angular velocity rad/sec
    we = spheroid(4)

    do j=1,num
        t(j) = tin(j)
        rho(j) = rho0+ws*(t(j)-t0)
        phi_p(j) = asin(cos(rho(j))*sin(oi))
        lam(j) = gamm0-pi/2.0d0+atan(tan(rho(j))/cos(oi))-we*(t(j)-t0)
        beta(j) = atan(-1.0d0/(tan(oi)*sin(rho(j))))
        skew(j) = atan(we*cos(phi_p(j))*cos(beta(j))/(ws+we* &
          cos(phi_p(j))*sin(beta(j))))
        call geo2metres_pixel_size(phi_p(j)*r2d,psx,psy,spheroid,&
               psx_out,psy_out,istat)
!        print*, 'psx_out,psy_out'
!        print*,psx_out,'|',psy_out
!        print*,'phi_p(j),r2d,psx,psy,spheroid,num,j'
!        print*,phi_p(j),'|',r2d,'|',psx,'|',psy,'|',spheroid,'|',num,'|',j
!        print*,'istat'
!        print*,istat
        hxy(j) = psx_out/psy_out
        if (j.gt.1) then
            mj(j-1) = (phi_p(j)-phi_p(j-1))/(lam(j)-lam(j-1))
        endif
    enddo

!   Store individual track variables into the track array
    track(:,1) = t
    track(:,2) = rho
    track(:,3) = phi_p
    track(:,4) = lam
    track(:,5) = beta
    track(:,6) = hxy
    track(:,7) = mj
    track(:,8) = skew

    return

END SUBROUTINE cal_track
