! subroutine set_satmod
SUBROUTINE set_satmod(lon0,lat0,spheroid,orb_elements,smodel,istat)

!   lon0,lat0 are the (geo) coordinates of a point on the track
!   orb_elements 1,2,3 are orb_i,orb_R,nodal_p
!   spheroid 1,2 are asph,finv

!   * Re-written as an independent subroutine by JS, Aug 2014

!   Inputs:
!       lon0
!       lat0
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
!
!   Outputs:
!       smodel
!       istat

    use sys_variables, only : pi, d2r, r2d

    implicit none

    double precision, dimension(4), intent(in) :: spheroid
    double precision, dimension(3), intent(in) :: orb_elements
    double precision, intent(in) :: lon0,lat0
    double precision oi,orad,ws
    double precision asph,finv,e2,we

!   smodel(phi0,phi0_p,rho0,t0,lam0,gamm0,beta0,rotn0,hxy0,N0,H0,th_ratio0)
    double precision, dimension(12), intent(out) :: smodel
    double precision phi0,phi0_p,rho0,t0,lam0,gamm0,beta0,rotn0,hxy0
    double precision N0,H0,th_ratio0
    double precision rn0,temp,psx,psy,psx_out,psy_out,lonin,latin

    integer, intent(out) :: istat

!   Initialise the return status
    istat = 0

!   Satellite orbital parameters
!   orbital inclination (degrees)
!   semi_major radius (m)
!   angular velocity (rad sec-1)
    oi = orb_elements(1)*d2r
    orad = orb_elements(2)
    ws = orb_elements(3)

!   Spheroid parameters
!   Spheroid major axis
!   Inverse flattening
!   Eccentricity squared
!   Earth rotational angular velocity rad/sec
    asph = spheroid(1)
    finv = spheroid(2)
    e2 = spheroid(3)
    we = spheroid(4)

!   now set up the track based on the central point
    phi0 = lat0*d2r
    temp = 1.0d0-e2*sin(phi0)**2
    Rn0 = asph/sqrt(temp)
    phi0_p = phi0-asin(Rn0*e2*sin(phi0)*cos(phi0)/orad)
    rho0 = acos(sin(phi0_p)/sin(oi))
    t0 = (rho0+pi/2.0d0)/ws
    lam0 = lon0*d2r
    gamm0 = lam0+pi/2.0d0-atan(tan(rho0)/cos(oi))
    beta0 = atan(-1.0d0/(tan(oi)*sin(rho0)))
    rotn0 = atan(we*cos(phi0_p)*cos(beta0)/ &
        (ws+we*cos(phi0_p)*sin(beta0)))
    N0 = asph/sqrt(1.0d0-e2*sin(phi0)**2)
    H0 = orad-N0
    th_ratio0 = N0/H0
    psx = 1.0d0/3600.0d0
    psy = 1.0d0/3600.0d0
!    psx = 100.0d0/3600.0d0
!    psy = 100.0d0/3600.0d0
    lonin = lon0
    latin = lat0
    call geo2metres_pixel_size(latin,psx,psy,spheroid,psx_out,&
                               psy_out,istat)
    hxy0 = psx_out/psy_out

!   smodel(phi0,phi0_p,rho0,t0,lam0,gamm0,beta0,rotn0,hxy0,N0,H0,th_ratio0)
    smodel(1) = phi0
    smodel(2) = phi0_p
    smodel(3) = rho0
    smodel(4) = t0
    smodel(5) = lam0
    smodel(6) = gamm0
    smodel(7) = beta0
    smodel(8) = rotn0
    smodel(9) = hxy0
    smodel(10) = N0
    smodel(11) = H0
    smodel(12) = th_ratio0

    return

END SUBROUTINE set_satmod
