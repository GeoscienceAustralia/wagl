! subroutine geo2metres_pixel_size
SUBROUTINE geo2metres_pixel_size(ycent,hx,hy,spheroid, &
             hx_out,hy_out,istat)

!   integer subroutine to get pixel size for geographic coordinates
!   depends only on latitude ycent
!   ycent, hx and hy are in degrees
!   hx_out and hy_out are in metres
!   spheroid is wgs84 which is appropriate here

!   * Re-written as an indepentent subroutine by JS, Aug 2014
!   * Probable overlap with subroutine pixelsize(rlat,dres,dx,dy)???
!   * pixelsize(rlat,dres,dx,dy) & pixelsize(rlat,dres,dx,dy,pia)
!     have now been replaced by this routine which itself has
!     been renamed to geo2metres_pixel_size which describes the
!     routine a little better.
!   * The routine now takes a 4 element spheroid parameter rather
!     than hard coded WGS84 parameters.

!   Inputs:
!       ycent
!       hx
!       hy
!       spheroid
!           1. Spheroid major axis
!           2. Inverse flattening
!           3. Eccentricity squared
!           4. Earth rotational angular velocity rad/sec
!
!   Outputs:
!       hx_out
!       hy_out
!       istat

    use sys_variables, only : pi, d2r, r2d

    implicit none

    double precision ycent
    double precision hx, hy, hx_out, hy_out

    double precision spheroid(4)
    double precision asph, finv

    double precision yc, a_e, b_e, p_e, q_e, Earth_Radius, d_phi

    integer istat

!   Spheroid parameters
!   Spheroid major axis
!   Inverse flattening
    asph = spheroid(1)
    finv = spheroid(2)

    yc = ycent*d2r
    a_e = asph
    b_e = a_e*(1.0d0-1.0d0/finv)
    p_e = a_e*cos(yc)
    q_e = b_e*sin(yc)
    Earth_Radius = sqrt(((a_e*p_e)**2+(b_e*q_e)**2)/(p_e**2+q_e**2))
    d_phi = (sin(yc))**2+cos(hx*d2r)*(cos(yc))**2
    d_phi = acos(d_phi)

    hx_out = d_phi*Earth_Radius
    hy_out = (hy*d2r)*Earth_Radius

    istat = 0

    return

END SUBROUTINE geo2metres_pixel_size
