subroutine set_borderf(set_border,phi_sun, zmax, zmin, sun_zen, hx, hy, &
    az_case, d, d0, k_max, h_offset, n_inc, m_inc, n_add, m_add, &
    k_setting, add_max, ierr)

    implicit none

!   set_border defines the buffer that is needed to find the
!   occluding terrain. This will be larger for low sun elevations

    integer k_max, k, n_add, m_add, k_setting, add_max, ierr, az_case
    real n_inc(k_setting), m_inc(k_setting)
    real h_offset(k_setting)
    real phi_sun, zmax, zmin, sun_zen, phc
    real*8 hx, hy
    real sinphc, cosphc, d, d0
    real pi, r2d, d2r, eps
    common/base/pi,r2d,d2r,eps
    logical set_border


    set_border=.true.
    ierr=0

!   the case is dependent on the sun azimuth
!   this defines which borders need to be available
!   NOTE: azimuth must be in degrees between 0 and 360 and clockwise from N

    if (phi_sun.ge.0.0 .and. phi_sun.le.90.0) then
        az_case=1
        phc=pi/2.0-phi_sun*d2r
    else if (phi_sun.ge.270.0 .and. phi_sun.le.360.0) then
        az_case=2
        phc=phi_sun*d2r-3.0*pi/2.0
    else if (phi_sun.ge.180.0 .and. phi_sun.le.270.0) then
        az_case=3
        phc=3.0*pi/2.0-phi_sun*d2r
    else if (phi_sun.ge.90.0 .and. phi_sun.le.180.0) then
        az_case=4
        phc=phi_sun*d2r-pi/2.0
    else
        ierr=61
        goto 99
    endif

!   calculate the border info
!   and the increments for the projection line

    sinphc=sin(phc)
    cosphc=cos(phc)

!   d is the plane distance from the pixel
!   with zmin where if the pixel at the
!   distance in the sun direction is at Zmax
!   then it can just occlude the first pixel
!
    d = (zmax-zmin)/tan(pi/2.0-sun_zen*d2r)
    if(cosphc*hy.gt.sinphc*hx) then
        d0 = 0.5*hx/cosphc
    else
        d0 = 0.5*hy/sinphc
    endif

!   d0 is the basic step in the d direction that increments
!   0.5 pixel in the larger of the two directions

!   k_max is the number of steps to get out to d
!   n_inc(k) and m_inc(k) are the sample and line increments
!   h_offset(k) on exit is the altitude above the start of the
!   sun beam as you move towards the sun

    k_max=ifix(d/d0+0.5)+1
    if (k_max.gt.k_setting) then
        ierr=62
        goto 99
    endif
    do k=1,k_max
        h_offset(k)=float(k)*d0

!       the increments to search for occlusions have sign depending
!       on the azimuth case

        if (az_case.eq.1.or.az_case.eq.4) then
            n_inc(k)=h_offset(k)*cosphc/hx
        else
            n_inc(k)=-h_offset(k)*cosphc/hx
        endif
        if (az_case.eq.3.or.az_case.eq.4) then
            m_inc(k)=h_offset(k)*sinphc/hy
        else
            m_inc(k)=-h_offset(k)*sinphc/hy
        endif
        h_offset(k)=h_offset(k)*tan(pi/2.0-sun_zen*d2r)
    enddo

!   n_add and m_add are the sizes of the pixel buffer
!   in sample and line dimensions

!   the two buffers will be on sides of the target image defined
!   by the azimuth case

    n_add=ifix(sngl(d*cosphc/hx+1.5))
    m_add=ifix(sngl(d*sinphc/hy+1.5))

    if ((n_add.gt.add_max .or. m_add.gt.add_max) .or. (n_add.lt.0.or.m_add.lt.0)) then
        ierr=63
        goto 99
    endif
    return
99  set_border=.false.
    return
end subroutine set_borderf

