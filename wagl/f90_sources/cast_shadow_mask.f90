SUBROUTINE get_proj_shadows(hx, hy, ns, nl, &
    htol, phi_sun, sun_zen, zmax, zmin, a, mask, h_offset, &
    n_inc, m_inc, aoff_x, aoff_y, nsA, nlA, k_setting, &
    dem_nr, dem_nc, nlA_ori, nsA_ori, ierr)

    implicit none

    integer ns, nl, dem_nr, dem_nc, nlA_ori, nsA_ori, k_setting
    integer k_max, n_add, m_add, ierr, az_case, i, j
!   NOTE: n_inc and m_inc are Floating Point arrays
    real n_inc(k_setting), m_inc(k_setting)
    real h_offset(k_setting)
    real phi_sun,zmax,zmin,sun_zen,htol
    real*8 hx, hy
    real d, d0, a(dem_nr, dem_nc)
    integer ncho, aoff_x, aoff_y, nsA, nlA
    integer t_aoff_x,t_aoff_y,t_nsA,t_nlA
    integer tval(8)
    logical rstatus
    integer*2 mask(nlA_ori, nsA_ori)
    integer*2 rmax, rmin

    real pi,r2d,d2r,eps
    common/base/pi,r2d,d2r,eps
!
!   calculate the border info for the sun position
!   In Australia and in the south in particular case=1
!   for Landsat since Landsat overpass is around 10:am
!   local time.
!   That is the sun azimuth is between East and North

    call set_borderf(rstatus,phi_sun,zmax,zmin,sun_zen,hx,hy,az_case, &
        d,d0,k_max,h_offset,n_inc,m_inc,n_add,m_add,k_setting,k_setting,ierr)

    if (.not.rstatus) then
        goto 99
    endif

!   define the maximum sized subset that can be processed
    if (az_case.eq.1) then
        t_aoff_x=0
        t_aoff_y=m_add
        t_nsA=ns-n_add
        t_nlA=nl-m_add
    else if (az_case.eq.2) then
        t_aoff_x=n_add
        t_aoff_y=m_add
        t_nsA=ns-n_add
        t_nlA=nl-m_add
    else if (az_case.eq.3) then
        t_aoff_x=n_add
        t_aoff_y=0
        t_nsA=ns-n_add
        t_nlA=nl-m_add
    else if (az_case.eq.4) then
        t_aoff_x=0
        t_aoff_y=0
        t_nsA=ns-n_add
        t_nlA=nl-m_add
    endif

!   Set the sub-matrix A where shade is to be found
!   aoff_x is the offset in samples
!   aoff_y is the offset in lines
!   pos in line in image for A(i,j) is aoff_x+j
!   pos in lines in image for A(i,j) is aoff_y+i

!   check the submatrix A is valid in various ways

!   first that the individual components are valid
    if (((aoff_x.lt.0) .or. (aoff_x.ge.ns)) .or. ((aoff_y.lt.0) .or. &
        (aoff_y.ge.nl)) .or. &
        ((nsA.lt.1) .or. (nsA.gt.ns)) .or. &
        ((nlA.lt.1) .or. (nlA.gt.nl))) then
            ierr = 71
            goto 99
    endif

!   Check A is embedded in the whole image
    if((aoff_x+nsA.gt.ns) .or. (aoff_y+nlA.gt.nl)) then
        ierr = 72
        goto 99
    endif

!   check the sub-image A plus the border area
!   needed to test A still fits inside the main image
!   with the buffer available
!   NOTE: treat the four cases in pairs

    if(az_case.eq.1.or.az_case.eq.2) then
        if (aoff_y.lt.m_add) then
            ierr = 73
            goto 99
        endif
    else if (az_case.eq.3.or.az_case.eq.4) then
        if (aoff_y+nlA+m_add.gt.nl) then
            ierr = 73
            goto 99
        endif
    endif

    if(az_case.eq.2.or.az_case.eq.3) then
        if (aoff_x.lt.n_add) then
            ierr = 74
            goto 99
        endif
    else if (az_case.eq.1.or.az_case.eq.4) then
        if (aoff_x+nsA+n_add.gt.ns) then
            ierr = 74
            goto 99
        endif
    endif
!    now set up the mask image to record shade pixels in
!    A NOTE: mask has the dimensions of A and not the
!    DEM - set as a 1-D array so that it can be indexed
!    in the subroutine without worrying about set bounds

!     Mask is a long integer (integer*4) due to recl issue
!     First set to 1 so zero will represent deep shadow

    do i=1,nlA
        do j=1,nsA
            mask(i,j)=1
        enddo
    enddo

    rmax=maxval(mask(1:nlA,1:nsA))
    rmin=minval(mask(1:nlA,1:nsA))

!   proj_terrain does the job of checking for occlusion
!   along the vector from a pixel to the sun
!   if any terrain obstructs the path the pixel is set to
!   zero in Mask

    call proj_terrain(ns, nl, nsA, nlA, a, mask, Aoff_x, Aoff_y, &
    k_max, n_inc, m_inc, h_offset, zmax, htol, dem_nr, dem_nc, nlA_ori, nsA_ori)

!   Description of proj_terrain
!   ===========================
!   call proj_terrain(n_max,m_max,n,m,z,mask,n_off,m_off,k_max,
!     . n_inc,m_inc,h_offset,zmax,htol,k_setting)
!
!   subroutine to construct mask of shade pixels
!
!   z(m_max,n_max) is the main array of heights
!   A(m,n) is the (sub-)matrix of target heights in:
!      z(m_off+1,n_off+1) to z(m_off+m,n_off+n)
!      mask(m,n) is the output mask
!      on input assumed to be 1 where valid data exist 0 else
!      k_max is the number of lags in the projection
!      n_inc(k_max) real increments column number for projection
!      m_inc(k_max)  real increments row number for the projection
!      h_offset(k_max) is the height of the projection hor each lag
!      zmax is the maximum altitude in the whole array
!      htol is a tolerance (m) for the test for a hit (RMS error in z)

99  continue
    return
END SUBROUTINE get_proj_shadows
