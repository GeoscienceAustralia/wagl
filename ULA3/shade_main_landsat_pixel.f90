SUBROUTINE shade_main_landsat_pixel( &
    dem_data, solar_data, sazi_data, &
    dres, alat1, alon1, &
    Aoff_x1, Aoff_x2, Aoff_y1, Aoff_y2, &
    nlA_ori, nsA_ori, &
    is_utm, &
    nrow, ncol, nl, ns, dem_nr, dem_nc, &
    a, solar, sazi, dem, alat, alon, mask, &
    ierr, mask_all)

    implicit none

!   NOTE THAT THIS CODE CANNOT BE COMPILED WITH -O3, AS IT PRODUCES
!   VERY DIFFERENT RESULTS WHEN IT IS. The only other optimisation
!   level I have tried is -O0.
!
!   Program to calculate cast shadow for a standard Landsat scene
!   the program was originally written by DLB Jupp in Oct. 2010
!   for a small sub_matrix and was modified by Fuqin Li in Oct.
!   2010 so that the program can be used for large landsat scene.
!
!   Basically, a sub-matrix A is embedded in a larger DEM image
!   and the borders must be large enough to find the shaded pixels.
!   If we assume the solar azimuth and zenith angles change very
!   little within the sub-matrix A, then the Landsat scene can be
!   divided into several sub_matrix.
!   For Australian region, with 0.00025 degree resolution, the
!   sub-marix A is set to 500x500
!
!   we also need to set extra DEM lines/columns to run the Landsat
!   scene. This will change with elevation difference within the
!   scene and solar zenith angle. For Australian region and Landsat
!   scene with 0.00025 degree resolution, the maximum extra lines
!   are set to 250 pixels/lines for each direction. This figure
!   shold be sufficient for everywhere and anytime in Australia.
!   thus the DEM image will be larger than landsat image for
!   500 lines x 500 columns
!
!   Current program operates in all 4 Azimuth cases (the four quadrants)
!
!   Arguments
!   ==========
!   dem_data is the dem data.
!   solar_data is the solar zenith angle data.
!   azi_data is the solar azimuth angle data.
!   nrow and ncol are the number of rows and columns in the region.
!   nl and ns are the number of rows and columns in the dem data
!   dem_nr and dem_nc are the number of rows and columns in 'one chunk' of a submatrix (includes the boundary padding).
!   dres is the cell size.
!   ierr provides a spot for the a return error code.
!   alat and alon are the lattitude and longitude of the origin of the region.
!   nlA_ori, nsA_ori are the sub-matrix lines and columns.
!   is_utm are the inputs in UTM (.true. == 'yes').
!   mask_all holds the result mask.
!   Aoff_x1 is the pixel number before the Landsat image and Aoff_x2 is pixel number after the Landsat image
!   Aoff_y1 is the line number before the Landsat image starts and Aoff_y2 is line number after teh Landsat image end
!
!   Errors
!   ======
!   in 20s to 50s: errors indexing various arrays.
!   in 60s: set_border
!   in 70s: get_proj_shadows
!   ierr = 61: 'azimuth case not possible - phi_sun must be in 0 to 360 deg'
!   ierr = 62: 'k_max gt k_setting'
!   ierr = 63: 'add outside add_max ranges'
!   ierr = 71: 'Parameters defining A are invalid'
!   ierr = 72: 'Matrix A not embedded in image'
!   ierr = 73: 'matrix A does not have sufficient y buffer'
!   ierr = 74: 'matrix A does not have sufficient x buffer'

    integer*4 k_setting
    parameter (k_setting=1500)

! arguments
    real*4 dem_data(nl, ns) !
    real*4 solar_data(nrow, ncol) !
    real*4 sazi_data(nrow, ncol) !
    real*8 dres, alat1, alon1
    integer*4 Aoff_x1, Aoff_x2, Aoff_y1, Aoff_y2
    integer*4 nlA_ori, nsA_ori
    logical is_utm
    integer*4 nrow, ncol, dem_nr, dem_nc, nl, ns
    real*4 a(dem_nr, dem_nc) !
    real*4 solar(nlA_ori, ncol) !
    real*4 sazi(nlA_ori, ncol) !
    real*4 dem(dem_nr, ns) !
    real*8 alat(nlA_ori) !
    real*8 alon(ncol) !
    integer*2 mask(nlA_ori, nsA_ori) !
    integer*2 mask_all(nrow, ncol) !
    integer*4 ierr

! internal variables
    integer*4 nsA, nlA, nchf, i, j, ii, jj
    integer*4 k, l, kkx, kky, nmax_sub, Mmax_sub
    real*4 n_inc(k_setting) !
    real*4 m_inc(k_setting) !
    real*4 h_offset(k_setting) !
    real*4 hx, hy
    real*4 zmax, zmin
    real*4 phi_sun
    real*4 sun_zen
    real*4 htol
    logical exists
    real pi, r2d, d2r, eps
    common/base/pi,r2d,d2r,eps

!f2py intent(in) dem_data, solar_data, sazi_data
!f2py intent(in) dres, alat1, alon1
!f2py intent(in) Aoff_x1, Aoff_x2, Aoff_y1, Aoff_y2
!f2py intent(in) nlA_ori, nsA_ori
!f2py intent(in) is_utm
!f2py integer intent(hide),depend(solar_data) :: nrow=shape(solar_data,0), ncol=shape(solar_data,1)
!f2py integer intent(hide),depend(dem_data) :: nl=shape(dem_data,0), ns=shape(dem_data,1)
!f2py integer intent(hide),depend(Aoff_y1,nlA_ori,Aoff_y2,Aoff_x1,nsA_ori,Aoff_x2) :: dem_nr=Aoff_y1+nlA_ori+Aoff_y2, dem_nc=Aoff_x1+nsA_ori+Aoff_x2
!f2py intent(hide) :: a, solar, sazi, dem, alat, alon, mask
!f2py intent(out) ierr
!f2py intent(out) mask_all

!   set basic constants
    pi=4.0*atan(1.0)
    r2d=180.0/pi
    d2r=pi/180.0
    eps=1.0e-7

!   set the tolerance for occlusion in metres
!   (usually >0 and less than 10.0m)
    htol=1.0

    if(is_utm) then
        hx = dres
        hy = dres
    else
        !       calculate longitude for each pixel of the line
        do j=1,ncol
            alon(j)=alon1+(j-1)*dres
        enddo
    endif

!--------------------------------------------------------------
!   kky for line and kkx for column
!   kky and kkx are the sub_marix number
    kky=int(nrow/nlA_ori)
    kkx=int(ncol/nsA_ori)

!write(*,*)"about to start main loop"
    do k=1, kky
!       calculate sub_marix DEM dimensions
        nlA=nlA_ori
        mmax_sub=nlA+Aoff_y1+Aoff_y2

!       not sure we really need to copy the array here - perhaps we could just index more cleverly below.
!write(*,*)"about to copy dem data"
!----------------------------bounds check----------------------------
        !dem(dem_nr, ns)
        !dem_data(nl, ns)
        if(mmax_sub .gt. dem_nr) then
            ierr = 20
            return
        endif
        if((k-1)*nlA_ori+mmax_sub .gt. nl) then
            ierr = 21
            return
        endif
!--------------------------end bounds check--------------------------
        do j=1,ns
            do i=1,mmax_sub
                dem(i, j) = dem_data((k-1)*nlA_ori+i, j)
            end do
        end do

        zmax=maxval(dem(1:mmax_sub,1:ns))
        zmin=minval(dem(1:mmax_sub,1:ns))

!write(*,*)"about to copy solar data"
!----------------------------bounds check----------------------------
        !solar(nlA_ori, ncol)
        !solar_data(nrow, ncol)
        if(nlA .gt. nlA_ori) then
            ierr = 22
            return
        endif
        if((k-1)*nlA_ori+nlA .gt. nrow) then
            ierr = 23
            return
        endif
!--------------------------end bounds check--------------------------
        do j=1,ncol
            do i=1,nlA
                solar(i, j) = solar_data((k-1)*nlA_ori+i, j)
            enddo
        enddo
        do j=1,ncol
            do i=1,nlA
                sazi(i, j) = sazi_data((k-1)*nlA_ori+i, j)
            enddo
        enddo

        if(.not.is_utm) then
!           calculate latitude for each line
            do i=1,nlA
                alat(i)=alat1-((k-1)*nlA_ori+i-1)*dres
            enddo
            ii=nlA/2
            call pixelsize(alat(ii),dres,hx,hy)
        endif


!       divide seveal sub_matrix according to columns
!write(*,*)"about to start cols: kkx = ",kkx
        do l=1,kkx
            nsA=nsA_ori
            nmax_sub=nsA+Aoff_x1+Aoff_x2

            jj=(l-1)*nsA_ori+nsA/2

            phi_sun=sazi(ii,jj)
!           NOTE zenith + 3 degrees
            sun_zen=solar(ii,jj)+3
!----------------------------bounds check----------------------------
            !dem(dem_nr, ns)
            !a(dem_nr, dem_nc)
            if(mmax_sub .gt. dem_nr) then
                ierr = 24
                return
            endif
            if((l-1)*nsA_ori+nmax_sub .gt. ns) then
                ierr = 25
                return
            endif
            if(nmax_sub .gt. dem_nc) then
                ierr = 26
                return
            endif
!--------------------------end bounds check--------------------------
            do j=1,nmax_sub
                do i=1,mmax_sub
                    a(i,j)=dem(i,(l-1)*nsA_ori+j)
                enddo
            enddo
!write(*,*)"about to call get_proj_shadows"
            call get_proj_shadows(hx, hy, nmax_sub, mmax_sub, &
            htol, phi_sun, sun_zen, zmax, zmin, a, mask, h_offset, &
            n_inc, m_inc, Aoff_x1, Aoff_y1, nsA, nlA, k_setting, &
            dem_nr, dem_nc, nlA_ori, nsA_ori, ierr)

!----------------------------bounds check----------------------------
            !mask(nlA_ori, nsA_ori)
            !mask_all(nrow, ncol)
            if((k-1)*nlA_ori+nlA .gt. nrow) then
                ierr = 27
                return
            endif

            if((l-1)*nsA_ori+nsA .gt. ncol) then
                ierr = 28
                return
            endif

            if(nlA .gt. nlA_ori) then
                ierr = 29
                return
            endif

            if(nsA .gt. nsA_ori) then
                ierr = 30
                return
            endif
!--------------------------end bounds check--------------------------
            do j=1,nsA
                do i=1,nlA
                    mask_all((k-1)*nlA_ori+i,(l-1)*nsA_ori+j)=mask(i,j)
                enddo
            enddo
        enddo

!       last column block (if required)
!write(*,*)"about to do last cols"
        if (ncol .gt. kkx*nsA_ori) then
            nsA=ncol-kkx*nsA_ori
            nmax_sub=nsA+Aoff_x1+Aoff_x2
            jj=kkx*nsA_ori+nsA/2
            phi_sun=sazi(ii,jj)
!           NOTE zenith + 3 degrees
            sun_zen=solar(ii,jj)+3

!----------------------------bounds check----------------------------
            !dem(dem_nr, ns)
            !a(dem_nr, dem_nc)
            if(mmax_sub .gt. dem_nr) then
                ierr = 31
                return
            endif
            if(nmax_sub .gt. dem_nc) then
                ierr = 32
                return
            endif
            if(kkx*nsA_ori+nmax_sub .gt. ns) then
                ierr = 33
                return
            endif
!--------------------------end bounds check--------------------------
            do i=1,mmax_sub
                do j=1,nmax_sub
                    a(i,j)=dem(i,kkx*nsA_ori+j)
                enddo
            enddo

            call get_proj_shadows(hx, hy, nmax_sub, mmax_sub, &
            htol, phi_sun, sun_zen, zmax, zmin, a, mask, h_offset, &
            n_inc, m_inc, Aoff_x1, Aoff_y1, nsA, nlA, k_setting, &
            dem_nr, dem_nc, nlA_ori, nsA_ori, ierr)

!----------------------------bounds check----------------------------
            !mask(nlA_ori, nsA_ori)
            !mask_all(nrow, ncol)
            if((k-1)*nlA_ori+nlA .gt. nrow) then
                ierr = 34
                return
            endif
            if(nlA .gt. nlA_ori) then
                ierr = 35
                return
            endif
            if(kkx*nsA_ori+nsA .gt. ncol) then
                ierr = 36
                return
            endif
            if(nsA .gt. nsA_ori) then
                ierr = 37
                return
            endif
!--------------------------end bounds check--------------------------
            do i=1,nlA
                do j=1,nsA
                    mask_all((k-1)*nlA_ori+i,kkx*nsA_ori+j)=mask(i,j)
                enddo
            enddo
        endif
    enddo

!write(*,*)"about to do last rows"
!   do the last rows (if required)
    if (nrow .gt. kky*nlA_ori) then
        nlA=nrow-kky*nlA_ori
        mmax_sub=nlA+Aoff_y1+Aoff_y2

!       not sure we really need to copy the array here - perhaps we could just index more cleverly below.
!----------------------------bounds check----------------------------
        !dem(dem_nr, ns)
        !dem_data(nl, ns)
        if(mmax_sub .gt. dem_nr) then
            ierr = 38
            return
        endif
        if(kky*nlA_ori+mmax_sub .gt. nl) then
            ierr = 39
            return
        endif
!--------------------------end bounds check--------------------------
        do i=1,mmax_sub
            do j=1,ns
                dem(i, j) = dem_data(kky*nlA_ori+i, j)
            end do
        end do

        zmax=maxval(dem(1:mmax_sub,1:ns))
        zmin=minval(dem(1:mmax_sub,1:ns))

!----------------------------bounds check----------------------------
        !solar(nlA_ori, ncol)
        !solar_data(nrow, ncol)
        if(nlA .gt. nlA_ori) then
            ierr = 40
            return
        endif
        if(kky*nlA_ori+nlA .gt. nrow) then
            ierr = 41
            return
        endif
!--------------------------end bounds check--------------------------
        do i=1,nlA
            do j=1,ncol
                solar(i,j) = solar_data(kky*nlA_ori+i, j)
            enddo
            do j=1,ncol
                sazi(i,j) = sazi_data(kky*nlA_ori+i, j)
            enddo
        enddo

        if(.not.is_utm) then
!           calculate latitude and longitude for sub_matrix
            do i=1,nlA
                alat(i)=alat1-(kky*nlA_ori+i-1)*dres
            enddo
        endif

        ii=nlA/2

        if(.not.is_utm) then
            call pixelsize(alat(ii),dres,hx,hy)
        endif

!       divide seveal sub_matrix according to columns
        do l=1,kkx
            nsA=nsA_ori
            nmax_sub=nsA+Aoff_x1+Aoff_x2

            jj=(l-1)*nsA_ori+nsA/2

            phi_sun=sazi(ii,jj)
!           NOTE zenith + 3 degrees
            sun_zen=solar(ii,jj)+3

!----------------------------bounds check----------------------------
            !dem(dem_nr, ns)
            !a(dem_nr, dem_nc)
            if(mmax_sub .gt. dem_nr) then
                ierr = 42
                return
            endif
            if(nmax_sub .gt. dem_nc) then
                ierr = 43
                return
            endif
            if((l-1)*nsA_ori+nmax_sub .gt. ns) then
                ierr = 44
                return
            endif
!--------------------------end bounds check--------------------------
            do i=1,mmax_sub
                do j=1,nmax_sub
                    a(i,j)=dem(i,(l-1)*nsA_ori+j)
                enddo
            enddo

            call get_proj_shadows(hx, hy, nmax_sub, mmax_sub, &
            htol, phi_sun, sun_zen, zmax, zmin, a, mask, h_offset, &
            n_inc, m_inc, Aoff_x1, Aoff_y1, nsA, nlA, k_setting, &
            dem_nr, dem_nc, nlA_ori, nsA_ori, ierr)

!----------------------------bounds check----------------------------
            !mask(nlA_ori, nsA_ori)
            !mask_all(nrow, ncol)
            if(kky*nlA_ori+nlA .gt. nrow) then
                ierr = 45
                return
            endif
            if(nlA .gt. nlA_ori) then
                ierr = 46
                return
            endif
            if((l-1)*nsA_ori+nsA .gt. ncol) then
                ierr = 47
                return
            endif
            if(nsA .gt. nsA_ori) then
                ierr = 48
                return
            endif
!--------------------------end bounds check--------------------------
            do i=1,nlA
                do j=1,nsA
                    mask_all(kky*nlA_ori+i,(l-1)*nsA_ori+j)=mask(i,j)
                enddo
            enddo
        enddo

        if (ncol .gt. kkx*nsA_ori) then
            nsA=ncol-kkx*nsA_ori
            nmax_sub=nsA+Aoff_x1+Aoff_x2
            jj=kkx*nsA_ori+nsA/2

            phi_sun=sazi(ii,jj)
!           NOTE zenith + 3 degrees
            sun_zen=solar(ii,jj)+3

!----------------------------bounds check----------------------------
            !dem(dem_nr, ns)
            !a(dem_nr, dem_nc)
            if(mmax_sub .gt. dem_nr) then
                ierr = 49
                return
            endif
            if(nmax_sub .gt. dem_nc) then
                ierr = 50
                return
            endif
            if(kkx*nsA_ori+nmax_sub .gt. ns) then
                ierr = 51
                return
            endif
!--------------------------end bounds check--------------------------
            do i=1,mmax_sub
                do j=1,nmax_sub
                    a(i,j)=dem(i,kkx*nsA_ori+j)
                enddo
            enddo

            call get_proj_shadows(hx, hy, nmax_sub, mmax_sub, &
            htol, phi_sun, sun_zen, zmax, zmin, a, mask, h_offset, &
            n_inc, m_inc, Aoff_x1, Aoff_y1, nsA, nlA, k_setting, &
            dem_nr, dem_nc, nlA_ori, nsA_ori, ierr)

!----------------------------bounds check----------------------------
            !mask(nlA_ori, nsA_ori)
            !mask_all(nrow, ncol)
            if(kky*nlA_ori+nlA .gt. nrow) then
                ierr = 52
                return
            endif
            if(nlA .gt. nlA_ori) then
                ierr = 53
                return
            endif
            if(kkx*nsA_ori+nsA .gt. ncol) then
                ierr = 54
                return
            endif
            if(nsA .gt. nsA_ori) then
                ierr = 55
                return
            endif
!--------------------------end bounds check--------------------------
            do i=1,nlA
                do j=1,nsA
                    mask_all(kky*nlA_ori+i,kkx*nsA_ori+j)=mask(i,j)
                enddo
            enddo
        endif
    endif


!write(*,*)"nrow = ",nrow
!write(*,*)"ncol = ",ncol
!write(*,*)"nlA_ori = ",nlA_ori
!write(*,*)"nsA_ori = ",nsA_ori
!write(*,*)"nl = ",nl
!write(*,*)"ns = ",ns
!write(*,*)"dem_nr = ",dem_nr
!write(*,*)"dem_nc = ",dem_nc
!write(*,*)"Aoff_x1 = ",Aoff_x1
!write(*,*)"Aoff_x2 = ",Aoff_x2
!write(*,*)"Aoff_y1 = ",Aoff_y1
!write(*,*)"Aoff_y2 = ",Aoff_y2
!write(*,*)"kky = ",kky
!write(*,*)"kkx = ",kkx

!write(*,*)"end of shade_main_landsat_pixel"
END SUBROUTINE shade_main_landsat_pixel





!---------------------------------------------------------------------------
subroutine set_borderf(set_border,phi_sun, zmax, zmin, sun_zen, hx, hy, &
    az_case, d, d0, k_max, h_offset, n_inc, m_inc, n_add, m_add, &
    k_setting, add_max, ierr)

    implicit none

!   set_border defines the buffer that is needed to find the
!   occluding terrain. This will be larger for low sun elevations

    integer k_max, k, n_add, m_add, k_setting, add_max, ierr, az_case
    real n_inc(k_setting), m_inc(k_setting)
    real h_offset(k_setting)
    real phi_sun, zmax, zmin, sun_zen, hx, hy, phc
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

    n_add=ifix(d*cosphc/hx+1.5)
    m_add=ifix(d*sinphc/hy+1.5)

    if ((n_add.gt.add_max .or. m_add.gt.add_max) .or. (n_add.lt.0.or.m_add.lt.0)) then
        ierr=63
        goto 99
    endif
    return
99  set_border=.false.
    return
end subroutine set_borderf





!-----------------------------------------------------------------
subroutine get_proj_shadows(hx, hy, ns, nl, &
    htol, phi_sun, sun_zen, zmax, zmin, a, mask, h_offset, &
    n_inc, m_inc, aoff_x, aoff_y, nsA, nlA, k_setting, &
    dem_nr, dem_nc, nlA_ori, nsA_ori, ierr)

    implicit none

    integer ns, nl, dem_nr, dem_nc, nlA_ori, nsA_ori, k_setting
    integer k_max, n_add, m_add, ierr, az_case, i, j
!   NOTE: n_inc and m_inc are Floating Point arrays
    real n_inc(k_setting), m_inc(k_setting)
    real h_offset(k_setting)
    real phi_sun,zmax,zmin,sun_zen,hx,hy,htol
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
    end






subroutine proj_terrain(n_max, m_max, n, m, z, mask, n_off, m_off, k_max, &
     n_inc, m_inc, h_offset, zmax, htol, dem_nr, dem_nc, nlA_ori, nsA_ori)

    implicit none

!   subroutine to construct mask of shade pixels

!   z(m_max,n_max) is the main array of heights
!   A(m,n) is the (sub-)matrix of target heights in:
!   z(m_off+1,n_off+1) to z(m_off+m,n_off+n)
!   mask(m,n) is the output mask
!   on input assumed to be 1 where valid data exist 0 else
!   k_max is the number of lags in the projection
!   n_inc(k_max) real increments column number for projection
!   m_inc(k_max) real increments row number for the projection
!   h_offset(k_max) is the height of the projection hor each lag
!   zmax is the maximum altitude in the whole array
!   htol is a tolerance (m) for the test for a hit (RMS error in z)

    integer dem_nr, dem_nc, nlA_ori, nsA_ori
    integer n_max, m_max, n, m, n_off, m_off, k_max
!   NOTE: n_inc and m_inc are Floating Point arrays
    real m_inc(k_max), n_inc(k_max)
    real h_offset(k_max)
    real z(dem_nr, dem_nc)
    integer*2 mask(nlA_ori, nsA_ori)
    real zmax, t, tt, test, htol, xd, yd, zpos
    integer i, j, ii, jj, k, ipos, jpos
    integer*4 itot
    real pi,r2d,d2r,eps
    common/base/pi,r2d,d2r,eps

!   loop over object matrix A and project out into buffer


!   The main 2 loops are over the pixels in the submatrix A

!   For given A(i,j) the search for occluding terrain is done along a "line"
!   in the sun direction. The search can stop when the terrain would
!   have to be higher than the maximum value to occlude the current test
!   pixel (i,j) in the sub-matrix A

    itot=0
    do 100 i=1,m
        ii=m_off+i
        do 110 j=1,n
            jj=n_off+j
            t=z(ii,jj)
            do 120 k=1,k_max
                tt=t+h_offset(k)
                if(tt .le. zmax+htol) then
                    ipos = ifix(float(ii)+m_inc(k))
                    jpos = ifix(float(jj)+n_inc(k))
                    yd = float(ii)+m_inc(k)-float(ipos)
                    xd = float(jj)+n_inc(k)-float(jpos)
                    zpos = z(ipos,jpos)*(1.0-yd)*(1.0-xd) + &
                        z(ipos,jpos+1)*xd*(1.0-yd) + &
                        z(ipos+1,jpos+1)*xd*yd + &
                        z(ipos+1,jpos)*(1.0-xd)*yd
                    test = tt-zpos
                    if (test.le.htol) then
                        mask(i,j) = 0
                        itot = itot+1
                        goto 125
                    endif
                else
                    go to 125
                endif
120         continue
125     continue
110 continue
100 continue
    return
end subroutine proj_terrain





subroutine pixelsize(rlat,dres,dx,dy)
!   subroutine is used to calculate pixel size (in meters) at latitude and longitude projection
    double precision aa,bb,cc,dd,ff,rlat,rlon,pi
    double precision pia,rr,dres,ddx,ddy
    real dx,dy

!   set projection parameters. here WGS84 is used
!   semi-major axis
    aa=6.3781370d6
!   flattening
    ff=2.98257223563d2
!   semi-minor axis
    bb=aa*(1.-1/ff)
    pi=4.0*atan(1.0)
    pia=pi/180.0
    cc=aa*cos(rlat*pia)
    dd=bb*sin(rlat*pia)
    rr=sqrt((aa**2*cc**2+bb**2*dd**2)/(cc**2+dd**2))
    ddy=dres*pia
    ddx=acos(sin(rlat*pia)**2+cos(rlat*pia)**2*cos(dres*pia))
    dy=rr*ddy
    dx=rr*ddx
    return
end subroutine pixelsize

