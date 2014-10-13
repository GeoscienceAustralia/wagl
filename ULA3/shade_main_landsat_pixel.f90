SUBROUTINE shade_main_landsat_pixel( &
    dem_data, solar_data, sazi_data, &
    dresx, dresy, spheroid, alat1, alon1, &
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
!   dresx is the x cell size.
!   dresy is the y cell size.
!   spheroid is the spheroidal parameters.
!       1. Spheroid major axis
!       2. Inverse flattening
!       3. Eccentricity squared
!       4. Earth rotational angular velocity rad/sec
!
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
    real*8 spheroid(4)
    real*8 dresx, dresy, alat1, alon1
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
    real*8 hx, hy
    real*4 zmax, zmin
    real*4 phi_sun
    real*4 sun_zen
    real*4 htol
    logical exists
    real pi, r2d, d2r, eps
    common/base/pi,r2d,d2r,eps

!f2py intent(in) dem_data, solar_data, sazi_data
!f2py intent(in) dres, spheroid, alat1, alon1
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
        hx = dresx
        hy = dresy
    else
        !       calculate longitude for each pixel of the line
        do j=1,ncol
            alon(j)=alon1+(j-1)*dresx
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
            call geo2metres_pixel_size(alat(ii), dresx, dresy, &
                                       spheroid, hx, hy, istat)
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
            call geo2metres_pixel_size(alat(ii), dresx, dresy, &
                                       spheroid, hx, hy, istat)
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
