SUBROUTINE slope_aspect(nrow, ncol, nrow_alloc, ncol_alloc, &
    dresx, dresy, spheroid, alat, is_utm, dem, &
    theta, phit, ierr)
! this program is used to calculate slope and aspect angles
! using Sobel filter
! note: the row and column of DEM data must be larger
! than the image (extra each line and column for the four sides.
! it is needed for sobel filter.

!   Inputs:
!       spheroid
!           1. Spheroid major axis
!           2. Inverse flattening
!           3. Eccentricity squared
!           4. Earth rotational angular velocity rad/sec

    integer*4, intent(in) :: nrow, ncol, nrow_alloc, ncol_alloc
    double precision, intent(in) :: dresx, dresy
    double precision, dimension(4), intent(in) :: spheroid
    double precision, intent(in) :: alat
    real, dimension(nrow, ncol), intent(in) :: dem
    real, dimension(nrow_alloc, ncol_alloc) :: theta, phit
    logical is_utm
    integer*4, intent(out) ::  ierr


!    real*8 dresx, dresy
!    real*8 spheroid(4)
!    real*8 alat(nrow) ! remember that row starts at 2
!    real*4 dem(nrow, ncol) !
!    real*4 theta(nrow_alloc, ncol_alloc) !
!    real*4 phit(nrow_alloc, ncol_alloc) !
!    logical is_utm
!    integer*4 nrow_alloc, ncol_alloc, nrow, ncol, row, col
!    integer*4 ierr

! internal variables
    real*8 dx, dy
    real*8 pi, pia, pib
    real*8 p, q
    integer istat
    integer*4 n_row, n_col, row, col

!f2py depend(nrow, ncol), dem
!f2py depend(nrow_alloc, ncol_alloc), theta, phit

!feff2py intent(in) dresx, dresy, spheroid, alat, dem, is_utm
!feff2py intent(inout) theta, phit
!feff2py intent(out) ierr
!feff2py integer intent(hide),depend(solar) :: nrow_alloc=shape(solar,0), ncol_alloc=shape(solar,1)
!feff2py integer intent(hide),depend(dem) :: nrow=shape(dem,0), ncol=shape(dem,1)

!---------------------------------------------------------------------
!   upper left coordinator for satellite image.
!   here we assume that Landsat and DEM have same spatial resolution
!   and coordinator. The information can be obtained from satellite
!   header file. DEM header should add an extra pixel.
    pi=4.0d0*atan(1.0d0)
    pia=pi/180.0d0
    pib=180.0d0/pi
    ierr = 0;
!   check that the arrays dimensions are correct (relatively)
    if(ncol_alloc .ne. (ncol-2)) then
        ierr = 10
        return
    endif

    if(nrow_alloc .ne. (nrow-2)) then
        ierr = 11
        return
    endif

    n_col = ncol-2
    n_row = nrow-2

!   calculate pixel size in meters
    if(is_utm) then
        dx = dresx
        dy = dresy
    endif
!-------------------------------------------------------------------
    do col=2,n_col+1
        if(.not. is_utm) then
            call geo2metres_pixel_size(alat, dresx, dresy, &
                                       spheroid, dx, dy, istat)
        endif

        do row=2,n_row+1
!            We are taking as input Transposed arrays, as such we need to
!            reorder the convolution window operator.
!            p = &
!                      (dble(dem(row-1, col+1))-dble(dem(row-1, col-1)) + &
!                2.0d0*(dble(dem(row  , col+1))-dble(dem(row  , col-1))) + &
!                       dble(dem(row+1, col+1))-dble(dem(row+1, col-1))) / (8.0d0*dx)
!            q = &
!                      (dble(dem(row-1, col-1))-dble(dem(row+1, col-1)) + &
!                2.0d0*(dble(dem(row-1, col  ))-dble(dem(row+1, col  ))) + &
!                       dble(dem(row-1, col+1))-dble(dem(row+1, col+1))) / (8.0d0*dy)
             p = &
                      ((dble(dem(row+1, col-1)) - dble(dem(row-1, col-1))) + &
                 2.0d0*(dble(dem(row+1, col)) - dble(dem(row-1, col))) + &
                       (dble(dem(row+1, col+1)) - dble(dem(row-1, col+1)))) / (8.0*dx)
             q = &
                      ((dble(dem(row-1, col-1)) - dble(dem(row-1, col+1))) + &
                 2.0d0*(dble(dem(row, col-1)) - dble(dem(row, col+1))) + &
                       (dble(dem(row+1, col-1)) - dble(dem(row+1, col+1)))) / (8.0*dy)
    !----------------------------------------------------------------
            theta(row-1, col-1) = sngl(atan(sqrt(p**2 + q**2)) * pib)
            phit(row-1, col-1) = sngl(atan2(-p, -q) * pib)
            if (phit(row-1, col-1) .le. -180.0) phit(row-1, col-1) = phit(row-1, col-1) + 360.0
            if (phit(row-1, col-1) .gt.  180.0) phit(row-1, col-1) = phit(row-1, col-1) - 360.0
        enddo
    enddo
    return
END SUBROUTINE slope_aspect
