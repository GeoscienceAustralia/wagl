SUBROUTINE slope_pixelsize_newpole( &
    dres, alat, is_utm, &
    dem, solar, view, sazi, azi, &
    mask, theta, phit, it, et, azi_it, azi_et, rela, &
    nrow_alloc, ncol_alloc, nrow, ncol, ierr)
! this program is used to calculate slope and aspect angles
! using Sobel filter and then calculate incident and
! exiting angles as well as their azimuth angles.
! note: the row and column of DEM data must be larger
! than the image (extra each line and column for the four sides.
! it is needed for sobel filter.
    real*8 dres
	real*8 alat(nrow) ! remember that row starts at 2
	real*4 dem(nrow, ncol) !
    real*4 solar(nrow_alloc, ncol_alloc) !
    real*4 view(nrow_alloc, ncol_alloc) !
    real*4 sazi(nrow_alloc, ncol_alloc) !
    real*4 azi(nrow_alloc, ncol_alloc) !
    integer*2 mask(nrow_alloc, ncol_alloc) !
    real*4 theta(nrow_alloc, ncol_alloc) !
    real*4 phit(nrow_alloc, ncol_alloc) !
    real*4 it(nrow_alloc, ncol_alloc) !
    real*4 et(nrow_alloc, ncol_alloc) !
    real*4 azi_it(nrow_alloc, ncol_alloc) !
    real*4 azi_et(nrow_alloc, ncol_alloc) !
    real*4 rela(nrow_alloc, ncol_alloc) !
    logical is_utm
    integer*4 nrow_alloc, ncol_alloc, nrow, ncol, row, col
    integer*4 ierr

! internal variables
    real*8 dx, dy
    real*8 pi, pia, pib
    real*8 p, q

!f2py intent(in) dres, alat, dem, solar, view, sazi, azi, is_utm
!f2py intent(out) mask, theta, phit, it, et, azi_it, azi_et, rela
!f2py intent(out) ierr
!f2py integer intent(hide),depend(solar) :: nrow_alloc=shape(solar,0), ncol_alloc=shape(solar,1)
!f2py integer intent(hide),depend(dem) :: nrow=shape(dem,0), ncol=shape(dem,1)

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

    ncol = ncol-2
    nrow = nrow-2

!   calculate pixel size in meters
    if(is_utm) then
        dx = dres
        dy = dres
    endif
!-------------------------------------------------------------------
    do row=2,nrow+1
        if(.not. is_utm) then
            call pixelsize(alat(row), dres, dx, dy, pia)
        endif

        do col=2,ncol+1
            mask(row-1, col-1)=1
            if (sazi(row-1, col-1) .le. -180.0) sazi(row-1, col-1)=sazi(row-1, col-1)+360.0
            if (sazi(row-1, col-1) .gt.  180.0) sazi(row-1, col-1)=sazi(row-1, col-1)-360.0
            if ( azi(row-1, col-1) .le. -180.0)  azi(row-1, col-1)= azi(row-1, col-1)+360.0
            if ( azi(row-1, col-1) .gt.  180.0)  azi(row-1, col-1)= azi(row-1, col-1)-360.0
            p = &
                      (dble(dem(row-1,col+1))-dble(dem(row-1,col-1)) + &
                2.0d0*(dble(dem(row  ,col+1))-dble(dem(row  ,col-1))) + &
                       dble(dem(row+1,col+1))-dble(dem(row+1,col-1))) / (8.0d0*dx)
            q = &
                      (dble(dem(row-1,col-1))-dble(dem(row+1,col-1)) + &
                2.0d0*(dble(dem(row-1,col  ))-dble(dem(row+1,col  ))) + &
                       dble(dem(row-1,col+1))-dble(dem(row+1,col+1))) / (8.0d0*dy)
    !----------------------------------------------------------------
            theta(row-1, col-1)=sngl(atan(sqrt(p**2+q**2))*pib)
            phit(row-1, col-1)=sngl(atan2(-p,-q)*pib)
            if (phit(row-1, col-1) .le. -180.0) phit(row-1, col-1)=phit(row-1, col-1)+360.0
            if (phit(row-1, col-1) .gt.  180.0) phit(row-1, col-1)=phit(row-1, col-1)-360.0
            call cal_pole(solar(row-1, col-1), sazi(row-1, col-1), &
            theta(row-1, col-1), phit(row-1, col-1), it(row-1, col-1), azi_it(row-1, col-1))
            call cal_pole( view(row-1, col-1),  azi(row-1, col-1), &
            theta(row-1, col-1), phit(row-1, col-1), et(row-1, col-1), azi_et(row-1, col-1))
            if (azi_it(row-1, col-1).le. -180.0) azi_it(row-1, col-1)=azi_it(row-1, col-1)+360
            if (azi_it(row-1, col-1).gt.  180.0) azi_it(row-1, col-1)=azi_it(row-1, col-1)-360
            if (azi_et(row-1, col-1).le. -180.0) azi_et(row-1, col-1)=azi_et(row-1, col-1)+360
            if (azi_et(row-1, col-1).gt.  180.0) azi_et(row-1, col-1)=azi_et(row-1, col-1)-360
            rela(row-1, col-1)=azi_it(row-1, col-1)-azi_et(row-1, col-1)
    !----------------------------------------------------------------
            if (rela(row-1, col-1) .le. -180.0) rela(row-1, col-1)=rela(row-1, col-1)+360.0
            if (rela(row-1, col-1) .gt. 180.0) rela(row-1, col-1)=rela(row-1, col-1)-360.0
            if (cos(it(row-1, col-1)*pia) .le. 0.0) then
                mask(row-1, col-1)=0
            endif
            if (cos(et(row-1, col-1)*pia) .le.0.0 ) then
                mask(row-1, col-1)=0
            endif
        enddo
    enddo
    return
END SUBROUTINE slope_pixelsize_newpole





SUBROUTINE pixelsize(rlat, dres, dx, dy, pia)
!   subroutine is used to calculate pixel size (in meters) at latitude and
!   longitude projection
    double precision aa,bb,cc,dd,ff,rlat,rlon,pi
    double precision pia,rr,dres,dx,dy,ddx,ddy
!   set projection parameters. here WGS84 is used
!   semi-major axis
    aa=6.3781370d6
!   flattening
    ff=2.98257223563d2
!   semi-minor axis
    bb=aa*(1.0d0-1.0d0/ff)
    cc=aa*cos(rlat*pia)
    dd=bb*sin(rlat*pia)
    rr=sqrt((aa**2*cc**2+bb**2*dd**2)/(cc**2+dd**2))
    ddy=dres*pia
    ddx=acos(sin(rlat*pia)**2+cos(rlat*pia)**2*cos(dres*pia))
    dy=rr*ddy
    dx=rr*ddx
    return
END SUBROUTINE pixelsize





SUBROUTINE cal_pole(theta, phi, theta_p, phi_p, thp, php)
! does some funky stuff to calculate various angles
    real theta, phi, theta_p, phi_p, thp, php, offset
    real pdiff, costhp, sinphp, cosphp, tnum, tden
    real eps, pi, d2r

    eps=1.0e-6
    pi=4.0*atan(1.0)
    d2r=pi/180.0

    if (abs(theta_p).le.eps) then
        thp=theta
        php=phi
        return
    endif

    offset = atan(tan(pi-d2r*phi_p)*cos(d2r*theta_p))
    pdiff = d2r*(phi-phi_p)
    costhp = cos(d2r*theta)*cos(d2r*theta_p)+sin(d2r*theta)*sin(d2r*theta_p)*cos(pdiff)

    if (costhp.ge.1.0-eps) then
        thp=0.0
        php=0.0-offset/d2r
        return
    else
        thp=acos(costhp)/d2r
        sinphp=sin(d2r*theta)*sin(pdiff)/sin(d2r*thp)
        cosphp=(cos(d2r*theta)*sin(d2r*theta_p)-sin(d2r*theta)*cos(d2r*theta_p)*cos(pdiff))/sin(d2r*thp)
        if (abs(sinphp).le.eps) then
            if (cosphp.gt.eps) then
                php=0.0-offset/d2r
            else
                php=180.0-offset/d2r
            endif
            return
        else if (abs(cosphp).le.eps) then
            if (sinphp.gt.eps) then
                php=90.0-offset/d2r
            else
                php=-90.0-offset/d2r
            endif
            return
        endif
    endif

    tnum=sin(d2r*theta)*sin(pdiff)
    tden=(cos(d2r*theta)*sin(d2r*theta_p)-sin(d2r*theta)*cos(d2r*theta_p)*cos(pdiff))
    php=(atan2(tnum,tden)-offset)/d2r
    return
end SUBROUTINE cal_pole
