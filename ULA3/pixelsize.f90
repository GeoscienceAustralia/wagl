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
