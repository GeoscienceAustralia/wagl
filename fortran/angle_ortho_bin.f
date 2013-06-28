      program angle
c     command line angle_ortho_bin <header_angle> <image_file_one_band>
c               <output_startend> <output_centreline> <view_angle>
c               <azi_angle> <solar_angle> <sazi_angle> <rela_angle>
c     program to calculate solar, view and azimuth angle from lat
c     lon projection if we only know the center point in the
c      original image.
c     In this calculation, we assumed that satellite track in the
c     scene as one straight line. Then we can calculate the view and
c     azimuth angle using sphere geometry
c     the program is written by Fuqin Li, May, 2008.
c-------------------------------------------------------------
      double precision alon(20000)
      real view(20000),azi(20000)
      real asol(20000),phi0(20000),rela_angle(20000)
      double precision beta(20000)
c      double precision skew(20000)
      double precision alat1,alon1,alat0,alon0,rr,oi,v0
      double precision aa,ee,w0,rnno,rlat0,rho0,time0,gamma0
      double precision rnn1,rlat1,rho1,time1,alon1c
      double precision alatn,rnnn,rlatn,rhon,timen,alonnc
      double precision alpha,alat,rnncx,rlatcx,rhocx,timecx
      double precision aloncx,tu,dres,pi,pia,pib,viewmax
      double precision dlon,dlat,brg,rnn,rlat,yy,dclat,dclon
      double precision alonc,rlatc,hh,thetap
      integer nrow,ncol,iday,line,istart,imid,iend,ii
      integer*2 dn(20000)
      integer*1 dn_1(20000)
      character fname*256
      pi=4.0*atan(1.0)
      pia=pi/180.0
      pib=180.0/pi
      if (IARGC() .ne. 9) then
        write(*,*) 'Error: Required parameters not specified properly!'
        write(*,*) 'Usage: angle_ortho_bin <header_angle>'
        write(*,*) '                       <image_file_one_band>'
        write(*,*) '                       <output_startend>'
        write(*,*) '                       <output_centreline>'
        write(*,*) '                       <output_view_angle>'
        write(*,*) '                       <output__azi_angle>'
        write(*,*) '                       <output_solar_angle>'
        write(*,*) '                       <output_sazi_angle>'
        write(*,*) '                       <output_rela_angle>'

        stop 1
      endif
c     open header file
c
      call GETARG(1, fname)
      open(1,file=fname,status='old')
c
c     read day and time (GMT)
      read(1,*)iday,tu
c     read row and column
      read(1,*)nrow,ncol
c     read image spatial resolution
      read(1,*)dres
c     read upperleft lat and lon
      read(1,*)alat1,alon1
c     read centre point lat and lon
      read(1,*)alat0,alon0
c     read satellite element, semi-major radius,orbital inclination,
c     and angular velocity
      read(1,*)rr,oi,v0
c     read maximum view angle
      read(1,*)viewmax
c
c    open raw image file
      call GETARG(2, fname)
       open(2,file=fname,status='old',access='direct',
     ! form='unformatted',recl=ncol)

c     output
c
c     open startend file
        call GETARG(3, fname)
        open(3,file=fname)
c
c      open centerline file
        call GETARG(4, fname)
         open(55,file=fname)
          write(55,*)viewmax
          write(55,*)nrow,ncol
c
c       open satellite view angle
        call GETARG(5, fname)
      open (11,file=fname,form='unformatted',
     ! access ='direct',recl=4*ncol)
c
c     open satellite azimuth angle file
      call GETARG(6, fname)
      open (12,file=fname,form='unformatted',
     ! access ='direct',recl=4*ncol)
c
c     open solar zenith angle file
      call GETARG(7, fname)
      open (13,file=fname,form='unformatted',
     ! access ='direct',recl=4*ncol)
c
c     open solar azimuth file
      call GETARG(8, fname)
      open (14,file=fname,form='unformatted',
     ! access ='direct',recl=4*ncol)
c
c    open relative azimuth angle
      call GETARG(9, fname)
      open (15,file=fname,form='unformatted',
     ! access ='direct',recl=4*ncol)

c     convert to radians
      alat1=alat1*pia
      alon1=alon1*pia
      alat0=alat0*pia
      alon0=alon0*pia
      oi=oi*pia
c
c     calculate satellite heading and assumed that it is a line
c     through the center point
c     spherod semi-major axis(m) (earth)
      aa= 6378137d0
c     shperoid eccentricity e2 (earth)
      ee=0.00669438d0
c     earth rotation w0 radians/second
      w0=0.000072722052d0
c     calculate geocentric latitude at center point
      rnn0=aa/sqrt(1-ee*sin(alat0)**2)
      rlat0=alat0-asin(rnn0*ee*sin(alat0)*cos(alat0)/rr)
c     calculate rho at centre point
      rho0=acos(sin(rlat0)/sin(oi))
c     calculate time at centre point of the image
       time0=(rho0+pi/2.0)/v0
c     calculate the longitude of the ascending node(in radians)
      gamma0=alon0+pi/2-atan(tan(rho0)/cos(oi))

c     calculate the centre point at first and last line of image

c     for the first line
      rnn1=aa/sqrt(1-ee*sin(alat1)**2)
      rlat1=alat1-asin(rnn1*ee*sin(alat1)*cos(alat1)/rr)
       rho1=acos(sin(rlat1)/sin(oi))
      time1=(rho1+pi/2)/v0
      alon1c=gamma0-pi/2+atan(tan(rho1)/cos(oi))-w0*(time1-time0)

c     for the last line
      alatn=alat1-(nrow-1)*dres*pia
      rnnn=aa/sqrt(1-ee*sin(alatn)**2)
      rlatn=alatn-asin(rnnn*ee*sin(alatn)*cos(alatn)/rr)
      rhon=acos(sin(rlatn)/sin(oi))
      timen=(rhon+pi/2)/v0
      alonnc=gamma0-pi/2+atan(tan(rhon)/cos(oi))-w0*(timen-time0)

c     calcualte centre line (alpha = heading+skew)

      alpha=atan(tan(alon1c-alonnc)/sin(rlat1-rlatn))
c      print*,alpha*pib

c     for simplify, we only use one hh and rnn for the whole scene
c     when calculate final view angle. the error is very small and the
c     calculation is much faster. We used the centre point.

       hh=rr*cos(rlat0)/cos(alat0)-rnn0
c       print*,alat0,rnn0,rlat0,hh

c     read raw image header

c     read raw digital number to decide the center point on the line
      do i=1,nrow
      read(2,rec=i)(dn_1(j),j=1,ncol)
      do j=1,ncol
      if (dn_1(j) .lt. 0) then
      dn(j)=dn_1(j)+256
      else
      dn(j)=dn_1(j)
      endif
      enddo
c
c     calculate startend pixels.
      istart=1
      iend=0
      ii=0
      do j=2,ncol
      if (dn(1) .le. 0) then
      if (dn(j) .gt. 0 .and. dn(j-1) .le.0) then
      if (istart.le.1) then
      istart=j
      else
      istart=istart
      endif
      endif
      else
      istart=1
      endif
      if (dn(j) .gt. 0 ) then
      iend=j
      endif
      enddo
      ii=iend-istart+1
      imid=ii/2+istart
      write(3,*)i,istart,imid,iend,ii

c     initial value
      do j=1,ncol
      view(j)=-999
      azi(j)=-999
      asol(j)=-999
      phi0(j)=-999
      rela_angle(j)=-999
      enddo


c     calculate centre line for each point as test

      alat=alat1-(i-1)*dres*pia
      rnncx=aa/sqrt(1-ee*sin(alat)**2)
      rlatcx=alat-asin(rnncx*ee*sin(alat)*cos(alat)/rr)
      rhocx=acos(sin(rlatcx)/sin(oi))
      timecx=(rhocx+pi/2)/v0
      aloncx=gamma0-pi/2+atan(tan(rhocx)/cos(oi))-w0*(timecx-time0)
      kk=(aloncx-alon1)/dres*pib+1
      write(55,*)i,kk,alat*pib,aloncx*pib

c      print*,i,alat*pib,aloncx*pib
      do j=1,ncol

      alon(j)=alon1+(j-1)*dres*pia

c     calculate solar angle
      call solar (alat,alon(j),iday,tu,asol(j),phi0(j))
c
      rnn=aa/sqrt(1-ee*sin(alat)**2)
      rlat=alat-asin(rnn*ee*sin(alat)*cos(alat)/rr)
c
c     calculate lat and lon at centre line correponspond to
c      the each grid point.
c      first, calculate arc between the anypoint at image to
c      centre point.

      dlon=alon(j)-alon0

      dlat=rlat-rlat0
      darc=acos(cos(dlon)*cos(dlat))
c     calculate bearing angle for each point
c     the procedure pretty follow slates, but in
c     lat, lon system using sphere geometry

           if (dlat.eq.0) then
             if(dlon.ge.0) brg =pi/2
               if(dlon.lt.0) brg =pi/2*3
                 else
                   if(dlat.gt.0) then
                     brg = atan(tan(dlon)/sin(dlat))
                     if(brg.lt.0) brg = brg + 2*pi
                     else
                    brg = atan(tan(dlon)/sin(dlat)) + pi
                     endif
                   endif

       if (abs(brg-alpha) .le. 0.0001 .or. abs(brg-alpha-pi)
     ! .le.0.0001)   then
       view(j)=0
       azi(j)=0
       rela_angle(j)=azi(j)-phi0(j)
       write(56,*)i,j,alat*pib,alon(j)*pib
c       print*,i,j,view(j),azi(j)
       endif

c      first section
       if (brg-alpha .gt.  0.0001 .and. brg-alpha .lt. pi/2) then
       yy=atan(tan(darc)*cos(brg-alpha))
       dclat=atan(tan(yy)*cos(alpha))
       dclon=asin(sin(yy)*sin(alpha))
       rlatc=rlat0+dclat
       alonc=alon0+dclon

       thetap=acos(sin(rlatc)*sin(rlat)+cos(rlatc)*
     ! cos(rlat)*cos(alon(j)-alonc))
       beta(j)=atan(-1/(tan(oi)*sqrt(1-(sin(rlatc)/sin(oi))**2)))
c       skew(j)=w0*cos(rlatc)*cos(beta(j))/(v0+w0*cos(rlatc)
c     !*sin(beta(j)))
        view(j)=rnn0*thetap/hh*pib
        azi(j)=beta(j)*pib+270
        rela_angle(j)=azi(j)-phi0(j)
        endif

c      second section
       if (brg-alpha .ge. pi/2 .and.
     ! brg-alpha .lt. pi-0.0001) then
       yy=atan(tan(darc)*cos(pi-brg+alpha))
       dclat=atan(tan(yy)*cos(alpha))
       dclon=asin(sin(yy)*sin(alpha))
       rlatc=rlat0-dclat
       alonc=alon0-dclon

       thetap=acos(sin(rlatc)*sin(rlat)+cos(rlatc)*
     ! cos(rlat)*cos(alon(j)-alonc))
       beta(j)=atan(-1/(tan(oi)*sqrt(1-(sin(rlatc)/sin(oi))**2)))
c       skew(j)=w0*cos(rlatc)*cos(beta(j))/(v0+w0*cos(rlatc)
c     ! *sin(beta(j)))
        view(j)=rnn0*thetap/hh*pib
        azi(j)=beta(j)*pib+270
        rela_angle(j)=azi(j)-phi0(j)
       endif

c      third section

       if (brg-alpha .gt. pi+0.0001 .and.
     ! brg-alpha .lt. pi*3/2) then
       yy=atan(tan(darc)*cos(brg-pi-alpha))
       dclat=atan(tan(yy)*cos(alpha))
       dclon=asin(sin(yy)*sin(alpha))
       rlatc=rlat0-dclat
       alonc=alon0-dclon

       thetap=acos(sin(rlatc)*sin(rlat)+cos(rlatc)*
     ! cos(rlat)*cos(alon(j)-alonc))
       beta(j)=atan(-1/(tan(oi)*sqrt(1-(sin(rlatc)/sin(oi))**2)))
c       skew(j)=w0*cos(rlatc)*cos(beta(j))/(v0+w0*cos(rlatc)
c     ! *sin(beta(j)))
        view(j)=rnn0*thetap/hh*pib
        azi(j)=beta(j)*pib+90.0
        rela_angle(j)=azi(j)-phi0(j)
        endif

c      forth section
       if (brg-alpha .ge. pi*3/2 .and.
     ! brg .lt.2*pi) then

       yy=atan(tan(darc)*cos((2*pi-brg+alpha)))
       dclat=atan(tan(yy)*cos(alpha))
       dclon=asin(sin(yy)*sin(alpha))
       rlatc=rlat0+dclat
       alonc=alon0+dclon

       thetap=acos(sin(rlatc)*sin(rlat)+cos(rlatc)*
     ! cos(rlat)*cos(alon(j)-alonc))
       beta(j)=atan(-1/(tan(oi)*sqrt(1-(sin(rlatc)/sin(oi))**2)))
c       skew(j)=w0*cos(rlatc)*cos(beta(j))/(v0+w0*cos(rlatc)
c     ! *sin(beta(j)))
        view(j)=rnn*thetap/hh*pib
        azi(j)=beta(j)*pib+90.0
        rela_angle(j)=azi(j)-phi0(j)
        endif

c      fifth section
       if (brg.ge.0 .and. brg .lt.alpha-0.00001) then
        yy=atan(tan(darc)*cos((alpha-brg)))
       dclat=atan(tan(yy)*cos(alpha))
       dclon=asin(sin(yy)*sin(alpha))
       rlatc=rlat0+dclat
       alonc=alon0+dclon

       thetap=acos(sin(rlatc)*sin(rlat)+cos(rlatc)*
     ! cos(rlat)*cos(alon(j)-alonc))
       beta(j)=atan(-1/(tan(oi)*sqrt(1-(sin(rlatc)/sin(oi))**2)))
c       skew(j)=w0*cos(rlatc)*cos(beta(j))/(v0+w0*cos(rlatc)
c     ! *sin(beta(j)))
        view(j)=rnn*thetap/hh*pib
        azi(j)=beta(j)*pib+90.0
        rela_angle(j)=azi(j)-phi0(j)
        endif
        enddo
       write(11,rec=i)(view(j),j=1,ncol)
       write(12,rec=i)(azi(j),j=1,ncol)
       write(13,rec=i)(asol(j),j=1,ncol)
       write(14,rec=i)(phi0(j),j=1,ncol)
       write(15,rec=i)(rela_angle(j),j=1,ncol)
        enddo
       close(1)
       close(2)
       close(3)
       close(55)
       close(11)
       close(12)
       close(13)
       close(14)
       close(15)
       stop
       end

      subroutine solar(xlat,xlon,iday,tu,asol,phi0)

c    xlat-latitude,xlon, longitude,iday, day of year,
c    xu GMT time, asol, solar zenith angle,phi0,solar azimuth
c    the code original is from Modtran and
c    modified  by Fuqin Li, 2005

      double precision   tu, xlat, tsm, xlon,xla, xj, tet,
     !   a1, a2, a3, a4, a5, et, tsv, ah, b1, b2, b3, b4,
     !   b5, b6, b7, delta, amuzero, elev, caz, azim, pi2
      real asol,phi0,az

      double precision pi,fac
      integer iday
      pi=4.0*atan(1.0)
      fac=pi/180.0

c     solar position (zenithal angle asol,azimuthal angle phi0
c                     in degrees)
c     j is the day number in the year
c

c    mean solar time (heure decimale)

      tsm=tu+xlon*180/pi/15
      xla=xlat
      xj=float(iday)
      tet=2.*pi*(xj)/365.

c    time equation (in mn.dec)
      a1=.000075d0
      a2=.001868d0
      a3=.032077d0
      a4=.014615d0
      a5=.040849d0
      et=a1+a2*cos(tet)-a3*sin(tet)-a4*cos(2.*tet)-a5*sin(2.*tet)
      et=et*12.*60./pi

c     true solar time

      tsv=tsm+et/60.
      tsv=(tsv-12.)

c     hour angle
      ah=tsv*15.*fac

c     solar declination   (in radian)

      b1=.006918d0
      b2=.399912d0
      b3=.070257d0
      b4=.006758d0
      b5=.000907d0
      b6=.002697d0
      b7=.001480d0
      delta=b1-b2*cos(tet)+b3*sin(tet)-b4*cos(2.*tet)+b5*sin(2.*tet)-
     &b6*cos(3.*tet)+b7*sin(3.*tet)

c     elevation,azimuth
      amuzero=sin(xla)*sin(delta)+cos(xla)*cos(delta)*cos(ah)
      elev=asin(amuzero)
      az=cos(delta)*sin(ah)/cos(elev)
      if ( (abs(az)-1.).gt.0.0) az = sign(1.00,az)
      caz=(-cos(xla)*sin(delta)+sin(xla)*cos(delta)*cos(ah))/cos(elev)
      azim=asin(az)
      if(caz.le.0.) azim=pi-azim
      if(caz.gt.0.and.az.le.0) azim=2*pi+azim
      azim=azim+pi
      pi2=2*pi
      if(azim.gt.pi2) azim=azim-pi2
      elev=elev*180./pi
c     conversion in degrees

      asol=90.-elev
      phi0=azim/fac
c      print*,elev,asol,phi0,cos(asol/180*pi)
      return
      end
