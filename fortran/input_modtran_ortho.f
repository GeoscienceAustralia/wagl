      program generate modtran input file
c     program to generate files needed MODTRAn input files
c     written by Fuqin Li in 2009 and update in 2010 
c     command line  input_modtran_ortho <input_modtraninput> 
c                <coordinator> <view_angle> <azi_angle> 
c               <output_TL_alb_0>  <output_TL_alb_1>
c               <output_TM_alb_0>  <output_TM_alb_1>
c               <output_TR_alb_0>  <output_TR_alb_1>
c               <output_ML_alb_0>  <output_ML_alb_1>
c               <output_MM_alb_0>  <output_MM_alb_1>
c               <output_MR_alb_0>  <output_MR_alb_1>
c               <output_BL_alb_0>  <output_BL_alb_1>
c               <output_BM_alb_0>  <output_BM_alb_1>
c               <output_BR_alb_0>  <output_BR_alb_1>
      character filt*80,name*80,fname*200
      double precision rlat(10),rlon(10),alat(20000),alon(20000)
      real view(20000),azi(20000)
      real view_cor(10),azi_cor(10)
      integer ix(10),iy(10),nrow,ncol,iday
      double precision dres
      double precision rlat0,rlon0,ozone,water,vis,dsea,sht,time
      
            if (IARGC() .ne. 24) then
        write(*,*) 'Error: Required parameters not specified properly!'
         stop 3
         endif
c
c        open modtran input file
         call GETARG(1, fname)
          open(1,file=fname,status='old')
c
c        open coordinator file     
         call GETARG(2, fname)
          open(2,file=fname,status='old')

c     read ozone data (ATM-cm)
      read(1,*)ozone
c     read water vapor data (g/cm2)
      read(1,*)water
c     read sensor response function
      read(1,'(a)')filt
c     read visibility data (km)
      read(1,*)vis
c     read above sea level in km
      read(1,*)dsea
c     annotation
      read(1,'(a)')name
c     satellite altitde above surface in km
      read(1,*)sht
c     day of year
      read(1,*)iday
c     decimal hours in UTC (24 h a day) 
      read(1,*)time
c
c     read row and column
      read(2,*)nrow,ncol
c
c     open satellite view angle file
      call GETARG(3, fname)
      open(11,file=fname,status='old',access='direct',
     ! form='unformatted',recl=4*ncol)
c
c     open satellite azimuth angle file
      call GETARG(4, fname)
      open(12,file=fname,status='old',access='direct',
     ! form='unformatted',recl=4*ncol)

c     open latitude file
      call GETARG(5, fname)
      open(13,file=fname,status='old',access='direct',
     ! form='unformatted',recl=8*ncol)
c     open longtitude file
c
      call GETARG(6, fname)
      open(14,file=fname,status='old',access='direct',
     ! form='unformatted',recl=8*ncol)

c     open output files for 9 coordinators and 2 albedo cases
       do i=1, 9 
c      when albedo is 0
       call GETARG(2*(i-1)+7, fname)
       open (50+i,file=fname)
c      when albedo is 1
       call GETARG(2*(i-1)+8, fname)
       open (60+i,file=fname)
       enddo

       do i=1,9
       read(2,*)ix(i),iy(i)
       read (11,rec=ix(i))(view(j),j=1,ncol)
       read (12,rec=ix(i))(azi(j),j=1,ncol)
       read (13,rec=ix(i))(alat(j),j=1,ncol)
       read (14,rec=ix(i))(alon(j),j=1,ncol)
       view_cor(i)=180-view(iy(i))
       azi_cor(i)=azi(iy(i))+180
       rlat(i)=alat(iy(i))
       rlon(i)=360-alon(iy(i))
c      if it is in western hemisphere
c      -------------------------------
       if (rlon(i) .ge.360) then
       rlon(i)=rlon(i)-360
       endif
       enddo
       do i=1,9
       if (180-view_cor(i) .lt. 0.1) then
       view_cor(i)=180
       azi_cor(i)=0
       endif
       if (azi_cor(i) .gt. 360) then
       azi_cor(i)=azi_cor(i)-360
       endif
       enddo
       do i=1,9
c      generate input file for 0 albedo
       write(50+i,*)0.0
       write(50+i,*)ozone
       write(50+i,*)water
       write(50+i,'(a)')filt
       write(50+i,*)vis
       write(50+i,*)dsea
       write(50+i,'(a)')name
       write(50+i,*)sht
       write(50+i,*) view_cor(i)
       write(50+i,*)iday
       write(50+i,*)rlat(i)
       write(50+i,*)rlon(i)
       write(50+i,*)time
       write(50+i,*)azi_cor(i)
c      generate input file for 1.0 albedo
       write(60+i,*)1.0
       write(60+i,*)ozone
       write(60+i,*)water
       write(60+i,'(a)')filt
       write(60+i,*)vis
       write(60+i,*)dsea
       write(60+i,'(a)')name
       write(60+i,*)sht
       write(60+i,*) view_cor(i)
       write(60+i,*)iday
       write(60+i,*)rlat(i)
       write(60+i,*)rlon(i)
       write(60+i,*)time
       write(60+i,*)azi_cor(i)
       enddo
       close(1)
       close(2)
       close(11)
       close(12)
       do i=1,9
       close(50+i)
       close(60+i)
       enddo
       stop
       end
 

