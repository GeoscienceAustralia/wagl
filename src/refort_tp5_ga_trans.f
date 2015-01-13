      program reformat
c       command line refort_tp5_trans <*_alb_*> <default.tp5>
c                               <output_*.tp5>

c       this program is used to modtran input file for tranmissant mode
c       written by Fuqin Li in 2009 and update in 2010
c       the input file does not use profile, but total water vapor
        character a(150)*110,filt*80,name*80
        character spa*10,fname*150
c        spa='          '
        spa=' '
        if (IARGC() .ne. 3) then
          write(*,*) 'Error: Required parameters not specified
     & properly!'
          stop 5
        endif

c       open modtran input file
        call GETARG(1, fname)
        open(1,file=fname,status='old')

c       open default .tp5 file
        call GETARG(2, fname)
        open(2,file=fname,status='old')

c       open output file
        call GETARG(3, fname)
        open(3,file=fname)

c       read albedo
        read(1,*)albedo
c       read input ozone data in ATM-cm
        read(1,*)ozone
c       read input  water vapor in g/cm2
        read(1,*)water
c       read instrument filt function
        read(1,'(a)')filt
c       read optical depth (negative value of
c       original optical depth at 550 nm)
        read(1,*)vis
c       read above sea level on teh site in km
        read(1,*)dsea
c       explantion of teh file
        read(1,'(a)')name
c       satellite height in km
        read(1,*)sht
c       view zenith angle in degree
        read(1,*)view
c       day of year
        read(1,*)iday
c       latitude in degree follow modtran
        read(1,*)rlat
c       longitude in degree follow modtran
        read(1,*)rlon
c       time decimal hour in UTC
        read(1,*)time
c       azimuth angle in degree
        read(1,*)azimuth

        ii=1
400     read(2,'(a)',end=300,err=300)a(ii)
        ii=ii+1
        go to 400
300     write(3,'(a73,f7.2)')a(1)(1:73),albedo
        write(3,'(a22,a1,f7.5,a5,f5.3,a70)')a(2)(1:22),'g',water,
     &    a(2)(31:35),ozone,a(2)(41:110)
        write(3,'(a)')filt
        write(3,'(a30,f10.5,a30,f10.3)')a(4)(1:30),vis,a(4)(41:70),dsea
        write(3,'(3f10.3,a50)')sht,dsea,view,a(5)(31:80)
        write(3,'(4i5)')2,0,iday,0
        write(3,'(2f10.3,5a10,f10.3)')0.0,180.0-view,spa,
     &    spa,spa,spa,spa,0.667
        write(3,'(a)')a(8)
        write(3,'(a)')a(9)
        close(1)
        close(2)
        close(3)
        stop
      end
