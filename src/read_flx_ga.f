      program to extract solar irradiance
c       program is used to extract information from MODTRAN
c       output file *.flx.
c       program is written by Fuqin Li in 2002 and
c       modified by Fuqin Li in 2008, 2009, 2010
c       command line read_flx_ga <*.flx> <filter_input>
c                               <output_*.dir>

        character a*150,fname*250,fname1*150
        real hgt(90),wave(2600),eup(2600,90),edo(2600,90)
        real esun(2600,90),res(2600,20),wave1(2600,20)
        real edo_band(20),esun_band(20),res_sum(20),esun_band_top(20)

        if (IARGC() .ne. 3) then
          write(*,*) 'Error: Required parameters not specified
     & properly!'
          stop 6
        endif
c       open modtran output *.flx file
        call GETARG(1, fname1)
        open(1,file=fname1,status='old')
c       open filter file
        call GETARG(2, fname1)
        open(2,file=fname1,status='old')
c       open output file
        call GETARG(3, fname1)
        open(9,file=fname1)
c
c       read satellite information, eg. bands and filter function

        read(2,*)nband
        read(2,'(a)')fname

c       read header from flux data
        do i=1,6
          read(1,'(a)')a
        enddo
c       read the profile layer number
c       assume that the profile layer number is less than 100
        read(1,'(a)')a
        read(a(2:3),'(i2)')inu
c       read header again
        do i=1,4
          read(1,'(a)')a
        enddo
        ii=inu/2
c       read profile height
        do i=1,ii
          read(1,'(a)')a
          read(a(20:28),'(f9.5)')hgt(i*2-1)
          read(a(56:64),'(f9.5)')hgt(i*2)
        enddo
        if((inu/2*2) .ne.inu) then
          read(1,'(a)')a
          read(a(20:28),'(f9.5)')hgt(ii*2+1)
        endif
        do i=1,4
          read(1,'(a)')a
        enddo
c       read modtran output for upward, downward, direct solar radiation
c       assume that teh wavelength is between 350 and 2600 in your
c        modtran input
        do k=2600,350,-1
          wave(k)=float(k)
          do i=1,ii
            read(1,'(a)')a
            read(a(9:20),'(e12.5)')eup(k,i*2-1)
            read(a(21:32),'(e12.5)')edo(k,i*2-1)
            read(a(33:44),'(e12.5)')esun(k,i*2-1)

c            read(a(9:44),'(3e12.5)')(eup(k,i*2-1),
c     &        edo(k,i*2-1),esun(k,i*2-1))

            read(a(45:56),'(e12.5)')eup(k,i*2)
            read(a(57:68),'(e12.5)')edo(k,i*2)
            read(a(69:80),'(e12.5)')esun(k,i*2)

c            read(a(45:80),,'(3e12.5))(eup(k,i*2),
c     &        edo(k,i*2),esun(k,i*2))
          enddo
          if(inu/2*2 .ne.inu) then
            read(1,'(a)')a
            read(a(9:20),'(e12.5)')eup(k,ii*2+1)
            read(a(21:32),'(e12.5)')edo(k,ii*2+1)
            read(a(33:44),'(e12.5)')esun(k,ii*2+1)

c            read(a(9:44),'(3e12.5)')(eup(k,ii*2+1),edo(k,ii*2+1)
c     &        ,esun(k,ii*2+1))
          endif
        enddo
c        do i=350,2600
c          write(13,*)wave(i),edo(i,1),esun(i,1)
c        enddo
c        open filter function. Note: increasement of filter function
c        must be the same.

        open(4,file=fname,status='old')
        do i=2600,350,-1
          do j=1,nband
            res(i,j)=0.0
          enddo
        enddo

        kk=0
        read(4,'(a)')a
10      read(4,'(a)',end=20,err=20)a
        if(a(1:1) .eq.'B') then
          kk=kk+1
          go to 10
        else
          read(a(3:6),'(i4)')jj
          read(a(1:17),'(f8.1,f9.5)')wave1(jj,kk),res(jj,kk)
c          write(12,*)jj,kk,wave1(jj,kk),res(jj,kk)
        endif
        go to 10
c       test data
c       landsat 7 has 7 bands, land5 has 6 bands
20      do i=1,nband
          do j=2600,350,-1
c            write(11,*)i,j,res(j,i)
          enddo
        enddo
        write(9,*)'band   ','diffuse   ','direct   ','directtop   '
        do j=1,nband
          edo_band(j)=(res(2600,j)*edo(2600,1)+res(350,j)*edo(350,1))/2
          esun_band(j)=(res(2600,j)*esun(2600,1)+res(350,j)*
     &      esun(350,1))/2
          esun_band_top(j)=(res(2600,j)*esun(2600,inu)+res(350,j)
     &      *esun(350,inu))/2
          res_sum(j)=(res(2600,j)+res(350,j))/2
          do i=2599,351,-1
            edo_band(j)=edo_band(j)+edo(i,1)*res(i,j)
            esun_band(j)=esun_band(j)+esun(i,1)*res(i,j)
            esun_band_top(j)=esun_band_top(j)+esun(i,inu)*res(i,j)
            res_sum(j)=res_sum(j)+res(i,j)
          enddo
          edo_band(j)=edo_band(j)/res_sum(j)
          esun_band(j)=esun_band(j)/res_sum(j)
          esun_band_top(j)=esun_band_top(j)/res_sum(j)
          write(9,*)j, edo_band(j),esun_band(j),esun_band_top(j)
        enddo
        stop
      end
