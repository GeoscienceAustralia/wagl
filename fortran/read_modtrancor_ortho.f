      program coordinator
c       program to find the line and column for 9 coordinators
c       written by Fuqin Li in 2009 and update in 2010
c       command line read_modtrancor_ortho <centreline>
c        <view_angle.bin> <output_coordinator>
        integer line(20000), ncenter(20000),istart(20000),iend(20000)
        real rlat(20000),rlon(20000)
        real view(20000)
        character fname*150
c       program determines 9 coordinators from centre line
c       and view angle
        if (IARGC() .ne. 4) then
          write(*,*) 'Error: Required parameters not specified
     & properly!'
          write(*,*) 'Usage: read_modtrancor_ortho <centreline_file>'
          write(*,*) '                             <view_angle_file>'
          write(*,*) '
     &<output_coordinator_file>'
          write(*,*) '
     &<output_boxline_file>'
          stop 2
        endif
c       open file
c       set the angle threshold for Landsat is 7.69 degree
c       this information can be built at input file.

c       get centreline file
        call GETARG(1, fname)
        open (1, file=fname,status='old')
        read(1,*)viewmax
        read(1,*)nrow,ncol

c       get view angle file
        call GETARG(2, fname)
        open (2, file=fname,status='old',access='direct',
     &    form='unformatted',recl=4*ncol)

c       get output file (coordinator)
        call GETARG(3, fname)
        open (3, file=fname)
        write(3,*)nrow,ncol
        print*,nrow,ncol
c       get output file (start box lines for bilinaer))
        call GETARG(4, fname)
        open (4, file=fname)
        do i=1,nrow
          read(1,*)line(i),ncenter(i),rlat(i),rlon(i)
          read(2,rec=i)(view(j), j=1,ncol)
          kk1=0
          kk2=0
          do j=1,ncol
            if (j .lt. ncenter(i)) then
              if (view(j) .le. viewmax) then
                if (kk1 .eq.0) then
                  istart(i)=j
                  kk1=kk1+1
                endif
              endif
            endif
            if (j .gt. ncenter(i)) then
              if (view(j) .gt. viewmax) then
                if (kk2 .eq.0) then
                  iend(i)=j-1
                  kk2=kk2+1
                endif
              endif
            endif
          enddo
          if (view(1) .lt. viewmax) then
            istart(i)=1
          endif
          if (view(ncol) .le. viewmax) then
            iend(i)=ncol
          endif
          write(4,*)i,istart(i),iend(i)
        enddo

        write(3,*)line(1),istart(1)
        write(3,*)line(1),ncenter(1)
        write(3,*)line(1),iend(1)
        write(3,*)line(nrow/2),istart(nrow/2)
        write(3,*)line(nrow/2),ncenter(nrow/2)
        write(3,*)line(nrow/2),iend(nrow/2)
        write(3,*)line(nrow),istart(nrow)
        write(3,*)line(nrow),ncenter(nrow)
        write(3,*)line(nrow),iend(nrow)
        close(1)
        close(2)
        close(3)
        close(4)
        stop
      end






