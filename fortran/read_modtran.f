      program read output
c     program reformat the previous atmospheric parameters to
c     parameters for four boxes. These aaaare needed to conduct
c     bilinear analysis.
c     written by Fuqin Li in 2009 and modified in 2010
c     command line read_modtran <filter_input> <tl.txt> <tm.txt>
c     <tr.txt> <ml.txt> <mm.txt> <mr.txt> <bl.txt> <bm.txt> <br.txt>
c     <output_*fv_b1> <fs_b1> <b_b1> <s_b1> <a_b1> <dir_b1> <dif_b1> <ts_b1>
c     <output_*fv_b2> <fs_b2> <b_b2> <s_b2> <a_b2> <dir_b1> <dif_b1> <ts_b1>
c     <output_*fv_b3> <fs_b3> <b_b3> <s_b3> <a_b3> <dir_b1> <dif_b1> <ts_b1>
c     <output_*fv_b4> <fs_b4> <b_b4> <s_b4> <a_b4> <dir_b1> <dif_b1> <ts_b1>
c     <output_*fv_b5> <fs_b5> <b_b5> <s_b5> <a_b5> <dir_b1> <dif_b1> <ts_b1>
c     <output_*fv_b7> <fs_b7> <b_b7> <s_b7> <a_b7> <dir_b1> <dif_b1> <ts_b1>

      character header*100,fname*50,fname1*150
      real x,fv(20,9),fs(20,9),b(20,9),s(20,9)
      real a(20,9),dir(20,9),dif(20,9),ts(20,9)

c      if (IARGC() .ne. 58) then
c      write(*,*) 'Error: Required parameters not specified properly!'
c      stop 9
c      endif

c
c    open sensor filter file
      call GETARG(1, fname1)
      open(2,file=fname1,status='old')
c
c     read satellite information, eg. bands and filter function

      read(2,*)nband
      read(2,*)fname
c
c     open input files from previous run
c     atmospheric parameters for 9 coordinators
      do i=1,9
      call GETARG(i+1, fname1)
      open(10+i,file=fname1,status='old')
      enddo
c
c     open output file for each atmospheric parameters and
c     individual bands
c
      do i=1,nband
      do j=1,8
      call GETARG((i-1)*8+j+10, fname1)
      open(i*10+j+10,file=fname1)
      enddo
      enddo
c     read modtran generated output
      do j=1,9
      read(10+j,*)header
      enddo
      do i=1,nband
      do j=1,9
      read(10+j,*)iband,fs(i,j),fv(i,j),a(i,j),b(i,j),s(i,j),
     # dir(i,j),dif(i,j),ts(i,j)
      enddo
      enddo
c
      do i=1,nband
      write(i*10+11,*)fv(i,1),fv(i,2),fv(i,4),fv(i,5)
      write(i*10+11,*)fv(i,2),fv(i,3),fv(i,5),fv(i,6)
      write(i*10+11,*)fv(i,4),fv(i,5),fv(i,7),fv(i,8)
      write(i*10+11,*)fv(i,5),fv(i,6),fv(i,8),fv(i,9)
c
      write(i*10+12,*)fs(i,1),fs(i,2),fs(i,4),fs(i,5)
      write(i*10+12,*)fs(i,2),fs(i,3),fs(i,5),fs(i,6)
      write(i*10+12,*)fs(i,4),fs(i,5),fs(i,7),fs(i,8)
      write(i*10+12,*)fs(i,5),fs(i,6),fs(i,8),fs(i,9)
c
      write(i*10+13,*)b(i,1),b(i,2),b(i,4),b(i,5)
      write(i*10+13,*)b(i,2),b(i,3),b(i,5),b(i,6)
      write(i*10+13,*)b(i,4),b(i,5),b(i,7),b(i,8)
      write(i*10+13,*)b(i,5),b(i,6),b(i,8),b(i,9)
c
      write(i*10+14,*)s(i,1),s(i,2),s(i,4),s(i,5)
      write(i*10+14,*)s(i,2),s(i,3),s(i,5),s(i,6)
      write(i*10+14,*)s(i,4),s(i,5),s(i,7),s(i,8)
      write(i*10+14,*)s(i,5),s(i,6),s(i,8),s(i,9)
c
      write(i*10+15,*)a(i,1),a(i,2),a(i,4),a(i,5)
      write(i*10+15,*)a(i,2),a(i,3),a(i,5),a(i,6)
      write(i*10+15,*)a(i,4),a(i,5),a(i,7),a(i,8)
      write(i*10+15,*)a(i,5),a(i,6),a(i,8),a(i,9)
c
      write(i*10+16,*)dir(i,1),dir(i,2),dir(i,4),dir(i,5)
      write(i*10+16,*)dir(i,2),dir(i,3),dir(i,5),dir(i,6)
      write(i*10+16,*)dir(i,4),dir(i,5),dir(i,7),dir(i,8)
      write(i*10+16,*)dir(i,5),dir(i,6),dir(i,8),dir(i,9)
c
      write(i*10+17,*)dif(i,1),dif(i,2),dif(i,4),dif(i,5)
      write(i*10+17,*)dif(i,2),dif(i,3),dif(i,5),dif(i,6)
      write(i*10+17,*)dif(i,4),dif(i,5),dif(i,7),dif(i,8)
      write(i*10+17,*)dif(i,5),dif(i,6),dif(i,8),dif(i,9)
c
      write(i*10+18,*)ts(i,1),ts(i,2),ts(i,4),ts(i,5)
      write(i*10+18,*)ts(i,2),ts(i,3),ts(i,5),ts(i,6)
      write(i*10+18,*)ts(i,4),ts(i,5),ts(i,7),ts(i,8)
      write(i*10+18,*)ts(i,5),ts(i,6),ts(i,8),ts(i,9)


      enddo
      close(2)
      do i=1,9
      close(10+i)
      enddo
      do i=1,nband
      do j=1,8
      close(i*10+j+10)
      enddo
      enddo

      stop
      end
