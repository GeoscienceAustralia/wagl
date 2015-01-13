        program calculate coefficients
c      program to computer the atmospheric parameters needed by
c      BRDF and atmospheric correction couple model
c      written by Fuqin Li in 2009 and modified 2010
c     command line coefficient <filter_input> <*_0.chn>
c                                 <*_1.chn>
c                      <*_0.dir> <*_1.dir> <*_t.dir>
c                             <output_*_alb.txt>

c      the program is used to calculate coefficients
c      which are used in teh BRDF and atmospheric correction
       parameter (pi=3.14159263)
       character header*150,fname*150,fname1*150

      if (IARGC() .ne. 7) then
      write(*,*) 'Error: Required parameters not specified properly!'
      stop 8
      endif
c     open sensor filter information
      call GETARG(1, fname1)
      open(2,file=fname1,status='old')
c     open MODTRAN output .chn file (albedo 0)
      call GETARG(2, fname1)
      open(11,file=fname1,status='old')
c     open MODTRAN output .chn file (albedo 1)
      call GETARG(3, fname1)
      open(13,file=fname1)
c     open solar radiation file (albedo 0)
      call GETARG(4, fname1)
      open(14,file=fname1)
c     open solar radiation file (albedo 1)
      call GETARG(5, fname1)
      open(15,file=fname1)
c     open open solar radiation file (transmittance mode)
      call GETARG(6, fname1)
      open(16,file=fname1)
c     open output file
      call GETARG(7, fname1)
      open(51,file=fname1)

c     read satellite information, eg. bands and filter function

      read(2,*)nband
      read(2,'(a)')fname

c      read header
       do i=1,4
       read(11,*)header
       read(13,*)header
       enddo
       read(14,*)header
       read(15,*)header
       read(16,*)header
       write(51,*)'band    ','fs     ','fv    ','a    ','b    ','s     '
     #,'direct  ','diffuse  ','ts  '
       do i=1,nband
       read(11,*)aa,iband,aa,path_0
       read(13,*)aa,iband,aa,path_1
       path_0=path_0*10000000
       path_1=path_1*10000000
       read(14,*)iband,diff_0,dir_0,dir0_top
       read(15,*)iband,diff_1,dir_1,dir1_top
       read(16,*)iband,diff_t,dir_t,aa,dirt_top,tv_total
       b=path_0
       diff_0=diff_0*10000000.0
       dir_0=dir_0*10000000.0
       dir0_top=dir0_top*10000000.0
       diff_1=diff_1*10000000.0
       dir_1=dir_1*10000000.0
       dir1_top=dir1_top*10000000.0
       ts_total=(diff_0+dir_0)/dir0_top
       ts_dir=dir_0/dir0_top
       fs=ts_dir/ts_total
       tv_dir=dir_t/dirt_top
       fv=tv_dir/tv_total
       s=1-(diff_0+dir_0)/(diff_1+dir_1)
       a=(diff_0+dir_0)/pi*tv_total
       write(51,*)iband,fs,fv,a,b,s,dir_0,diff_0,ts_dir
       enddo
       close(2)
       do i=1,5
       close(i+10)
       enddo
       close(51)
       stop
       end
