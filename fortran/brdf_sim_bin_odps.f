      program calculate brdf
c     the program is used to conduct BRDF and atmospheric
c     correction using a couple model
c     see Li et al., 2010
c     written by Fuqin Li in 2008 and modified in 2010
c     command line  brdf_sim_bin <brdf_modis> <startend>
c                   <coordinator> <one_image_band>
c                   <solar_angle> <view_angle> <rela_angle>
c                   <fv> <fs> <a> <b> <s>
c                  <output_ref_top> <ref_nobrdf> <ref_wbrdf>
c
      integer nrow,ncol
      real pi
      parameter (pi=3.14159263)
      character fname*150
      integer i,j,n_factor
      real solar_angle(20000),v_angle(20000),rela_angle(20000)
      real fv(20000),fs(20000),a_mod(20000),b_mod(20000),
     & s_mod(20000)
      real ref_top(20000),ref_lm(20000),ref_brdf(20000)
      real brdf0,brdf1,brdf2,bias,slope,esun,dd,lt
      real ann,aa_view,aa_solar,aa_white,solar,
     $ view,ra,norm_1,norm_2
      integer*2 dn(20000),iref_brdf(20000)
      integer*2 iref_lm(20000),iref_top(20000)
      integer line,istart,imid,iend,ii
      integer*1 dn_1(20000)
      double precision a_eq,b_eq,c_eq,ref_bar,aa_final
      real hb,br,pib,fnn
      real RL_brdf,black_sky,white_sky
      external RL_brdf,black_sky,white_sky
      if (IARGC() .ne. 15) then
      write(*,*) 'Error: Required parameters not specified properly!'
      stop 11
      endif

c
      hb=2.0
      br=1.0
      pib=pi/180.0
c
c     input data
c     open input file (BRDF parameter, startend and coordinator)
      do i=1,3
      call GETARG(i, fname)
      open(i,file=fname,status='old')
      enddo

      read(3,*)nrow,ncol

c     open image file
      call GETARG(4, fname)

      open(4,file=fname,status='old',
     !access='direct',form='unformatted',recl=ncol)
c
c     open solar angle, view angle, relative azimuth angle,
c     fv,fs,a,b,s files
c
      do i=1,8
      call GETARG(i+4, fname)
      open(i+10,file=fname,status='old',
     !access='direct',recl=4*ncol)
      enddo
c
c    ouptput data
c    open output for TOA
      call GETARG(13, fname)
      open(52,file=fname, access='direct',recl=2*ncol)
c
c    open output for Lambertian reflectance
      call GETARG(14, fname)
      open(53,file=fname, access='direct',recl=2*ncol)

c     open output for BRDF adjusted surface reflectance
      call GETARG(15, fname)
      open(54,file=fname, access='direct',recl=2*ncol)
c
      read(1,*)brdf0,brdf1,brdf2
      read(1,*)bias,slope,esun,dd
c
      norm_1=brdf1/brdf0
      norm_2=brdf2/brdf0
c      calculate white sky albedo
      aa_white=white_sky(1.0,norm_1,norm_2)
      n_factor=1.0
      if (n_factor.eq.0) then
        fnn=1.0
      else
        fnn=RL_brdf(45*pib,0.0,0.0,hb,br,1.0,norm_1,norm_2)
      endif
      print*,fnn
c
      do i=1,nrow
        read (2,*)line,istart,imid,iend,ii
        read (4,rec=i)(dn_1(j),j=1,ncol)
        read(11,rec=i)(solar_angle(j),j=1,ncol)
        read(12,rec=i)(v_angle(j),j=1,ncol)
        read(13,rec=i)(rela_angle(j),j=1,ncol)
        read(14,rec=i)(fv(j),j=1,ncol)
        read(15,rec=i)(fs(j),j=1,ncol)
        read(16,rec=i)(a_mod(j),j=1,ncol)
        read(17,rec=i)(b_mod(j),j=1,ncol)
        read(18,rec=i)(s_mod(j),j=1,ncol)

         do j=1,ncol
         if (dn_1(j) .lt. 0) then
         dn(j)=dn_1(j)+256
          else
         dn(j)=dn_1(j)
         endif
         enddo

c
        do j=1,ncol

       if (j .ge.istart .and. j .le. iend) then

c     convert angle to radians
            solar=solar_angle(j)*pib
            view=v_angle(j)*pib
            ra=rela_angle(j)*pib
c     calculate radiance at top atmosphere
            lt=bias+dn(j)*slope
c     calculate reflectance at top atmosphere
            ref_top(j)=pi*lt*dd*dd/(esun*cos(solar))
            iref_top(j)=10000*ref_top(j)+0.5
c     calcualte lambetian reflectance with bilinear average
            ref_lm(j)=(lt-b_mod(j))/(a_mod(j)+s_mod(j)*(lt-b_mod(j)))
            iref_lm(j)=ref_lm(j)*10000+0.5
c      set as zero if atmospheric corrected reflectance is zero
c

             if (ref_top(j).lt. 0) then
              ref_top(j)=0.001
              iref_top(j)=10
              endif

            if (ref_lm(j).lt. 0.001) then
             ref_lm(j)=0.001
             ref_brdf(j)=0.001
             iref_lm(j)=10
             iref_brdf(j)=10
           else
c
c     calculate normalized BRDF shape function
            ann=RL_brdf(solar,view,ra,hb,br,1.0,norm_1,norm_2)
c     calculate black sky albedo for sloar angle
            aa_solar=black_sky(1.0,norm_1,norm_2,solar)
c      calculate black sky albedo for view angle
            aa_view=black_sky(1.0,norm_1,norm_2,view)
c
            aa_final=(fv(j)*(fs(j)*ann+(1.0-fs(j))*aa_view)+
     !        (1.0-fv(j))*(fs(j)*aa_solar+(1.0-fs(j))*aa_white))
     !        /aa_white

              a_eq=(1-aa_final)*s_mod(j)
              b_eq=aa_final
              c_eq=-ref_lm(j)

             if (abs(a_eq) .lt. 0.0000001) then
             ref_brdf(j)=ref_lm(j)
             iref_brdf(j)=ref_brdf(j)*10000+0.5
             else
            ref_bar=(-b_eq+sqrt(b_eq**2-4*a_eq*c_eq))/(2*a_eq)
            ref_brdfreal=ann*ref_bar/aa_white

            if (abs(ref_brdfreal-ref_top(j)) .lt. 0.5) then
            ref_brdf(j)=ref_bar*fnn/aa_white
            iref_brdf(j)=ref_brdf(j)*10000+0.5
            else
            ref_brdf(j)=-9.9
            iref_brdf(j)=-999
            endif

            endif
            endif

          if (ref_brdf(j) .ge.1) then
          ref_brdf(j)=1.0
           iref_brdf(j)=ref_brdf(j)*10000
          endif

          else
          ref_top(j)=-999.0
          ref_lm(j)=-999.0
          ref_brdf(j)=-999.0
          iref_brdf(j)=-999
          iref_lm(j)=-999
          iref_top(j)=-999
          endif

        enddo
        write(52,rec=i)(iref_top(j),j=1,ncol)
        write(53,rec=i)(iref_lm(j),j=1,ncol)
        write(54,rec=i)(iref_brdf(j),j=1,ncol)
      enddo
      close(1)
      close(2)
      close(3)
      close(4)
      close(52)
      close(53)
      close(54)
c
      do i=11,18
        close(i)
      enddo
c
      stop
      end
      real function RL_brdf (solar,view,ra,hb,br,brdf0,brdf1,brdf2)
c
      real pi
      parameter (pi=3.14159263)
      real solar, view, ra, hb, br
      real brdf0,brdf1,brdf2
      real rs_thick,li_sparse,secsolar,secvia
      real cossolar,cosvia,cosra,sinsolar,sinvia,sinra
      real cosxi,xi,tansolar,tanvia,theta_new_v,theta_new_s
      real d_li2,x_li,cosl,l_li,o_li
c
      RL_brdf=0.0
c
c     calculate Ross-thick kernel
c
      cossolar=cos(solar)
      cosvia=cos(view)
      cosra=cos(ra)
      sinsolar=sin(solar)
      sinvia=sin(view)
      sinra=sin(ra)
      cosxi=cossolar*cosvia+sinsolar*sinvia*cosra
      if (cosxi .ge.1) then
      cosxi=1
      endif
      xi=acos(cosxi)
      rs_thick=((pi/2.0-xi)*cos(xi)+sin(xi))/(cossolar+cosvia)
     ! -pi/4.0
c
c     calculate Li-sparse
c
      tansolar=sinsolar/cossolar
      tanvia=sinvia/cosvia
      theta_new_v=atan(br*tanvia)
      theta_new_s=atan(br*tansolar)
      cosxi=cos(theta_new_s)*cos(theta_new_v)+sin(theta_new_s)
     ! *sin(theta_new_v)*cosra
       if (cosxi.ge.1) then
       cosxi=1
       endif
      secsolar=1.0/cos(theta_new_s)
      secvia=1.0/cos(theta_new_v)
      d_li2=abs(tan(theta_new_s)**2+tan(theta_new_v)**2-
     ! 2.0*tan(theta_new_s)*tan(theta_new_v)*cosra)
      x_li=tan(theta_new_s)*tan(theta_new_v)*sinra
      cosl=hb*sqrt(d_li2+x_li**2)/(secsolar+secvia)
      if (cosl .ge. 1.0) then
        o_li=0.0
      else
         l_li=acos(cosl)
         o_li=(l_li-sin(l_li)*cos(l_li))*(secsolar+secvia)/pi
      endif
       li_sparse=o_li-(secsolar+secvia)+0.5*(1.0+cosxi)
     !   *secsolar*secvia
      RL_brdf=brdf0+brdf1*rs_thick+brdf2*li_sparse
      return
c
      end
      real function black_sky(brdf0,brdf1,brdf2,theta)
c
      real brdf0,brdf1,brdf2
      real theta
c
c     theta should be in radians
c
      black_sky=brdf0+brdf1*(-0.007574-0.070987*theta**2+
     !        0.307588*theta**3)+brdf2*(-1.284909-0.166314*
     !        theta**2+0.041840*theta**3)
      return
      end
      real function white_sky(brdf0,brdf1,brdf2)
c
      real brdf0,brdf1,brdf2
      white_sky=brdf0+brdf1*0.189184-brdf2*1.377622
      return
      end

