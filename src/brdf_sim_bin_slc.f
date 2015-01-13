      program CALCULATE_BRDF_SLC

c     command line  brdf_sim_bin <brdf_modis> <startend>
c                   <coordinator> <one_image_band>
c                   <solar_angle> <view_angle> <rela_angle>
c                   <fv> <fs> <a> <b> <s>
c                  <output_ref_top> <ref_nobrdf> <ref_wbrdf>
c

!DEC$ ATTRIBUTES FORCEINLINE :: black_sky
!DEC$ ATTRIBUTES FORCEINLINE :: white_sky
!DEC$ ATTRIBUTES FORCEINLINE :: RL_brdf

      INTEGER, PARAMETER :: NCDIM = 16000
      REAL, PARAMETER    :: pi = 3.14159263

      CHARACTER fname*150

      INTEGER nrow, ncol
      INTEGER i, j, n_factor
      INTEGER   line, istart, imid, iend, ii

      REAL solar_angle(NCDIM), v_angle(NCDIM), rela_angle(NCDIM)
      REAL fv(NCDIM), fs(NCDIM)
      REAL a_mod(NCDIM), b_mod(NCDIM), s_mod(NCDIM)
      REAL ref_top(NCDIM), ref_lm(NCDIM), ref_brdf(NCDIM)
      REAL brdf0, brdf1, brdf2
      REAL bias, slope, esun, dd, lt
      REAL ann, aa_view, aa_solar, aa_white, solar
      REAL view, ra, norm_1, norm_2
      REAL hb, br, pib, fnn

      INTEGER*2 dn(NCDIM), iref_brdf(NCDIM)
      INTEGER*2 iref_lm(NCDIM), iref_top(NCDIM)
      INTEGER*1 dn_1(NCDIM)

      DOUBLE PRECISION a_eq, b_eq, c_eq, ref_bar, aa_final

      REAL     RL_brdf, black_sky, white_sky
      EXTERNAL RL_brdf, black_sky, white_sky


      if (IARGC() .ne. 15) then
      write(*,*) 'Error: Required parameters not specified properly!'
      stop 11
      endif

      hb=2.0
      br=1.0
      pib=pi/180.0

      do i=1,3
      call GETARG(i, fname)
      open(i,file=fname,status='old')
      enddo

      read(3,*)nrow,ncol
      call GETARG(4, fname)

      open(4,file=fname,status='old',
     &       access='direct',form='unformatted',recl=ncol)

      do i=1,8
      call GETARG(i+4, fname)
      open(i+10,file=fname,status='old',
     &          access='direct',recl=4*ncol)
      enddo
c
c    ouptput data
c
      call GETARG(13, fname)
! ULA -- TOA reflectance is not written to disk
!      open(52,file=fname, access='direct',recl=2*ncol)

      call GETARG(14, fname)
! ULA -- Lambertian reflectance is not written to disk
!      open(53,file=fname, access='direct',recl=2*ncol)

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

        dn = dn_1
        WHERE (dn_1 .lt. 0) dn = dn + 256

        do j=1,ncol

       if (dn(j) .gt.0.0 .and. a_mod(j) .gt. 0.0) then

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
     &        (1.0-fv(j))*(fs(j)*aa_solar+(1.0-fs(j))*aa_white))
     &        /aa_white

            a_eq=(1-aa_final)*s_mod(j)
            b_eq=aa_final
            c_eq=-ref_lm(j)

            if (abs(a_eq) .lt. 0.0000001) then
              ref_bar=-c_eq/b_eq
            else
              ref_bar=(-b_eq+sqrt(b_eq**2-4*a_eq*c_eq))/(2*a_eq)
            endif

            ref_brdfreal=ann*ref_bar/aa_white

            if (abs(ref_brdfreal-ref_top(j)) .lt. 0.5) then
              ref_brdf(j)=ref_bar*fnn/aa_white
              iref_brdf(j)=ref_brdf(j)*10000+0.5
            else
              ref_brdf(j)=-9.9
              iref_brdf(j)=-999
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

! ULA -- TOA and Lambertian reflectance are not written to disk
!        write(52,rec=i)(iref_top(j),j=1,ncol)
!        write(53,rec=i)(iref_lm(j),j=1,ncol)

        write(54,rec=i)(iref_brdf(j),j=1,ncol)
      enddo
      close(1)
      close(2)
      close(3)

! ULA -- TOA and Lambertian reflectance are not written to disk
!      close(52)
!      close(53)

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
     &          -pi/4.0
c
c     calculate Li-sparse
c
      tansolar=sinsolar/cossolar
      tanvia=sinvia/cosvia
      theta_new_v=atan(br*tanvia)
      theta_new_s=atan(br*tansolar)
      cosxi=cos(theta_new_s)*cos(theta_new_v)+sin(theta_new_s)
     &      * sin(theta_new_v)*cosra
       if (cosxi.ge.1) then
       cosxi=1
       endif
      secsolar=1.0/cos(theta_new_s)
      secvia=1.0/cos(theta_new_v)
      d_li2=abs(tan(theta_new_s)**2+tan(theta_new_v)**2-
     &          2.0*tan(theta_new_s)*tan(theta_new_v)*cosra)
      x_li=tan(theta_new_s)*tan(theta_new_v)*sinra
      cosl=hb*sqrt(d_li2+x_li**2)/(secsolar+secvia)
      if (cosl .ge. 1.0) then
        o_li=0.0
      else
         l_li=acos(cosl)
         o_li=(l_li-sin(l_li)*cos(l_li))*(secsolar+secvia)/pi
      endif
       li_sparse=o_li-(secsolar+secvia)+0.5*(1.0+cosxi)
     &           *secsolar*secvia
      RL_brdf=brdf0+brdf1*rs_thick+brdf2*li_sparse
      return
c
      end



      real function black_sky(brdf0,brdf1,brdf2,theta)
      real brdf0,brdf1,brdf2
      real theta
c     theta should be in radians
      black_sky=brdf0+brdf1*(-0.007574-0.070987*theta**2+
     &           0.307588*theta**3)+brdf2*(-1.284909-0.166314*
     &           theta**2+0.041840*theta**3)
      return
      end



      real function white_sky(brdf0,brdf1,brdf2)
      real brdf0,brdf1,brdf2
      white_sky=brdf0+brdf1*0.189184-brdf2*1.377622
      return
      end


