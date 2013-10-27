!program terrain_correction
SUBROUTINE terrain_correction( &
    nrow, ncol, &
    rori, &
    brdf0, brdf1, brdf2, &
    bias, slope_ca, esun, dd, &
    ref_adj, &
    dn_1, &
    mask_self, &
    mask_castsun, &
    mask_castview, &
    solar_angle, &
    sazi_angle, &
    view_angle, &
    rela_angle, &
    slope_angle, &
    aspect_angle, &
    it_angle, &
    et_angle, &
    rela_slope, &
    a_mod, &
    b_mod, &
    s_mod, &
    fs, &
    fv, &
    ts, &
    edir_h, &
    edif_h, &
    dn, &
    ref_lm, &
    ref_brdf, &
    ref_terrain, &
    iref_lm, &
    iref_brdf, &
    iref_terrain &
)
!   input parameters
    integer nrow, ncol
    real*4 rori ! threshold for terrain correction
    real*4 brdf0, brdf1, brdf2 ! BRDF parameters
    real*4 bias, slope_ca, esun, dd ! satellite calibration coefficients
    real*4 ref_adj ! average reflectance for terrain correction
    integer*2 dn_1(nrow, ncol) ! raw image
    integer*2 mask_self(nrow, ncol) ! mask
    integer*2 mask_castsun(nrow, ncol) ! self shadow mask
    integer*2 mask_castview(nrow, ncol) ! cast shadow mask
    real*4 solar_angle(nrow, ncol) ! solar zenith angle
    real*4 sazi_angle(nrow, ncol) ! solar azimuth angle
    real*4 view_angle(nrow, ncol) ! view angle (for flat surface)
    real*4 rela_angle(nrow, ncol) ! relative azimuth angle (for flat surface)
    real*4 slope_angle(nrow, ncol) ! slop angle
    real*4 aspect_angle(nrow, ncol) ! aspect angle
    real*4 it_angle(nrow, ncol) ! incident angle (for inclined surface)
    real*4 et_angle(nrow, ncol) ! exiting angle (for inclined surface)
    real*4 rela_slope(nrow, ncol) ! relative angle (for inclined surface)
    real*4 a_mod(nrow, ncol) ! modtran output (a)
    real*4 b_mod(nrow, ncol) ! modtran output (b)
    real*4 s_mod(nrow, ncol) ! modtran output (s)
    real*4 fv(nrow, ncol) ! modtran output (fs)
    real*4 fs(nrow, ncol) ! modtran output (fv)
    real*4 ts(nrow, ncol) ! modtran output (ts)
    real*4 edir_h(nrow, ncol) ! modtran output (direct irradiance)
    real*4 edif_h(nrow, ncol) ! modtran output (diffuse irradiance)
    integer*2 iref_lm(nrow, ncol) ! atmospheric corrected lambertial reflectance
    integer*2 iref_brdf(nrow, ncol) ! atmospheric and brdf corrected reflectance
    integer*2 iref_terrain(nrow, ncol) ! atmospheric and brdf and terrain corrected reflectance

!internal parameters passed as arrays.
    integer*4 dn(nrow)
    real*4 ref_lm(nrow)
    real*4 ref_brdf(nrow)
    real*4 ref_terrain(nrow)
!f2py integer intent(hide),depend(dn_1) :: nrow=shape(dn_1,0), ncol=shape(dn_1,1)
!f2py intent(in) rori, brdf0, brdf1, brdf2, bias, slop_ca, esun, dd, ref_adj
!f2py intent(in) dn_1, mask_self, mask_castsun, mask_castview, solar_angle,
!f2py intent(in) sazi_angle, view_angle, rela_angle, slope_angle, aspect_angle
!f2py intent(in) it_angle, et_angle, rela_slope, a_mod, b_mod, s_mod, fv, fs, ts, edir_h, edif_h
!f2py intent(hide) ref_lm, ref_brdf, ref_terrain, dn
!f2py intent(out) iref_lm, iref_brdf, iref_terrain

! equivalent to the following must be done in Python.
!----------------------------------------------------------------------
!     open brdf parameter file, startend file and coordinator file
!----------------------------------------------------------------------
!    do i=1,3
!        call GETARG(i, fname)
!        open(i,file=fname,status='old')
!    enddo
!    read(3,*)nrow,ncol
!----------------------------------------------------------------------
!   read threshold for terrain correction
!    read(1,*)rori
!   read BRDF parameters
!    read(1,*)brdf0,brdf1,brdf2
!   read satellite calibration coefficients
!    read(1,*)bias,slope_ca,esun,dd
!   read average reflectance for terrain correction
!    read(1,*)ref_adj

!   internal parameters
    integer i, j
    real pi, pib
    real ann_f, aa_viewf, aa_solarf, aa_white
    real ann_s, aa_views, aa_solars
    real lt, norm_1, norm_2
    real solar, sazi, view, slope, aspect, ra_lm, ra_sl, it, et

    double precision a_eqf,b_eqf,c_eqf,ref_barf,aa_flat
    double precision a_eqs,b_eqs,c_eqs,ref_bars,aa_slope
    double precision ref_brdfrealf,ref_brdfreals
    double precision aa_angle,bb_angle,cc_angle,tttt
    real hb,br,fnn,it_brdf,it_bk,et_brdf,et_bk
    real angle_th
    real vt,vd,angle,edir_t,eadj,edif_t,rdir,rdif,rtotal
    real RL_brdf,black_sky,white_sky
    real cosslope,ff,rth
    external RL_brdf,black_sky,white_sky

!   li-sparse parameters
    hb=2.0
    br=1.0
    pi=4*atan(1.0)
    pib=pi/180.0

    norm_1=brdf1/brdf0
    norm_2=brdf2/brdf0
    print*,rori

!   calculate white sky albedo
    aa_white=white_sky(1.0,norm_1,norm_2)
!   calcualte BRDF at 45 solar angle and 0 view angle
    fnn=RL_brdf(45*pib,0.0,0.0,hb,br,1.0,norm_1,norm_2)
!    print*,fnn

!   Now loop over the colss of the images
    do j=1,ncol
!----------------------------------------------------------------------
!   convert byte to integer for the raw image
        do i=1,nrow
            if (dn_1(i, j) .lt. 0) then
                dn(i) = dn_1(i, j)+65536
            else
                dn(i) = dn_1(i, j)
            endif
        enddo

!       now loop over the rows of the images
        do i=1,nrow
!           if valid masks and valid digital number then do the calcs
            if (a_mod(i, j) .ge. 0 .and. dn(i) .gt. 0) then
                if (rela_angle(i, j) .gt. 180) rela_angle(i, j) = rela_angle(i, j)-360
                if (rela_slope(i, j) .gt. 180) rela_slope(i, j) = rela_slope(i, j)-360

!               convert angle to radians
                solar = solar_angle(i, j)*pib
                sazi = sazi_angle(i, j)*pib
                view = view_angle(i, j)*pib
                slope = slope_angle(i, j)*pib
                aspect = aspect_angle(i, j)*pib
                ra_lm = rela_angle(i, j)*pib
                ra_sl = rela_slope(i, j)*pib
                it = it_angle(i, j)*pib
                et = et_angle(i, j)*pib

                if (it_angle(i, j) .ge.70.0) then
                    it_brdf=70.0*pib
                else
                    it_brdf = it_angle(i, j)*pib
                endif

                if (it_angle(i, j) .ge.80.0) then
                    it_bk = 80.0*pib
                else
                    it_bk = it_angle(i, j)*pib
                endif

                if (et_angle(i, j) .ge.60.0) then
                    et_brdf = 60.0*pib
                else
                    et_brdf = et_angle(i, j)*pib
                endif

                if (et_angle(i, j) .ge.80.0) then
                    et_bk = 80.0*pib
                else
                et_bk = et_angle(i, j)*pib
            endif
!------------------------------------------------
!           for flat surface

!           calculate radiance at top atmosphere
            lt = bias+dn(i)*slope_ca

!           calcualte lambetian reflectance with bilinear average

            ref_lm(i) = (lt-b_mod(i, j))/(a_mod(i, j)+s_mod(i, j)*(lt-b_mod(i, j)))
            iref_lm(i, j) = ref_lm(i)*10000+0.5

!           set as small number if atmospheric corrected reflectance
!           below 0.001

            if (ref_lm(i).lt. 0.001) then
                ref_lm(i) = 0.001
                ref_brdf(i) = 0.001
                iref_lm(i, j) = 10
                iref_brdf(i, j) = 10
            else
!               calculate normalized BRDF shape function
                ann_f = RL_brdf(solar,view,ra_lm,hb,br,1.0,norm_1,norm_2)
!               calculate black sky albedo for sloar angle
                aa_solarf = black_sky(1.0,norm_1,norm_2,solar)
!               calculate black sky albedo for view angle
                aa_viewf = black_sky(1.0,norm_1,norm_2,view)
!
                aa_flat = &
                    (fv(i, j)     *(fs(i, j)*ann_f     + (1.0-fs(i, j))*aa_viewf) + &
                    (1.0-fv(i, j))*(fs(i, j)*aa_solarf + (1.0-fs(i, j))*aa_white)) / aa_white

                a_eqf = (1-aa_flat)*s_mod(i, j)*(1-s_mod(i, j)*ref_lm(i))
                b_eqf = aa_flat+ref_lm(i)*(1-aa_final)*s_mod(i, j)
                c_eqf = -ref_lm(i)

                if (abs(a_eqf) .lt. 0.0000001) then
                    ref_barf=-c_eqf/b_eqf
                else
                    ref_barf = (-b_eqf+sqrt(b_eqf**2-4*a_eqf*c_eqf))/(2*a_eqf)
                endif
                    ref_brdfrealf = ann_f*ref_barf/aa_white
                    ref_brdf(i) = ref_barf*fnn/aa_white
                    iref_brdf(i, j) = ref_brdf(i)*10000+0.5

!               this is to ensure that the brdf correction
!               is the same as (or as close as possible to) the original NBAR version
                if (ref_brdf(i) .ge. 1) then
                  ref_brdf(i)=1.0
                  iref_brdf(i, j)=ref_brdf(i)*10000
                endif

            endif
!-------------------------------------------------------------------
!           conduct terrain correction
            if ((mask_self(i, j) .gt. 0) .and. &
                (mask_castsun(i, j).gt. 0) .and. &
                (mask_castview(i, j) .gt. 0) .and. (it_angle(i, j) .lt. 90.0) .and. &
                (et_angle(i, j) .lt. 90.0)) then
!----------------------------------------------------------
                cosslope = cos(slope)
!               calculate vd and vt
                vd = 0.5*(1.0+cosslope)
                vt = 1.0-vd
!---------------------------------------------------------
!               calculate direct irradiance
!               Note the account taken of threshold

                edir_t = edir_h(i, j)*cos(it)/cos(solar)
!               calcualte adjacent irradiance for anisotropical surface
!               see Iqbal, 1983 "an introduction to solar
!               radiation"
                eadj = (edir_h(i, j)+edif_h(i, j))*vt*ref_adj*(1.0+sin(solar/2.0)**2)*abs(cos(aspect-sazi))
!---------------------------------------------------------------
!               sky diffuse irradiation for anisotropical surface
!               see Iqbal, 1983 "an introduction to solar
!               radiation" Hay model
                edif_t = edif_h(i, j)*(ts(i, j)*cos(it)/cos(solar)+vd*(1-ts(i, j)))+eadj
                rdir = edir_t/(edir_h(i, j)+edif_h(i, j))
                rdif = edif_t/(edir_h(i, j)+edif_h(i, j))
                rtotal = (edir_t+edif_t)/(edir_h(i, j)+edif_h(i, j))

                rth = (rori-s_mod(i, j)*ref_lm(i)) / (1-s_mod(i, j)*ref_lm(i))

                if (rtotal .le. rth) then
                    bb_angle = fs(i, j)/cos(solar)+(1-fs(i, j))*ts(i, j)/cos(solar)
                    cc_angle = -rth+(1.0-fs(i, j))*vd*(1.0-ts(i, j)) + eadj/(edir_h(i, j)+edif_h(i, j))
                    tttt = -cc_angle/bb_angle
                    if (tttt .gt. 1.0) tttt = 1.0
                    if (tttt .lt. -1.0) tttt = -1.0
                    angle_th = acos(tttt)*180.0/pi
                    angle = 90.0 - it_angle(i, j) + angle_th
                    edir_t = edir_h(i, j)*(cos(it)+cos(angle*pib)) / (cos(solar)+cos(angle*pib))
                    rdir = edir_t/(edir_h(i, j)+edif_h(i, j))
                    rtotal = (edir_t+edif_t)/(edir_h(i, j)+edif_h(i, j))
                endif

!----------------------------------------------------------------
!               calculate normalized BRDF shape function for sloping surface
                ann_s = RL_brdf(it_brdf,et_brdf,ra_sl,hb,br,1.0,norm_1,norm_2)
!----------------------------------------------------------------
!               calculate black sky albedo for sloar angle
                aa_solars = black_sky(1.0,norm_1,norm_2,it_bk)
!--------------------------------------------------------------
!               calculate black sky albedo for view angle
                aa_views = black_sky(1.0,norm_1,norm_2,et_bk)
!-------------------------------------------------------------
                aa_slope = &
                    (rdir*(fv(i, j)*ann_s    + (1.0-fv(i, j))*aa_solars) + &
                     rdif*(fv(i, j)*aa_views + (1.0-fv(i, j))*aa_white)) / aa_white

                a_eqs = (rtotal-aa_slope)*s_mod(i, j)*(1-s_mod(i, j)*ref_lm(i))
                b_eqs = aa_slope+ref_lm(i)*(1-aa_slope)*s_mod(i, j)
                c_eqs = -ref_lm(i)

                if (abs(a_eqs) .lt. 0.00001) then
                    ref_bars = -c_eqs/b_eqs
                else
                    ref_bars = (-b_eqs+sqrt(b_eqs**2-4*a_eqs*c_eqs))/(2*a_eqs)
                endif
                ref_brdfreals = ann_s*ref_bars/aa_white
                ref_terrain(i) = ref_bars*fnn/aa_white
                iref_terrain(i, j) = int(ref_terrain(i)*10000.0+0.5)

                if (ref_terrain(i) .ge. 1) then
                    ref_terrain(i) = 1.0
                    iref_terrain(i, j) = int(ref_terrain(i)*10000.0+0.5)
                endif

!               Should test for these cases in initial tests! (ie must be lt these)
!               presently comments as test for ge 90 in initial one
!
                if (it_angle(i, j) .ge. 85.0) then
                    iref_terrain(i, j)=-999
                endif

                if (et_angle(i, j) .ge. 85.0) then
                    iref_terrain(i, j)=-999
                endif

!     if in masked or otherwise unused area you get here
!     put in -999 values to indicate null data
                else
                    iref_terrain(i, j)=-999
                endif
            else
                ref_lm(i)=-999.0
                ref_brdf(i)=-999.0
                ref_terrain(i)=-999.0
                iref_lm(i, j)=-999
                iref_brdf(i, j)=-999
                iref_terrain(i, j)=-999
            endif
        enddo
        !write(55,rec=i)(ref_lm(i, j),j=1,ncol)
        !write(56,rec=i)(iref_brdf(i, j),j=1,ncol)
        !write(57,rec=i)(iref_terrain(i, j),j=1,ncol)
    enddo
END SUBROUTINE terrain_correction





real function RL_brdf (solar,view,ra,hb,br,brdf0,brdf1,brdf2)
    real pi
    real solar, view, ra, hb, br
    real brdf0,brdf1,brdf2
    real rs_thick,li_sparse,secsolar,secvia
    real cossolar,cosvia,cosra,sinsolar,sinvia,sinra
    real cosxi,xi,tansolar,tanvia,theta_new_v,theta_new_s
    real d_li2,x_li,cosl,l_li,o_li

    RL_brdf=0.0

    pi = 4 * atan(1.0)
!   calculate Ross-thick kernel
    cossolar = cos(solar)
    cosvia = cos(view)
    cosra = cos(ra)
    sinsolar = sin(solar)
    sinvia = sin(view)
    sinra = sin(ra)
    cosxi = cossolar*cosvia+sinsolar*sinvia*cosra
    if (cosxi .ge.1) then
        cosxi = 1
    endif
    xi = acos(cosxi)
    rs_thick = ((pi/2.0-xi)*cos(xi)+sin(xi)) / (cossolar+cosvia) - pi/4.0

!   alculate Li-sparse
    tansolar = sinsolar/cossolar
    tanvia = sinvia/cosvia
    theta_new_v = atan(br*tanvia)
    theta_new_s = atan(br*tansolar)
    cosxi = cos(theta_new_s)*cos(theta_new_v)+sin(theta_new_s)*sin(theta_new_v)*cosra
    if (cosxi.ge.1) then
        cosxi = 1
    endif
    secsolar = 1.0/cos(theta_new_s)
    secvia = 1.0/cos(theta_new_v)
    d_li2 = abs(tan(theta_new_s)**2+tan(theta_new_v)**2 - 2.0*tan(theta_new_s)*tan(theta_new_v)*cosra)
    x_li = tan(theta_new_s)*tan(theta_new_v)*sinra
    cosl = hb*sqrt(d_li2+x_li**2)/(secsolar+secvia)
    if (cosl .ge. 1.0) then
        o_li=0.0
    else
        l_li=acos(cosl)
        o_li=(l_li-sin(l_li)*cos(l_li))*(secsolar+secvia)/pi
    endif
       li_sparse=o_li-(secsolar+secvia)+0.5*(1.0+cosxi) * secsolar*secvia
        RL_brdf=brdf0+brdf1*rs_thick+brdf2*li_sparse
    return
end





real function black_sky(brdf0,brdf1,brdf2,theta)
      real brdf0,brdf1,brdf2
      real theta

!   theta should be in radians
    black_sky = brdf0 + &
        brdf1*(-0.007574-0.070987*theta**2 + 0.307588*theta**3) + &
        brdf2*(-1.284909-0.166314*theta**2 + 0.041840*theta**3)
    return
end





real function white_sky(brdf0, brdf1, brdf2)
    real brdf0, brdf1, brdf2
    white_sky = brdf0 + brdf1*0.189184 - brdf2*1.377622
    return
end

