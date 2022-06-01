!program reflectance
SUBROUTINE reflectance( &
    nrow, ncol, &
    rori, &
    norm_1, norm_2, &
    ref_adj, &
    no_data, &
    radiance, &
    shadow_mask, &
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
    ref_lm, &
    ref_brdf, &
    ref_terrain, &
    iref_lm, &
    iref_brdf, &
    iref_terrain, &
    norm_solar_zenith)

!   Calculates lambertian, brdf corrected and terrain corrected surface
!   reflectance.

!   input parameters
    integer nrow, ncol ! we should be passing in transposed arrays cols = rows and vice versa
    real*4 rori ! threshold for terrain correction
    real*4 ref_adj ! average reflectance for terrain correction
    real*4 no_data ! input & output no data value
    real*4 radiance(nrow, ncol) ! at sensor radiance image
    integer*1 shadow_mask(nrow, ncol) ! shadow mask
    real*4 solar_angle(nrow, ncol) ! solar zenith angle
    real*4 sazi_angle(nrow, ncol) ! solar azimuth angle
    real*4 view_angle(nrow, ncol) ! view angle (for flat surface)
    real*4 rela_angle(nrow, ncol) ! relative azimuth angle (for flat surface)
    real*4 slope_angle(nrow, ncol) ! slope angle
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
    real norm_solar_zenith ! solar zenith to normalize surface reflectance to

!internal parameters passed as arrays.
    real*4 ref_lm(nrow)
    real*4 ref_brdf(nrow)
    real*4 ref_terrain(nrow)

!f2py depend(nrow), ref_lm, ref_brdf, ref_terrain
!f2py depend(nrow, ncol), radiance, shadow_mask
!f2py depend(nrow, ncol), solar_angle, sazi_angle, view_angle, rela_angle
!f2py depend(nrow, ncol), slope_angle,, aspect_angle, it_angle, et_angle
!f2py depend(nrow, ncol), rela_slope, a_mod, b_mod, s_mod, fv, fs, ts
!f2py depend(nrow, ncol), edir_h, edif_h
!f2py depend(nrow, ncol), iref_lm, iref_brdf, iref_terrain
!f2py intent(in) rori, norm_1, norm_2, ref_adj, no_data
!f2py intent(in) dn_1, shadow_mask, solar_angle,
!f2py intent(in) sazi_angle, view_angle, rela_angle, slope_angle, aspect_angle
!f2py intent(in) it_angle, et_angle, rela_slope, a_mod, b_mod, s_mod, fv, fs, ts, edir_h, edif_h
!f2py intent(in) ref_lm, ref_brdf, ref_terrain, dn
!f2py intent(inout) iref_lm, iref_brdf, iref_terrain

!   internal parameters
    integer i, j, i_no_data
    real pi, pib
    real ann_f, aa_viewf, aa_solarf, aa_white
    real ann_s, aa_views, aa_solars
    real lt, norm_1, norm_2
    real solar, sazi, view, slope, aspect, ra_lm, ra_sl, it, et

    double precision a_eqf, b_eqf, c_eqf, ref_barf, aa_flat
    double precision a_eqs, b_eqs, c_eqs, ref_bars, aa_slope
    double precision ref_brdfrealf, ref_brdfreals
    double precision bb_angle, cc_angle, tttt
    real hb, br, fnn, it_brdf, it_bk, et_brdf, et_bk
    real angle_th
    real vt, vd,angle, edir_t, eadj, edif_t, rdir, rdif, rtotal
    real RL_brdf, black_sky, white_sky
    real cosslope, rth
    external RL_brdf, black_sky, white_sky

!   li-sparse parameters
    hb = 2.0
    br = 1.0
    pi = 4 * atan(1.0)
    pib = pi / 180.0

!   integer version of the no_data value
    i_no_data = int(no_data)

!   calculate white sky albedo
    aa_white = white_sky(1.0, norm_1, norm_2)
!   calcualte BRDF at 45 solar angle and 0 view angle
    fnn = RL_brdf(norm_solar_zenith * pib, 0.0, 0.0, hb, br, 1.0, norm_1, norm_2)
!    print*,fnn

!   Now loop over the cols of the images
    do j=1,ncol
!----------------------------------------------------------------------

!       now loop over the rows of the images
        do i=1,nrow
!           radiance value for the pixel
            lt = radiance(i, j)

!           TODO: check radiance against a better null than 0
!           TODO: some at sensor radiance values are < 0
!           if valid masks and valid digital number then do the calcs
            if (a_mod(i, j) .ge. 0 .and. lt .ne. no_data) then
                if (rela_angle(i, j) .gt. 180) rela_angle(i, j) = rela_angle(i, j) - 360
                if (rela_slope(i, j) .gt. 180) rela_slope(i, j) = rela_slope(i, j) - 360

!               convert angle to radians
                solar = solar_angle(i, j) * pib
                sazi = sazi_angle(i, j) * pib
                view = view_angle(i, j) * pib
                slope = slope_angle(i, j) * pib
                aspect = aspect_angle(i, j) * pib
                ra_lm = rela_angle(i, j) * pib
                ra_sl = rela_slope(i, j) * pib
                it = it_angle(i, j) * pib
                et = et_angle(i, j) * pib

                if (it_angle(i, j) .ge. 70.0) then
                    it_brdf = 70.0 * pib
                else
                    it_brdf = it_angle(i, j) * pib
                endif

                if (it_angle(i, j) .ge. 80.0) then
                    it_bk = 80.0 * pib
                else
                    it_bk = it_angle(i, j) * pib
                endif

                if (et_angle(i, j) .ge. 60.0) then
                    et_brdf = 60.0 * pib
                else
                    et_brdf = et_angle(i, j) * pib
                endif

                if (et_angle(i, j) .ge. 80.0) then
                    et_bk = 80.0 * pib
                else
                et_bk = et_angle(i, j) * pib
            endif
!------------------------------------------------
!           for flat surface

!           calcualte lambetian reflectance with bilinear average

            ref_lm(i) = (lt - b_mod(i, j)) / (a_mod(i, j) + s_mod(i, j) * &
                        (lt - b_mod(i, j)))
            iref_lm(i, j) = ref_lm(i) * 10000 + 0.5

!           set as small number if atmospheric corrected reflectance
!           below 0.0001

            if (ref_lm(i).lt. 0.0001) then
                ref_lm(i) = 0.0001
                ref_brdf(i) = 0.0001
                iref_lm(i, j) = 1
                iref_brdf(i, j) = 1
            else
!               calculate normalized BRDF shape function
                ann_f = RL_brdf(solar, view, ra_lm, hb, br, 1.0, norm_1, norm_2)
!               calculate black sky albedo for sloar angle
                aa_solarf = black_sky(1.0, norm_1, norm_2, solar)
!               calculate black sky albedo for view angle
                aa_viewf = black_sky(1.0, norm_1, norm_2, view)
!
                aa_flat = (fv(i, j) * (fs(i, j) * ann_f + (1.0 - fs(i, j)) * &
                          aa_viewf) + (1.0 - fv(i, j)) * (fs(i, j) * &
                          aa_solarf + (1.0 - fs(i, j)) * aa_white)) / aa_white

                a_eqf = (1 - aa_flat) * s_mod(i, j) * (1 - s_mod(i, j) * &
                        ref_lm(i))
                b_eqf = aa_flat + ref_lm(i) * (1 - aa_flat) * s_mod(i, j)
                c_eqf = -ref_lm(i)

                if (abs(a_eqf) .lt. 0.0000001) then
                    ref_barf = -c_eqf / b_eqf
                else
                    ref_barf = (-b_eqf + &
                                sqrt(b_eqf**2 - 4 * a_eqf * c_eqf)) / &
                               (2 * a_eqf)
                endif
                    ref_brdfrealf = ann_f * ref_barf / aa_white
                    ref_brdf(i) = ref_barf * fnn / aa_white
                    iref_brdf(i, j) = ref_brdf(i) * 10000 + 0.5

!               this is to ensure that the brdf correction
!               is the same as (or as close as possible to) the original NBAR version
                if (ref_brdf(i) .ge. 1) then
                  ref_brdf(i) = 1.0
                  iref_brdf(i, j) = ref_brdf(i) * 10000
                endif

            endif
!-------------------------------------------------------------------
!           conduct terrain correction
            if ((shadow_mask(i, j) .gt. 0) .and. &
                (it_angle(i, j) .lt. 90.0) &
                .and. (et_angle(i, j) .lt. 90.0)) then
!----------------------------------------------------------
                cosslope = cos(slope)
!               calculate vd and vt
                vd = 0.5 * (1.0 + cosslope)
                vt = 1.0 - vd
!---------------------------------------------------------
!               calculate direct irradiance
!               Note the account taken of threshold

                edir_t = edir_h(i, j) * cos(it) / cos(solar)
!               calculate adjacent irradiance for anisotropical surface
!               see Iqbal, 1983 "an introduction to solar
!               radiation"
                eadj = (edir_h(i, j) + edif_h(i, j)) * vt * ref_adj * &
                       (1.0 + sin(solar / 2.0)**2) * abs(cos(aspect - sazi))
!---------------------------------------------------------------
!               sky diffuse irradiation for anisotropical surface
!               see Iqbal, 1983 "an introduction to solar
!               radiation" Hay model
                edif_t = edif_h(i, j) * (ts(i, j) * cos(it) / cos(solar) + &
                         vd * (1 - ts(i, j))) + eadj
                rdir = edir_t / (edir_h(i, j) + edif_h(i, j))
                rdif = edif_t / (edir_h(i, j) + edif_h(i, j))
                rtotal = (edir_t + edif_t) / (edir_h(i, j) + edif_h(i, j))

                rth = (rori - s_mod(i, j) * ref_lm(i)) / (1 - s_mod(i, j) * &
                      ref_lm(i))

                if (rtotal .le. rth) then
                    bb_angle = fs(i, j) / cos(solar) + (1 - fs(i, j)) * &
                               ts(i, j) / cos(solar)
                    cc_angle = -rth + (1.0 - fs(i, j)) * vd * &
                               (1.0 - ts(i, j)) + &
                               eadj / (edir_h(i, j) + edif_h(i, j))
                    tttt = -cc_angle / bb_angle
                    if (tttt .gt. 1.0) tttt = 1.0
                    if (tttt .lt. -1.0) tttt = -1.0
                    angle_th = acos(tttt) * 180.0 / pi
                    angle = 90.0 - it_angle(i, j) + angle_th
                    edir_t = edir_h(i, j) * (cos(it) + cos(angle * pib)) / &
                             (cos(solar) +cos(angle * pib))
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
                    (rdir * (fv(i, j) * ann_s    + (1.0-fv(i, j)) * &
                     aa_solars) + rdif * (fv(i, j) * aa_views + &
                     (1.0-fv(i, j)) * aa_white)) / aa_white

                a_eqs = (rtotal - aa_slope) * s_mod(i, j) * &
                        (1 - s_mod(i, j) * ref_lm(i))
                b_eqs = aa_slope + ref_lm(i) * (1 - aa_slope) * s_mod(i, j)
                c_eqs = -ref_lm(i)

                if (abs(a_eqs) .lt. 0.00001) then
                    ref_bars = -c_eqs / b_eqs
                else
                    ref_bars = (-b_eqs + sqrt(b_eqs**2 -4 * a_eqs * c_eqs)) / &
                               (2*a_eqs)
                endif
                ref_brdfreals = ann_s * ref_bars / aa_white
                ref_terrain(i) = ref_bars * fnn / aa_white
                iref_terrain(i, j) = int(ref_terrain(i) * 10000.0 + 0.5)

                if (ref_terrain(i) .ge. 1) then
                    ref_terrain(i) = 1.0
                    iref_terrain(i, j) = int(ref_terrain(i) * 10000.0 + 0.5)
                endif

!               set terrain corrected reflectance less than 0.0001 to 0.0001
                if (ref_terrain(i) .lt. 0.0001) then 
                    ref_terrain(i) = 0.0001
                    iref_terrain(i, j) = 1 
                endif 

!               Should test for these cases in initial tests! (ie must be lt these)
!               presently comments as test for ge 90 in initial one
!
                if (it_angle(i, j) .ge. 85.0) then
                    iref_terrain(i, j) = i_no_data
                endif

                if (et_angle(i, j) .ge. 85.0) then
                    iref_terrain(i, j) = i_no_data
                endif

!     if in masked or otherwise unused area you get here
!     put in no_data values to indicate null data
                else
                    iref_terrain(i, j) = i_no_data
                endif
            else
                ref_lm(i) = no_data
                ref_brdf(i) = no_data
                ref_terrain(i) = no_data
                iref_lm(i, j) = i_no_data
                iref_brdf(i, j) = i_no_data
                iref_terrain(i, j) = i_no_data
            endif
        enddo
    enddo
END SUBROUTINE reflectance
