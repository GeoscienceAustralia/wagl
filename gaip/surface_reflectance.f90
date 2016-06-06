!program reflectance
SUBROUTINE reflectance( &
    nrow, ncol, &
    rori, &
    brdf0, brdf1, brdf2, &
    bias, slope_ca, &
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
    iref_terrain)

!   Calculates lambertian, brdf corrected and terrain corrected surface
!   reflectance.
  use reflectance_mod

  implicit none
!   input parameters
    integer nrow, ncol ! we should be passing in transposed arrays cols = rows and vice versa
    real*4 rori ! threshold for terrain correction
    real*4 brdf0, brdf1, brdf2 ! BRDF parameters
    real*4 bias, slope_ca ! satellite calibration coefficients
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

!f2py depend(nrow), ref_lm, ref_brdf, ref_terrain, dn
!f2py depend(nrow, ncol), dn_1, mask_self, mask_castsun, mask_castview
!f2py depend(nrow, ncol), solar_angle, sazi_angle, view_angle, rela_angle
!f2py depend(nrow, ncol), slope_angle,, aspect_angle, it_angle, et_angle
!f2py depend(nrow, ncol), rela_slope, a_mod, b_mod, s_mod, fv, fs, ts
!f2py depend(nrow, ncol), edir_h, edif_h
!f2py depend(nrow, ncol), iref_lm, iref_brdf, iref_terrain
!f2py intent(in) rori, brdf0, brdf1, brdf2, bias, slop_ca, ref_adj
!f2py intent(in) dn_1, mask_self, mask_castsun, mask_castview, solar_angle,
!f2py intent(in) sazi_angle, view_angle, rela_angle, slope_angle, aspect_angle
!f2py intent(in) it_angle, et_angle, rela_slope, a_mod, b_mod, s_mod, fv, fs, ts, edir_h, edif_h
!f2py intent(in) ref_lm, ref_brdf, ref_terrain, dn
!f2py intent(inout) iref_lm, iref_brdf, iref_terrain

!   internal parameters
    real, parameter :: hb = 2.0
    real, parameter :: br = 1.0
    integer, parameter :: veclen = 16

    !!! Tiled input quantities
    real*4 tile_dn_1(veclen) ! raw image (in real4 as replacement for dn array)
    integer*2 tile_mask_self(veclen) ! mask
    integer*2 tile_mask_castsun(veclen) ! self shadow mask
    integer*2 tile_mask_castview(veclen) ! cast shadow mask
!    real*4 tile_solar_angle(veclen) ! solar zenith angle
!    real*4 tile_sazi_angle(veclen) ! solar azimuth angle
!    real*4 tile_view_angle(veclen) ! view angle (for flat surface)
    real*4 tile_rela_angle(veclen) ! relative azimuth angle (for flat surface)
!    real*4 tile_slope_angle(veclen) ! slop angle
!    real*4 tile_aspect_angle(veclen) ! aspect angle
!    real*4 tile_it_angle(veclen) ! incident angle (for inclined surface)
!    real*4 tile_et_angle(veclen) ! exiting angle (for inclined surface)
    real*4 tile_rela_slope(veclen) ! relative angle (for inclined surface)
    real*4 tile_a_mod(veclen) ! modtran output (a)
    real*4 tile_b_mod(veclen) ! modtran output (b)
    real*4 tile_s_mod(veclen) ! modtran output (s)
    real*4 tile_fv(veclen) ! modtran output (fs)
    real*4 tile_fs(veclen) ! modtran output (fv)
    real*4 tile_ts(veclen) ! modtran output (ts)
    real*4 tile_edir_h(veclen) ! modtran output (direct irradiance)
    real*4 tile_edif_h(veclen) ! modtran output (diffuse irradiance)
    integer*2 tile_iref_lm(veclen) ! atmospheric corrected lambertial reflectance
    integer*2 tile_iref_brdf(veclen) ! atmospheric and brdf corrected reflectance
    integer*2 tile_iref_terrain(veclen) ! atmospheric and brdf and terrain correc

    real*4 tile_ref_lm(veclen)
    real*4 tile_ref_brdf(veclen)
    real*4 tile_ref_terrain(veclen)
!!! Tile sizing
    integer :: istart, iend, tilelen

    integer i, j
    real ann_f(veclen), aa_viewf(veclen), aa_solarf(veclen), aa_white
    real ann_s(veclen), aa_views(veclen), aa_solars(veclen)
    real lt(veclen), norm_1, norm_2
    real solar(veclen), sazi(veclen), view(veclen), slope(veclen), aspect(veclen), ra_lm(veclen), ra_sl(veclen), it(veclen), et(veclen)

    double precision a_eqf(veclen), b_eqf(veclen), c_eqf, ref_barf(veclen), aa_flat(veclen)
    double precision a_eqs(veclen), b_eqs(veclen), c_eqs, ref_bars(veclen), aa_slope(veclen)
    double precision ref_brdfrealf(veclen), ref_brdfreals(veclen)
    double precision bb_angle(veclen), cc_angle(veclen), tttt(veclen)
    real fnn, fnn_a(1), it_brdf(veclen), it_bk(veclen), et_brdf(veclen), et_bk(veclen)
    real angle_th(veclen)
    real vt(veclen), vd(veclen),angle(veclen), edir_t(veclen), eadj(veclen), edif_t(veclen), rdir(veclen), rdif(veclen), rtotal(veclen), rtotal_tmp(veclen), rdir_tmp(veclen)
    real cossolar(veclen),cosit(veclen),rth(veclen)
    
    logical :: tilemask(veclen)

!   li-sparse parameters


    norm_1 = brdf1 / brdf0
    norm_2 = brdf2 / brdf0

!   calculate white sky albedo
    aa_white = white_sky(1.0, norm_1, norm_2)
!   calcualte BRDF at 45 solar angle and 0 view angle
    fnn_a = RL_brdf( (/ 45 * pib /), (/ 0.0 /), (/ 0.0 /), hb, br, 1.0, norm_1, norm_2,1)
    fnn = fnn_a(1)

!   Now loop over the cols of the images
    do j=1,ncol
!----------------------------------------------------------------------
!   convert byte to integer for the raw image
!        do i=1,nrow
!            if (dn_1(i, j) .lt. 0) then
!                dn(i) = dn_1(i, j) + 65536
!            else
!                dn(i) = dn_1(i, j)
!            endif
!        enddo
       istart = 0
       iend = 0
       do while ( iend < nrow )
          istart = iend + 1
          iend = iend + veclen
          if ( iend > nrow ) iend = nrow
          tilelen = iend - istart + 1

          where ( dn_1(istart:iend,j) .lt. 0 ) 
             tile_dn_1(1:tilelen) = real(dn_1(istart:iend,j) + 65536,4)
          elsewhere
             tile_dn_1(1:tilelen) = real(dn_1(istart:iend,j),4)
          end where


          if( ANY( a_mod(istart:iend,j) .gt. 0 .AND. tile_dn_1(1:tilelen) .gt. 0.0 ) ) then

             tile_a_mod(1:tilelen) = a_mod(istart:iend,j)
             tile_b_mod(1:tilelen) = b_mod(istart:iend,j)
             tile_s_mod(1:tilelen) = s_mod(istart:iend,j)

             tile_fs(1:tilelen) = fs(istart:iend,j)
             tile_fv(1:tilelen) = fv(istart:iend,j)
             tile_ts(1:tilelen) = ts(istart:iend,j)

             tile_edir_h(1:tilelen) = edir_h(istart:iend,j)
             tile_edif_h(1:tilelen) = edif_h(istart:iend,j)
!!! Mask tiles
             tile_mask_self(1:tilelen) = mask_self(istart:iend,j)
             tile_mask_castsun(1:tilelen) = mask_castsun(istart:iend,j)
             tile_mask_castview(1:tilelen) = mask_castview(istart:iend,j)


!       now loop over the rows of the images
!        do i=1,nrow
!           if valid masks and valid digital number then do the calcs
!!! Mask after calculations
!            if (a_mod(i, j) .ge. 0 .and. dn(i) .gt. 0) then
             
!          if (rela_angle(i, j) .gt. 180) rela_angle(i, j) = rela_angle(i, j) - 360
!          if (rela_slope(i, j) .gt. 180) rela_slope(i, j) = rela_slope(i, j) - 360
             where ( rela_angle(istart:iend,j) .gt. 180 )
                tile_rela_angle(1:tilelen) = rela_angle(istart:iend,j) - 360
             elsewhere
                tile_rela_angle(1:tilelen) = rela_angle(istart:iend,j)
             end where
             
             where ( rela_slope(istart:iend,j) .gt. 180 )
                tile_rela_slope(1:tilelen) = rela_slope(istart:iend,j) - 360
             elsewhere
                tile_rela_slope(1:tilelen) = rela_slope(istart:iend,j)
             end where
             
!               convert angle to radians
!          solar = solar_angle(i, j) * pib
             solar(1:tilelen) = solar_angle(istart:iend,j)
             solar = solar * pib
          ! Trig is expensive
             cossolar = cos(solar)
!          sazi = sazi_angle(i, j) * pib
             sazi(1:tilelen) = sazi_angle(istart:iend,j)
             sazi = sazi * pib
 !         view = view_angle(i, j) * pib
             view(1:tilelen) = view_angle(istart:iend,j)
             view = view * pib
 !         slope = slope_angle(i, j) * pib
             slope(1:tilelen) = slope_angle(istart:iend,j)
             slope = slope * pib
 !         aspect = aspect_angle(i, j) * pib
             aspect(1:tilelen) = aspect_angle(istart:iend,j)
             aspect = aspect * pib
 !         ra_lm = rela_angle(i, j) * pib
             ra_lm(1:tilelen) = rela_angle(istart:iend,j)
             ra_lm = ra_lm * pib
 !         ra_sl = rela_slope(i, j) * pib
             ra_sl(1:tilelen) = rela_slope(istart:iend,j)
             ra_sl = ra_sl * pib
 !         it = it_angle(i, j) * pib
             it(1:tilelen) = it_angle(istart:iend,j)
 !         et = et_angle(i, j) * pib
             et(1:tilelen) = et_angle(istart:iend,j)
          
!          if (it_angle(i, j) .ge. 70.0) then
!             it_brdf= 70.0 * pib
!          else
!             it_brdf = it_angle(i, j) * pib
!          endif
             it_brdf = min(70.0,it) * pib
          
!          if (it_angle(i, j) .ge. 80.0) then
!             it_bk = 80.0 * pib
!          else
!             it_bk = it_angle(i, j) * pib
!          endif
             it_bk = min(80.0, it) * pib

             it = it * pib 
          ! Trig is expensive
             cosit = cos(it)

 !         if (et_angle(i, j) .ge. 60.0) then
 !            et_brdf = 60.0 * pib
 !         else
 !            et_brdf = et_angle(i, j) * pib
 !         endif
             et_brdf = min(60.0,et) * pib
                 
!          if (et_angle(i, j) .ge. 80.0) then
!             et_bk = 80.0 * pib
!          else
!             et_bk = et_angle(i, j) * pib
!          endif
             et_bk = min(80.0,et) * pib

             et = et * pib
!------------------------------------------------
!           for flat surface

!           calculate radiance at top atmosphere
!          lt = bias + dn(i) * slope_ca
             lt = bias + tile_dn_1 * slope_ca

!           calcualte lambetian reflectance with bilinear average

!          ref_lm(i) = (lt - b_mod(i, j)) / (a_mod(i, j) + s_mod(i, j) * &
!               (lt - b_mod(i, j)))
!          iref_lm(i, j) = ref_lm(i) * 10000 + 0.5
             tile_ref_lm = (lt - tile_b_mod) / ( tile_a_mod + tile_s_mod * &
                  ( lt - tile_b_mod ))
             tile_iref_lm = nint(tile_ref_lm * 10000.0,2)

!           set as small number if atmospheric corrected reflectance
!           below 0.001

!!!          Set these at the end
!            if (ref_lm(i).lt. 0.001) then
!                ref_lm(i) = 0.001
!                ref_brdf(i) = 0.001
!                iref_lm(i, j) = 10
!                iref_brdf(i, j) = 10
          !            else
!               calculate normalized BRDF shape function
             if( all(tile_ref_lm(1:tilelen) .lt.0.001) ) then
                tile_ref_lm = 0.001
                tile_ref_brdf = 0.001
                tile_iref_lm = 10
                tile_iref_brdf = 10
             else
                ann_f = RL_brdf(solar, view, ra_lm, hb, br, 1.0, norm_1, norm_2,veclen)
!               calculate black sky albedo for sloar angle
                aa_solarf = black_sky(1.0, norm_1, norm_2, solar,veclen)
!               calculate black sky albedo for view angle
                aa_viewf = black_sky(1.0, norm_1, norm_2, view,veclen)
!
                aa_flat = (tile_fv * (tile_fs * ann_f + (1.0 - tile_fs) * &
                     aa_viewf) + (1.0 - tile_fv) * (tile_fs * &
                     aa_solarf + (1.0 - tile_fs) * aa_white)) / aa_white

                a_eqf = (1 - aa_flat) * tile_s_mod * (1 - tile_s_mod * &
                     tile_ref_lm)
                b_eqf = aa_flat + tile_ref_lm * (1 - aa_flat) * tile_s_mod
!          c_eqf = -ref_lm(i)

!          if (abs(a_eqf) .lt. 0.0000001) then
!             ref_barf = tile_ref_lm / b_eqf
!          else
                ref_barf = (-b_eqf + &
                     sqrt(b_eqf**2 + 4 * a_eqf * tile_ref_lm)) / &
                     (2 * a_eqf)
                where( abs(a_eqf) .lt. 0.0000001 ) ref_barf = tile_ref_lm / b_eqf
!          endif
                ref_brdfrealf = ann_f * ref_barf / aa_white
                tile_ref_brdf = min(ref_barf * fnn / aa_white, 1.0)
                tile_iref_brdf = nint(tile_ref_brdf * 10000.0,2)
!          tile_iref_brdf = int(tile_ref_brdf * 10000.0 + 0.5,2)

!               this is to ensure that the brdf correction
          !               is the same as (or as close as possible to) the original NBAR version
!!! Condition from the top of this section.
                where( tile_ref_lm .lt.0.001 )
                   tile_ref_lm = 0.001
                   tile_ref_brdf = 0.001
                   tile_iref_lm = 10
                   tile_iref_brdf = 10
                end where
             end if
          
!          if (ref_brdf(i) .ge. 1) then
!             ref_brdf(i) = 1.0
!             iref_brdf(i, j) = ref_brdf(i) * 10000
!          endif

!       endif
!-------------------------------------------------------------------
!           conduct terrain correction
!!! Tile-wide mask at the start, per-element mask at the end
             tilemask(1:tilelen) = ((tile_mask_self(1:tilelen) .gt. 0) .and. &
                  (tile_mask_castsun(1:tilelen) .gt. 0)                .and. &
                  (tile_mask_castview(1:tilelen) .gt. 0)               .and. &
                  (it(1:tilelen) .lt. 85.0*pib )                       .and. &
                  (et(1:tilelen) .lt. 85.0*pib ) )
             if ( ANY( tilemask(1:tilelen) ) ) then
!----------------------------------------------------------
!               calculate vd and vt
                vd = 0.5 * (1.0 + cos(slope))
                vt = 1.0 - vd
!---------------------------------------------------------
!               calculate direct irradiance
!               Note the account taken of threshold
                edir_t = tile_edir_h * cosit / cossolar
!               calcualte adjacent irradiance for anisotropical surface
!               see Iqbal, 1983 "an introduction to solar
!               radiation"
                eadj = (tile_edir_h + tile_edif_h) * vt * ref_adj * &
                     (1.0 + sin(solar / 2.0)**2) * abs(cos(aspect - sazi))
! sin(0.5*theta)**2 = 0.5 + 0.5*cos(theta)
!---------------------------------------------------------------
!               sky diffuse irradiation for anisotropical surface
!               see Iqbal, 1983 "an introduction to solar
!               radiation" Hay model
                edif_t = tile_edif_h * (tile_ts * cosit / cossolar + &
                     vd * (1 - tile_ts)) + eadj
                rdir = edir_t / (tile_edir_h + tile_edif_h)
                rdif = edif_t / (tile_edir_h + tile_edif_h)
                rtotal = (edir_t + edif_t) / (tile_edir_h + tile_edif_h)
                
                rth = (rori - tile_s_mod * tile_ref_lm) / (1 - tile_s_mod * &
                     tile_ref_lm)
                

!                if (rtotal .le. rth) then
                bb_angle = tile_fs / cossolar + (1 - tile_fs) * &
                     tile_ts / cossolar
                cc_angle = -rth + (1.0 - tile_fs) * vd * &
                     (1.0 - tile_ts) + &
                     eadj / (tile_edir_h + tile_edif_h)
!                    tttt = -cc_angle / bb_angle
                tttt = max(min(-cc_angle / bb_angle,1.0),-1.0)
!                    if (tttt .gt. 1.0) tttt = 1.0
!                    if (tttt .lt. -1.0) tttt = -1.0
!                    angle_th = acos(tttt) * 180.0 / pi
                angle_th  = acos(tttt) 
!!! Do this in radians instead of degrees
!                    angle = 90.0 - it_angle(i, j) + angle_th
                angle = pi / 2.0 - it + angle_th
!                    edir_t = edir_h(i, j) * (cos(it) + cos(angle * pib)) / &
!                             (cos(solar) +cos(angle * pib))
                edir_t = tile_edir_h * ( cosit + cos(angle) ) /  &
                     ( cossolar + cos(angle) )
                rdir_tmp = edir_t/(tile_edir_h+tile_edif_h)
                rtotal_tmp = (edir_t+edif_t)/(tile_edir_h+tile_edif_h)
                where( rtotal .le. rth )
                   rtotal = rtotal_tmp
                   rdir = rdir_tmp
                end where
!                endif

!----------------------------------------------------------------
!               calculate normalized BRDF shape function for sloping surface
                ann_s = RL_brdf(it_brdf,et_brdf,ra_sl,hb,br,1.0,norm_1,norm_2,veclen)

!----------------------------------------------------------------
!               calculate black sky albedo for sloar angle
                aa_solars = black_sky(1.0,norm_1,norm_2,it_bk,veclen)
!--------------------------------------------------------------
!               calculate black sky albedo for view angle
                aa_views = black_sky(1.0,norm_1,norm_2,et_bk,veclen)
!-------------------------------------------------------------
                aa_slope = &
                     (rdir * (tile_fv * ann_s    + (1.0-tile_fv) * &
                     aa_solars) + rdif * (tile_fv * aa_views + &
                     (1.0-tile_fv) * aa_white)) / aa_white
                
                a_eqs = (rtotal - aa_slope) * tile_s_mod * &
                     (1 - tile_s_mod * tile_ref_lm)
                b_eqs = aa_slope + tile_ref_lm * (1 - aa_slope) * tile_s_mod
!                c_eqs = -ref_lm(i)

!                if (abs(a_eqs) .lt. 0.00001) then
!                    ref_bars = -c_eqs / b_eqs
!                else
                ref_bars = (-b_eqs + sqrt(b_eqs**2 +4 * a_eqs * tile_ref_lm)) / &
                     (2*a_eqs)
                where( abs(a_eqs) .lt. 0.00001 ) ref_bars = tile_ref_lm / b_eqs
!                endif
                ref_brdfreals = ann_s * ref_bars / aa_white
          
!          ref_terrain(i) = ref_bars * fnn / aa_white
!          iref_terrain(i, j) = int(ref_terrain(i) * 10000.0 + 0.5)
                tile_ref_terrain = min(ref_bars * fnn / aa_white, 1.0)
                tile_iref_terrain = nint(tile_ref_terrain * 10000.0,2)
!          if (ref_terrain(i) .ge. 1) then
!             ref_terrain(i) = 1.0
!             iref_terrain(i, j) = int(ref_terrain(i) * 10000.0 + 0.5)
!          endif

!               Should test for these cases in initial tests! (ie must be lt these)
!               presently comments as test for ge 90 in initial one
!
!          if (it_angle(i, j) .ge. 85.0) then
!             iref_terrain(i, j) = -999
!          endif
!          
!          if (et_angle(i, j) .ge. 85.0) then
!             iref_terrain(i, j) = -999
!          endif

!     if in masked or otherwise unused area you get here
!     put in -999 values to indicate null data
!       else
!          iref_terrain(i, j) = -999
!       endif
                where( it .ge. 85.0*pib          .OR.  &
                       et .ge. 85.0*pib          .OR.  &
                       tile_mask_self .le. 0     .OR.  &
                       tile_mask_castview .le. 0 .OR.  &
                       tile_mask_castsun .le. 0 )      &
                     tile_iref_terrain = -999
             else
                tile_iref_terrain = -999
             end if
         
             where (tile_a_mod(1:tilelen) .ge. 0 .and. tile_dn_1(1:tilelen) .gt. 0 )
                iref_lm(istart:iend,j) = tile_iref_lm(1:tilelen)
                iref_brdf(istart:iend,j) = tile_iref_brdf(1:tilelen)
                iref_terrain(istart:iend,j) = tile_iref_terrain(1:tilelen)
             elsewhere
!             ref_lm(i) = -999.0
!             ref_brdf(i) = -999.0
!             ref_terrain(i) = -999.0
                iref_lm(istart:iend, j) = -999
                iref_brdf(istart:iend, j) = -999
                iref_terrain(istart:iend, j) = -999
             end where
          else
             iref_lm(istart:iend, j) = -999
             iref_brdf(istart:iend, j) = -999
             iref_terrain(istart:iend, j) = -999
          end if
       enddo
    enddo
  END SUBROUTINE reflectance
