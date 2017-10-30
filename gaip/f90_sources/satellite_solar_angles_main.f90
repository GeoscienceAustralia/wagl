! subroutine angle
subroutine angle(ncol,nrow,row_num,col_offset,alat,alon,spheroid,orb_elements, &
             hours,century,ntpoints,smodel,track, &
             view,azi,asol,soazi,rela_angle,tim,X_cent,N_cent,istat)

!   program to calculate solar, view and azimuth angle from
!   both UTM and lat/lon projection if we only know one point in the
!   satellite track. the program is written by Fuqin Li, Mar., 2014.
!   the subroutines cal_angle are from David Jupp

!   * Re-written as an F2Py subroutine by JS, Aug 2014

!   Inputs:
!       ncol
!       nrow
!       row_num
!       col_offset
!       alat
!       alon
!       spheroid
!           1. Spheroid major axis
!           2. Inverse flattening
!           3. Eccentricity squared
!           4. Earth rotational angular velocity rad/sec
!       orb_elements
!           1. Orbital inclination (degrees)
!           2. Semi_major radius (m)
!           3. Angular velocity (rad sec-1)
!       hours
!       century
!       ntpoints (number of time points created in determining the satellite track)
!       smodel
!           1. phi0
!           2. phi0_p
!           3. rho0
!           4. t0
!           5. lam0
!           6. gamm0
!           7. beta0
!           8. rotn0
!           9. hxy0
!           10. N0
!           11. H0
!           12. th_ratio0
!       track
!           1. t
!           2. rho
!           3. phi_p
!           4. lam
!           5. beta
!           6. hxy
!           7. mj
!           8. skew
!       view
!       azi
!       asol
!       soazi
!       rela_angle
!       tim
!       X_cent
!       N_cent
!
!   Outputs
!       view
!       azi
!       asol
!       soazi
!       rela_angle
!       tim
!       X_cent
!       N_cent
!       istat

    use sys_variables, only : pi, d2r, r2d

    implicit none

    integer, intent(in) :: ncol, nrow, row_num, col_offset
    double precision, dimension(ncol), intent(in) :: alat, alon
    double precision, dimension(4), intent(in) :: spheroid
    double precision, dimension(3), intent(in) :: orb_elements
    double precision, intent(in) :: hours, century
    integer, intent(in) :: ntpoints
    double precision, dimension(12), intent(in) :: smodel(12)
    double precision, dimension(12,8), intent(in) :: track(12,8)
    real, dimension(ncol), intent(inout) :: view, azi, asol, soazi, rela_angle, tim
    real, dimension(nrow), intent(inout) ::  X_cent, N_cent
    integer, intent(out) :: istat
    double precision delxx, tol_lam
    double precision xout, yout
    double precision phip_p, lam_p
    double precision timet, theta_p, azimuth
    integer j

!f2py depend(ncol), alat, alon, view, azi, asol, soazi, rela_angle, tim
!f2py depend(nrow), X_cent, N_cent

    do j=1,ncol
        xout = alon(j)
        yout = alat(j)

!       calculate pixel size (half)
        if (j .gt. 1) then
            delxx = (alon(j)-alon(j-1))/2
        else
            delxx = (alon(j+1)-alon(j))/2
        endif

        tol_lam = delxx*d2r*1.2

!       calculate solar angle
        call solar(yout*d2r, xout*d2r, century, hours, asol(j), &
               soazi(j))

!       go through the base sequence used in the test examples
        lam_p = xout*d2r
        call geod2geo(yout, orb_elements, spheroid, phip_p, istat)

        call cal_angles(lam_p, phip_p, tol_lam, orb_elements, &
               spheroid, smodel, track, ntpoints, timet, theta_p, &
               azimuth, istat)

!       We need a better mechanism to terminate upon error
!       and report the failure location
!       We'll try and use istat and let Python create an error
!       report if it found where any istat != 0
        if (istat .ne. 0) then
            print*, 'cal_angles failed at i,j=', row_num, j + col_offset
            goto 99
        endif

        tim(j) = timet
        view(j) = theta_p*r2d
        azi(j) = azimuth*r2d
        rela_angle(j) = azi(j)-soazi(j)
        if ((abs(timet) .gt. 1.0e-5) .and. (abs(view(j)) .lt. 1.0e-7) &
          .and. (abs(azi(j)) .lt. 1.0e-7)) then
            X_cent(row_num) = X_cent(row_num)+real(j + col_offset)
            N_cent(row_num) = N_cent(row_num)+1.0
        endif
    enddo

99  continue
    return

end subroutine angle
