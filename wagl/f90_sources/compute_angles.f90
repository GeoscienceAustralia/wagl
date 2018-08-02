! subroutine cal_angles
SUBROUTINE cal_angles(lam_p,phip_p,tol_lam,orb_elements,spheroid, &
             smodel,track,num_tpts,timet,theta_p,azimuth,istat)

!     Cal_Angles is the main routine in the set
!     It finds the segment of the track model where a line from a point
!     hits the track line at an angle approx perpendicular to
!     the heading. It then refines the estimated crossing point
!     The angle from the track point to the given point
!     assumes a great circle relationship
!     Returned data are:
!          timet (time of track point -t0) secs
!          theta_p (look angle from the track to the given point) radians
!          azimuth (azimuth from track point to given point) radians
!          (The azimuth uses the BRDF convention)

!   * Re-written as an independent subroutine by JS, Aug 2014

!   Inputs:
!       lam_p
!       phip_p
!       tol_lam
!       orb_elements
!           1. Orbital inclination (degrees)
!           2. Semi_major radius (m)
!           3. Angular velocity (rad sec-1)
!       spheroid
!           1. Spheroid major axis
!           2. Inverse flattening
!           3. Eccentricity squared
!           4. Earth rotational angular velocity rad/sec
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
!   Outputs
!       timet
!       theta_p
!       azimuth
!       istat

    use sys_variables, only : pi, d2r, r2d

    implicit none

    double precision, intent(in) :: lam_p, phip_p, tol_lam
    double precision, dimension(4), intent(in) :: spheroid
    double precision, dimension(3), intent(in) :: orb_elements
!   smodel(phi0,phi0_p,rho0,t0,lam0,gamm0,beta0,rotn0,hxy0,N0,H0,th_ratio0)
    double precision, dimension(12), intent(in) :: smodel

    integer, intent(in) :: num_tpts
!   track(t,rho,phi_p,lam,beta,hxy,mj,skew)
    double precision, dimension(num_tpts,8), intent(in) :: track
    double precision, dimension(num_tpts) :: phi_p, lam, beta
    double precision, dimension(num_tpts) :: hxy, mj

!f2py depend(num_tpts), track, phi_p, lam, beta, hxy, mj

    double precision orad
    double precision asph, e2
    double precision t0, hxy0

    double precision, intent(out) :: timet, theta_p, azimuth
    integer, intent(out) :: istat

    integer west, neg, hit, j, negp
    integer left, right
    double precision rho_p, t_p, lamcal, betacal, adiff
    double precision delta_mu, tbet_p, phi_pt, lamt, dmu_prev, mjt, lam_test
    double precision phip_test, rho_c, t_c
    double precision theta_e, rn_p, ratio, temp

!   smodel parameters
    t0 = smodel(4)
    hxy0 = smodel(9)

!   track parameters
    phi_p(:) = track(:,3)
    lam(:) = track(:,4)
    beta(:) = track(:,5)
    hxy(:) = track(:,6)
    mj(:) = track(:,7)

!   Satellite orbital parameters
!   semi_major radius (m)
    orad = orb_elements(2)

!   Spheroid parameters
!   Spheroid major axis
!   Eccentricity squared
    asph = spheroid(1) !6378137d0
    e2 = spheroid(3)

!   Initialise the return values
    istat = 0
    timet = 0.0
    theta_p = 0.0
    azimuth = 0.0

    temp = phip_p
    call q_cal(temp,orb_elements,spheroid,smodel,rho_p,t_p,lamcal, &
           betacal,istat)

!   find where the azimuth of the track crossing at the same phi_p is
!   relative to input and use it to decide what side of the line
!   the point is or if it is on the track
    adiff = lam_p-lamcal

!   if on the track easy - done and go back
    if (abs(adiff) .lt. tol_lam) then
        timet = t_p-t0
        theta_p = 0.0d0
        azimuth = 0.0d0
        go to 99
    endif

!   look to see if the point is west or east of the track
    if (adiff .lt. 0.0) then
        west = 1
    else
        west = 0
    endif

!   now look over the sample points on the track
!   to find a bounding interval (the time interval between two track points)
!   for the solution
!
!   method is to define a line through given point expressed as
!   L(lam,phi)=0 and evaluate L at the track points
!   slope base on the above East/West point on the
!   the track
    delta_mu = 0.0d0

!   tbet_p is tan(beta)*hy/hx
    tbet_p = hxy0*tan(betacal)

    do j=1,num_tpts
        hit = 0
        phi_pt = dble(phi_p(j))
        lamt = dble(lam(j))
        if (j .gt. 1) then
            dmu_prev = delta_mu
        endif

!   Equation L(lamt,phi_pt)
        delta_mu = (phi_pt-phip_p)+tbet_p*(lamt-lam_p)
        if (j .eq. 1) then
            if (delta_mu .lt. 0.0) then
                neg = 1
            else
                neg = 0
            endif
            dmu_prev = delta_mu
        else
            if (delta_mu .lt. 0.0) then
                negp = 1
            else
                negp = 0
            endif
            if (negp .ne. neg) then
                hit = 1
                left = j-1
                right = j
                goto 10
            endif
        endif
    enddo

!   if no segment found it is bad news! Something is wrong.
    print*, 'error in cal_angles - no interval found in track'
    istat = 99
    go to 99

10  continue

!   Now solve for the crossing of a line from the given point
!   to the track interval (assumed a line joining ends)
!   using the local heading and scaling information
    tbet_p = (hxy(left)*tan(beta(left))+hxy(right)* &
      tan(beta(right)))/2.0d0
    mjt = mj(left)

    phi_pt = dble(phi_p(left))
    lamt = dble(lam(left))

    lam_test = ((phip_p-phi_pt)+tbet_p*lam_p+mjt*lamt)/(tbet_p+mjt)
    phip_test = phi_pt+mjt*(lam_test-lamt)

!   get info at the solution
    call q_cal(phip_test,orb_elements,spheroid,smodel,rho_c,t_c, &
           lamt,betacal,istat)

!   calculate the return results
!   theta_e is angle at Earth centre
!   theta_p is the look angle
    phi_pt = phip_test
    theta_e = acos(sin(phip_p)*sin(phi_pt)+ &
      cos(phip_p)*cos(phi_pt)*cos(lam_p-lamt))

    RN_p = dble(asph)/sqrt(1.0d0-dble(e2)*sin(phi_pt)**2)
    ratio = (RN_p*cos(theta_e))/(orad-RN_p*cos(theta_e))
    theta_p = atan(ratio*tan(theta_e))
    if (west.gt.0) then
        azimuth = betacal+(dble(pi)/2.0d0)
    else
        azimuth = 3.0d0*dble(pi)/2.0d0+betacal
    endif
    timet = t_c-t0

!   return to caller
99  continue
    return

END SUBROUTINE cal_angles
