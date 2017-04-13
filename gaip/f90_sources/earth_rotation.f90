SUBROUTINE cal_pole(theta, phi, theta_p, phi_p, thp, php)
! does some funky stuff to calculate various angles
    real theta, phi, theta_p, phi_p, thp, php, offset
    real pdiff, costhp, sinphp, cosphp, tnum, tden
    real eps, pi, d2r

    eps=1.0e-6
    pi=4.0*atan(1.0)
    d2r=pi/180.0

    if (abs(theta_p).le.eps) then
        thp=theta
        php=phi
        return
    endif

    offset = atan(tan(pi-d2r*phi_p)*cos(d2r*theta_p))
    pdiff = d2r*(phi-phi_p)
    costhp = cos(d2r*theta)*cos(d2r*theta_p)+sin(d2r*theta)*sin(d2r*theta_p)*cos(pdiff)

    if (costhp.ge.1.0-eps) then
        thp=0.0
        php=0.0-offset/d2r
        return
    else
        thp=acos(costhp)/d2r
        sinphp=sin(d2r*theta)*sin(pdiff)/sin(d2r*thp)
        cosphp=(cos(d2r*theta)*sin(d2r*theta_p)-sin(d2r*theta)*cos(d2r*theta_p)*cos(pdiff))/sin(d2r*thp)
        if (abs(sinphp).le.eps) then
            if (cosphp.gt.eps) then
                php=0.0-offset/d2r
            else
                php=180.0-offset/d2r
            endif
            return
        else if (abs(cosphp).le.eps) then
            if (sinphp.gt.eps) then
                php=90.0-offset/d2r
            else
                php=-90.0-offset/d2r
            endif
            return
        endif
    endif

    tnum=sin(d2r*theta)*sin(pdiff)
    tden=(cos(d2r*theta)*sin(d2r*theta_p)-sin(d2r*theta)*cos(d2r*theta_p)*cos(pdiff))
    php=(atan2(tnum,tden)-offset)/d2r
    return
end SUBROUTINE cal_pole
