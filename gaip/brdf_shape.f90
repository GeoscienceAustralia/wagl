real function RL_brdf (solar,view,ra,hb,br,brdf0,brdf1,brdf2)
!   calculate the normalised BRDF shape function
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
