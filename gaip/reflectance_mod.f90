module reflectance_mod

  real, parameter :: pi = 3.14159265358979
  real, parameter :: pib = 0.017453292519943 ! pi / 180

contains
    
  function RL_brdf (solar,view,ra,hb,br,brdf0,brdf1,brdf2, veclen)
!   calculate the normalised BRDF shape function
    integer :: veclen
    real solar(veclen), view(veclen), ra(veclen), hb, br
    real brdf0,brdf1,brdf2
    real RL_brdf(veclen)

    real rs_thick(veclen),li_sparse(veclen),secsolar(veclen),secvia(veclen)
    real cossolar(veclen),cosvia(veclen),cosra(veclen),sinsolar(veclen),sinvia(veclen),sinra(veclen)
    real cosxi(veclen),xi(veclen),tansolar(veclen),tanvia(veclen),theta_new_v(veclen),theta_new_s(veclen)
    real d_li2(veclen),x_li(veclen),cosl(veclen),l_li(veclen),o_li(veclen)

    real tan_theta_new_v(veclen), cos_theta_new_v(veclen), sin_theta_new_v(veclen)
    real tan_theta_new_s(veclen), cos_theta_new_s(veclen), sin_theta_new_s(veclen)

    RL_brdf=0.0

!    pi = 4 * atan(1.0)
!   calculate Ross-thick kernel
    cossolar = cos(solar)
    cosvia = cos(view)
    cosra = cos(ra)
    sinsolar = sin(solar)
    sinvia = sin(view)
    sinra = sin(ra)
    cosxi = min(cossolar*cosvia+sinsolar*sinvia*cosra,1.0)
!    if (cosxi .ge.1) then
!        cosxi = 1
!    endif
    xi = acos(cosxi)
    rs_thick = ((pi/2.0-xi)*cosxi+sin(xi)) / (cossolar+cosvia) - pi/4.0

!   alculate Li-sparse
    tansolar = sinsolar/cossolar
    tanvia = sinvia/cosvia
!    theta_new_v = atan(br*tanvia)
!    theta_new_s = atan(br*tansolar)
    tan_theta_new_v = br*tanvia
    ! Trig identities
    cos_theta_new_v = 1.0 / sqrt(1.0 + tan_theta_new_v * tan_theta_new_v)
    sin_theta_new_v = tan_theta_new_v * cos_theta_new_v
    
    tan_theta_new_s = br*tansolar
    ! Trig identities
    cos_theta_new_s = 1.0 / sqrt(1.0 + tan_theta_new_s * tan_theta_new_s)
    sin_theta_new_s = tan_theta_new_s * cos_theta_new_s

!    cosxi = cos(theta_new_s)*cos(theta_new_v)+sin(theta_new_s)*sin(theta_new_v)*cosra
    cosxi = min( cos_theta_new_s*cos_theta_new_v + sin_theta_new_s*sin_theta_new_v*cosra, 1.0 )
!    if (cosxi.ge.1) then
!        cosxi = 1
!    endif
    secsolar = 1.0/cos_theta_new_s
    secvia = 1.0/cos_theta_new_v
!    d_li2 = abs(tan(theta_new_s)**2+tan(theta_new_v)**2 - 2.0*tan(theta_new_s)*tan(theta_new_v)*cosra)
    d_li2 = abs(tan_theta_new_s**2+tan_theta_new_v**2 - 2.0*tan_theta_new_s*tan_theta_new_v*cosra)
    x_li = tan_theta_new_s*tan_theta_new_v*sinra
    cosl = hb*sqrt(d_li2+x_li**2)/(secsolar+secvia)
!    if (cosl .ge. 1.0) then
!    o_li=0.0
!    else
    l_li=acos(cosl)
    o_li=(l_li-sin(l_li)*cosl)*(secsolar+secvia)/pi
    where( cosl .ge. 1.0 ) o_li = 0.0
!    endif
    li_sparse=o_li-(secsolar+secvia)+0.5*(1.0+cosxi) * secsolar*secvia
    RL_brdf=brdf0+brdf1*rs_thick+brdf2*li_sparse
    return
  end function RL_brdf

  function black_sky(brdf0,brdf1,brdf2,theta,veclen)
    integer veclen
    real brdf0,brdf1,brdf2
    real theta(veclen)
    real black_sky(veclen)

!   theta should be in radians
    black_sky = brdf0 + &
         brdf1*(-0.007574-0.070987*theta**2 + 0.307588*theta**3) + &
         brdf2*(-1.284909-0.166314*theta**2 + 0.041840*theta**3)
    return
  end function black_sky

  
  function white_sky(brdf0, brdf1, brdf2)
    real brdf0, brdf1, brdf2
    real white_sky
    white_sky = brdf0 + brdf1*0.189184 - brdf2*1.377622
    return
  end function white_sky

end module reflectance_mod
