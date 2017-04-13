real function black_sky(brdf0,brdf1,brdf2,theta)
      real brdf0,brdf1,brdf2
      real theta

!   theta should be in radians
    black_sky = brdf0 + &
        brdf1*(-0.007574-0.070987*theta**2 + 0.307588*theta**3) + &
        brdf2*(-1.284909-0.166314*theta**2 + 0.041840*theta**3)
    return
end
