real function white_sky(brdf0, brdf1, brdf2)
    real brdf0, brdf1, brdf2
    white_sky = brdf0 + brdf1*0.189184 - brdf2*1.377622
    return
end
