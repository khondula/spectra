# putting the model all together

# reflectance just below the water surface
# rs[wl] = C * ( backs[wl] / (abs[wl] + backs[wl]))

# C is a function of the solar zenith angle
# C = -0.629 * mu0 + 0.975
# mu0 is the cosine of the zenith angle of refracted photons
# sun at zenith C = 0.346

# reflectance just above water surface
# 0.544 factor for surface to just below surface
# rs[wl] = 0.544 * C * ( backs[wl] / (abs[wl] + backs[wl]))
