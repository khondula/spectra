# putting the model all together

# reflectance just BELOW the water surface
# rs[wl] = C * ( backs[wl] / (abs[wl] + backs[wl]))

# C is a function of the solar zenith angle
# C = -0.629 * mu0 + 0.975
# mu0 is the cosine of the zenith angle of refracted photons
# sun at zenith C = 0.346

# reflectance just ABOVE water surface
# 0.544 factor for surface to just below surface
# rs[wl] = 0.544 * C * ( backs[wl] / (abs[wl] + backs[wl]))

calc_c_zenith <- function(zenith_angle){
  c_zenith <- -0.629 * cos(zenith_angle) + 0.975
  return(c_zenith)
}
# cdom absorbence
calc_cdom_absorb <- function(wl, cdom_absREF, cdom_wlREF = 400, SF_cdom = 0.0017){
  abs_wl <- cdom_absREF * exp(-SF_cdom * (wl - cdom_wlREF))
  return(abs_wl)
}
# specific absorbance of suspended matter
calc_sm_absorb_spec <- function(wl, sm_absREF, sm_wlREF = 400, SF_sm = 0.0017){
  abs_wl <- sm_absREF * exp(-SF_sm * (wl - sm_wlREF))
  return(abs_wl)
}
# specific absorbance of chl a
# calc_chla_absorb_spec <- function(wl, )

# total absorption at wavelength
calc_abs_wl <- function(wl, 
                        h2o_absorb,
                        cdom_absorb, 
                        sm_absorb_spec, sm_conc,
                        chl_absorb_spec, chla_conc){
  abs_wl <- h2o_absorb + cdom_absorb + 
    (sm_absorb_spec * sm_conc) + (chl_absorb_spec * chla_conc)
  return(abs_wl)
}

# total backscattering at wavelenggth
calc_backs_wl <- function(wl, 
                        h2o_backs,
                        sm_backs_spec, sm_conc,
                        chl_backs_spec, chla_conc){
  backs_wl <- h2o_backs + (sm_backs_spec * sm_conc) + (chl_backs_spec * chla_conc)
  return(backs_wl)
}

calc_rs_wl <- function(wl, c_zenith, abs_wl, backs_wl){
  rs_wl <- 0.544 * c_zenith * (backs_wl / (abs_wl + backs_wl))
  return(rs_wl)
}

