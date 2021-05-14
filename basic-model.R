# INPUTS

absCHL443 <- 0.143
absNAP443 <- 0.139
absCDOM443 <- 5.01e-1

sNAP <- 1.3e-2
sCDOM <- 1.7e-2

scatterCHL443 <- 1.87
scatterSED443 <- 1.85
bbCHL <- 1.73e-2
bbSED <- 4.64e-2

bratioCHL <- bbCHL/scatterCHL443
bratioSED <- bbSED/scatterSED443

# sbPART <- 1
sbCHL <- 1
sbSED <- 1

sbPART <- (scatterCHL443 * sbCHL + scatterSED443 * sbSED)/(scatterCHL443 + scatterSED443)
sbPART

# MODELS

calc_absorb <- function(wl, absREF, aSF =  1.3e-2, wlREF = 443){
  abs_wl <- absREF * exp(-SF * (wl - wlREF))
  return(abs_wl)
}

calc_scatter <- function(wl, scatterREF, bSF = 1, wlREF = 443){
  scatter_wl <- scatterREF * (wl/wlREF)^(bSF)
  return(scatter_wl)
}

# calc_scatter(400, scatterCHL443, 1)

calc_scatter_spectra <- function(scatterREF, bSF =  1,
                                 bratio = 0.01, wlREF = 443, 
                                spectra_min = 400, spectra_max = 800){
  wl_range <- spectra_min:spectra_max
  my_spectra <- wl_range %>% 
    purrr::map_dbl(~calc_scatter(.x, scatterREF, bSF, wlREF)) %>%
    as_tibble() %>%
    mutate(backscatter_m1 = value*bratio) %>%
    mutate(wavelength = wl_range) %>%
    mutate(wavelength = as.numeric(wavelength)) %>%
    mutate(wlREF = wlREF, SF = bSF)
  return(my_spectra)
}

calc_absorb_spectra <- function(absREF,  =  1.3e-2, wlREF = 443, 
                         spectra_min = 400, spectra_max = 800){
  wl_range <- spectra_min:spectra_max
  my_spectra <- wl_range %>% 
    purrr::map_dbl(~calc_absorb(.x, absREF, aSF, wlREF)) %>%
    as_tibble() %>%
    rename(abs_m1 = value) %>%
    mutate(wavelength = wl_range) %>%
    mutate(wavelength = as.numeric(wavelength)) %>%
    mutate(wlREF = wlREF, SF = SF)
  return(my_spectra)
}



ABSnap_spectra <- calc_absorb_spectra(absNAP443, aSF = sNAP)
ABScdom_spectra <- calc_absorb_spectra(absCDOM443, bSF = sCDOM)

BSchl_spectra <- calc_scatter_spectra(scatterCHL443, sbCHL, bratioCHL)
BSsed_spectra <- calc_scatter_spectra(scatterSED443, sbSED, bratioSED)



C_mu0 = -0.629 * mu0 + 0.975
surface_factor = 0.544
