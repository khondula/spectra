
# FUNCTIONS FOR SHINY APP
library(tidyverse)

calc_absorb <- function(wl, absREF, aSF =  1.3e-2, wlREF = 443){
  abs_wl <- absREF * exp(-aSF * (wl - wlREF))
  return(abs_wl)
}

calc_scatter <- function(wl, scatterREF, bSF, wlREF = 443){
  scatter_wl <- scatterREF * (wl/wlREF)^(bSF)
  return(scatter_wl)
}

calc_scatter_spectra <- function(scatterREF, bSF,
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

calc_absorb_spectra <- function(absREF,  aSF =  1.3e-2, wlREF = 443, 
                                spectra_min = 400, spectra_max = 800){
  wl_range <- spectra_min:spectra_max
  my_spectra <- wl_range %>% 
    purrr::map_dbl(~calc_absorb(.x, absREF, aSF, wlREF)) %>%
    as_tibble() %>%
    rename(abs_m1 = value) %>%
    mutate(wavelength = wl_range) %>%
    mutate(wavelength = as.numeric(wavelength)) %>%
    mutate(wlREF = wlREF, SF = aSF)
  return(my_spectra)
}