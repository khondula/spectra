# chapter 5

# ABSORPTION
# eqn 5.1 - absorption coefficient of Non-algal particles
# abs[wl] = abs[wlREF] * exp((-SF) * [wl] - [wlREF])

# abs: absorption coefficient (m^-1)
# wl: wavelength (lambda)
# wlREF: reference wavelength
#   commonly in blue (460 nm)
# SF: slope factor
#   can be calculated using a linear regression between 
#   wavelength and natural log abs[wl] for 380-730
#   excluding 400-480 and 620-710 (pigments)

# absREF is absorbance at reference wavelength
absREF <- 50
# wl <- 500
spectra_min <- 300
spectra_max <- 800

calc_sm_absorb <- function(wl, absREF, wlREF = 400, SF = 0.0017){
  abs_wl <- absREF * exp(-SF * (wl - wlREF))
  return(abs_wl)
}

calc_sm_spectra <- function(absREF, 
                              spectra_id, 
                              spectra_min = 300, 
                              spectra_max = 800){
  wl_range <- spectra_min:spectra_max
  sm_spectra <- wl_range %>% 
    purrr::map_dbl(~calc_sm_absorb(.x, absREF)) %>%
    as_tibble() %>%
    mutate(wavelength = wl_range) %>%
    mutate(wavelength = as.numeric(wavelength)) %>%
    mutate(spectra_id = spectra_id)
  return(sm_spectra)
}

my_abs <- c(0.1, 1, 10, 20, 50)
my_ids <- c('abs01', 'abs1', 'abs10', 'abs20', 'abs50')

my_spectra <- purrr::map2_dfr(my_abs, my_ids, ~calc_sm_spectra(.x, .y))

my_spectra %>%
  ggplot(aes(x = wavelength, y = value)) +
  geom_line(aes(col = spectra_id)) + 
  theme_bw()

# BACKSCATTERING

# eqn 5.2
# wl - wavelength
# backsREF - backscattering at reference wavelength 
# wlREF - reference wavelength
# gammabb - power law exponent of spectral slope of scattering coefficient
#   low values ~ 0.3 are larger particles
#   high values ~ 1.7 are smaller particles
#   also an indicator of composition
# returns backs_wl - backscattering coeff at given wlb


calc_sm_backs <- function(wl, backsREF, wlREF = 400, gammabb = 1){
  backs_wl <- backsREF * (backsREF/wl)^(-gammabb)
  return(backs_wl)
}

calc_sm_spectra_backs <- function(backsREF, 
                            spectra_id, 
                            spectra_min = 300, 
                            spectra_max = 800){
  wl_range <- spectra_min:spectra_max
  sm_spectra <- wl_range %>% 
    purrr::map_dbl(~calc_sm_backs(.x, backsREF)) %>%
    as_tibble() %>%
    mutate(wavelength = wl_range) %>%
    mutate(wavelength = as.numeric(wavelength)) %>%
    mutate(spectra_id = spectra_id)
  return(sm_spectra)
}



## Equation 2.42 and 2.43 in Gege ch 2

backs_ratio = 0.03
bstarTSM = 1
wlREF = 555
angstrom = 0.1

calc_backs_wl <- function(wl, bstarTSM = 1, wlREF = 555, angstrom = 0.1){
  backs_wl <- bstarTSM * (wl/wlREF)^(-1 *angstrom)
  return(backs_wl)
}
calc_backs_wl(500)

calc_tsm_bspectra <- function(spectra_min = 300, 
                              spectra_max = 800,
                              bstarTSM = 1, 
                              wlREF = 555, 
                              angstrom = 0.1,
                              backs_ratio = 0.03){
  wl_range <- spectra_min:spectra_max
  my_spectra <- wl_range %>% 
    purrr::map_dbl(~calc_backs_wl(.x, bstarTSM, wlREF, angstrom)) %>%
    as_tibble() %>%
    mutate(wavelength = wl_range) %>%
    mutate(wavelength = as.numeric(wavelength))
  return(my_spectra)
}

calc_tsm_bspectra(angstrom = 0.01) %>%
  ggplot(aes(x = wavelength, y = value)) +
  geom_line()
