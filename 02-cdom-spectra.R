# chapter 4

# eqn 4.1 - absorption coefficient of CDOM
# abs[wl] = abs[wlREF] * exp((-SF) * [wl] - [wlREF])

# abs: absorption coefficient (m^-1)
# wl: wavelength (lambda)
# wlREF: reference wavelength
#   commonly 380, 400, 420
#   or 412, 443 (central wavelength of satellite bands)
# SF: slope factor
#   **indicator of MW, source, etc**

# abs at 400 nm range 0.64 to 63.05 m-1 in estonian lakes
# alpine lakes abs 400 around 0.1 m-1

# slope factor in lakes usually 0.017-0.018 m-1
# varies between 0.006 and 0.036 m-1
library(tidyverse)
# SUVA 254 = abs[254]/DOC
# absREF is absorbance at reference wavelength
absREF <- 50
# wl <- 500
spectra_min <- 300
spectra_max <- 800

calc_cdom_absorb <- function(wl, absREF, wlREF = 400, SF = 0.0017){
  abs_wl <- absREF * exp(-SF * (wl - wlREF))
  return(abs_wl)
}

calc_cdom_spectra <- function(absREF, 
                              spectra_id, 
                              spectra_min = 300, 
                              spectra_max = 800,
                              wlREF = 400, SF = 0.0017){
  wl_range <- spectra_min:spectra_max
  cdom_spectra <- wl_range %>% 
    purrr::map_dbl(~calc_cdom_absorb(.x, absREF, wlREF, SF)) %>%
    as_tibble() %>%
    mutate(wavelength = wl_range) %>%
    mutate(wavelength = as.numeric(wavelength)) %>%
    mutate(spectra_id = spectra_id, wlREF = wlREF, SF = SF)
  return(cdom_spectra)
}

my_suva254 = my_samps$SUVA254$analyteConcentration
my_suva280 = my_samps$SUVA280$analyteConcentration
my_doc = my_samps$doc$analyteConcentration

abs254 = my_suva254 * my_doc
abs254
abs280 = my_suva280 * my_doc
abs280

# abs[wl] = abs[wlREF] * exp((-SF) * [wl] - [wlREF])

my_sf = log(abs280/abs254)/(-(280-254))

my_spectra <- calc_cdom_spectra(abs254, 
                                spectra_min = 250,
                                spectra_id = 'bla2017', 
                                wlREF = 254, SF = my_sf)
my_spectra %>%
  mutate(wlREF = forcats::as_factor(wlREF)) %>%
  ggplot(aes(x = wavelength, y = value, group = wlREF)) +
  geom_line(aes(col = wlREF)) + 
  geom_point(aes(x = 254, y = abs254)) +
  geom_point(aes(x = 280, y = abs280)) +
  theme_bw()

# Find slope factor based on 254 and 280 absorbence values


my_abs <- c(0.1, 1, 10, 20, 50)
my_ids <- c('abs01', 'abs1', 'abs10', 'abs20', 'abs50')

my_spectra <- purrr::map2_dfr(my_abs, my_ids, ~calc_cdom_spectra(.x, .y))


