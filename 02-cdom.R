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
                              spectra_max = 800){
  wl_range <- spectra_min:spectra_max
  cdom_spectra <- wl_range %>% 
    purrr::map_dbl(~calc_cdom_absorb(.x, absREF)) %>%
    as_tibble() %>%
    mutate(wavelength = wl_range) %>%
    mutate(wavelength = as.numeric(wavelength)) %>%
    mutate(spectra_id = spectra_id)
  return(cdom_spectra)
}

my_abs <- c(0.1, 1, 10, 20, 50)
my_ids <- c('abs01', 'abs1', 'abs10', 'abs20', 'abs50')

my_spectra <- purrr::map2_dfr(my_abs, my_ids, ~calc_cdom_spectra(.x, .y))

my_spectra %>%
  ggplot(aes(x = wavelength, y = value)) +
  geom_line(aes(col = spectra_id)) + 
  theme_bw()

