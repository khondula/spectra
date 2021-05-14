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

# WATER data
water3 <- read_excel("../../data/water-spectra/water-spectra-buiteveld94.xlsx") %>%
  dplyr::filter(lambda_nm >= 400)
water_wls <- water3$lambda_nm

# water3 %>% write_csv("bo-model-app/data/water.csv")
# MODELS

calc_absorb <- function(wl, absREF, aSF =  1.3e-2, wlREF = 443){
  abs_wl <- absREF * exp(-aSF * (wl - wlREF))
  return(abs_wl)
}

calc_scatter <- function(wl, scatterREF, bSF = 1, wlREF = 443){
  scatter_wl <- scatterREF * (wl/wlREF)^(bSF)
  return(scatter_wl)
}

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


# ABSORBANCE
abs_water = data.frame(wl = water3$lambda_nm,
                       abs_water = water3$absorption_m1)

ABSnap_spectra <- calc_absorb_spectra(absNAP443, aSF = sNAP) %>%
  rename(wl = wavelength, abs_nap = abs_m1) %>%
  dplyr::select(wl, abs_nap)
ABScdom_spectra <- calc_absorb_spectra(absCDOM443, aSF = sCDOM) %>%
  rename(wl = wavelength, abs_cdom = abs_m1) %>%
  dplyr::select(wl, abs_cdom)

abs_df <- abs_water %>%
  left_join(ABSnap_spectra) %>%
  left_join(ABScdom_spectra) %>%
  mutate(abs_total = abs_water + abs_nap + abs_cdom)

# BACKSCATTERING
backs_water = data.frame(wl = water3$lambda_nm,
                         backs_water = 0.5*(water3$scattering_m1))
BSchl_spectra <- calc_scatter_spectra(scatterCHL443, sbCHL, bratioCHL) %>%
  rename(wl = wavelength, bb_chl = backscatter_m1) %>%
  dplyr::select(wl, bb_chl)
BSsed_spectra <- calc_scatter_spectra(scatterSED443, sbSED, bratioSED) %>%
  rename(wl = wavelength, bb_sed = backscatter_m1) %>%
  dplyr::select(wl, bb_sed)

backs_df <- backs_water %>%
  left_join(BSchl_spectra) %>%
  left_join(BSsed_spectra) %>%
  mutate(backs_total = backs_water + bb_chl + bb_sed)

mu0 = 1
C_mu0 = -0.629 * mu0 + 0.975
surface_factor = 0.544

rs_df <- abs_df %>% left_join(backs_df) %>%
  mutate(rD = 0.544 * C_mu0 * (backs_total/ (abs_total + backs_total)))

rs1 <- rs_df %>% 
  dplyr::filter(wl %in% water_wls) %>%
  ggplot(aes(x = wl, y = rD)) +
  geom_line() +
  theme_bw() +
  ggtitle("Reflectance")
rs2 <- rs_df %>% 
  dplyr::filter(wl %in% water_wls) %>%
  ggplot(aes(x = wl, y = abs_total)) +
  geom_line(lwd = 1) +
  geom_line(aes(y = abs_nap), col = "brown") +
  geom_line(aes(y = abs_cdom), col = "orange") +
  geom_line(aes(y = abs_water), col = "blue") +
  theme_bw() +
  ggtitle("Absorbance")
rs3 <- rs_df %>% 
  dplyr::filter(wl %in% water_wls) %>%
  ggplot(aes(x = wl, y = backs_total)) +
  geom_line(lwd = 1) +
  geom_line(aes(y = bb_chl), col = "green") +
  geom_line(aes(y = bb_sed), col = "orange") +
  geom_line(aes(y = backs_water), col = "blue") +
  theme_bw() +
  ggtitle("Backscattering")

gg <- cowplot::plot_grid(rs1, rs2, rs3, 
                   ncol = 1, rel_heights = c(1, 0.8, 0.8))

gg1 <- rs1 +
  geom_point_interactive(aes(tooltip = wl, data_id = wl), size = 1) 
x <- girafe(code = print(gg1))
x

