library(lhs)
library(cowplot)

# get measured spectra
my_aq_site <- "PRLA"
my_aop_yr <- "2017"
my_loc_type <- 'buoy.c0'
spectra_dir <- "data/spectra_buff5m"
spectra_file <- fs::dir_ls(spectra_dir, regexp = glue("{my_aq_site}_{my_aop_yr}"))
spectra_df <- read_csv(spectra_file)
my_spectra_pt <- spectra_df %>% filter(my_cellid) %>% 
  filter(wavelength < 700, wavelength > 400)
# probably need to chop off < 400 nm? 
my_spectra_pt
my_spectra_pt$reflectance %>% plot()
my_spectra_pt$wavelength
my_wls <- my_spectra_pt$wavelength

# read in water absorb, water backs, phyto* absorb
water_df <- read_csv('data/water-approx1nm.csv') %>%
  dplyr::filter(wl %in% my_wls)
wasi_df <- read_csv('data/gege-wasi-spectra.csv') %>% 
  dplyr::filter(nm %in% my_wls)

water_df$wl

absorb_water <- water_df$absorb_m1
backs_water <- water_df$backs_m1
absorb_phyto <- wasi_df$phyto_mixA
i = 1
beta <- c(abs440_cdom[i],
          slope_cdom[i],
          ref_cdom[i],
          abs440_nap[i],
          slope_nap[i],
          ref_napA[i],
          chl_ugL[i],
          backs440_parts[i],
          ref_napBb[i],
          bbp_gamma[i],
          runif(1, 0, 1)) # random error
# Forward Spectral models
y = my_spectra_pt$reflectance/pi
plot(y)
spectra_error_model <- function(beta, y){
  ## CDOM Absorption
  abs_cdom_spectra <- my_wls %>%
    purrr::map_dbl(~beta[1] * exp(-beta[2] * (.x - beta[3])))
  ## NAP Absorption
  abs_nap_spectra <- my_wls %>%
    purrr::map_dbl(~beta[4] * exp(-beta[5] * (.x - beta[6])))
  ## Phyto Absorption
  abs_phy_spectra <- absorb_phyto * beta[7]
  ## Total Absorption
  abs_tot <- absorb_water + abs_cdom_spectra + abs_nap_spectra + abs_phy_spectra
  ## Particulate backscattering
  backs_nap_spectra <- my_wls %>%
    purrr::map_dbl(~beta[8] * (beta[9]/.x)^-beta[10])
  ## Total backscattering
  backs_tot <- backs_water + backs_nap_spectra
  ## u
  u_spectra <- (backs_tot/(abs_tot + backs_tot))
  g1 = 0.949 # eqn 8.1
  g2 = 0.794 # eqn 8.1
  rrs_model_spectra <- (g1*u_spectra) + (g2 * u_spectra^2)
  # my_zeta = 0.52
  # my_gamma = 1.6 
  Rrs_model_spectra <- (0.52 * rrs_model_spectra)/(1-1.6*rrs_model_spectra) + beta[11]
  sum(((y - Rrs_model_spectra)^2)/length(my_wls))
  # max(abs(y - Rrs_model_spectra))
  # return(Rrs_model_spectra)
}

spectra_params_model <- function(beta){
  ## CDOM Absorption
  abs_cdom_spectra <- my_wls %>%
    purrr::map_dbl(~beta[1] * exp(-beta[2] * (.x - beta[3])))
  ## NAP Absorption
  abs_nap_spectra <- my_wls %>%
    purrr::map_dbl(~beta[4] * exp(-beta[5] * (.x - beta[6])))
  ## Phyto Absorption
  abs_phy_spectra <- absorb_phyto * beta[7]
  ## Total Absorption
  abs_tot <- absorb_water + abs_cdom_spectra + abs_nap_spectra + abs_phy_spectra
  ## Particulate backscattering
  backs_nap_spectra <- my_wls %>%
    purrr::map_dbl(~beta[8] * (beta[9]/.x)^-beta[10])
  ## Total backscattering
  backs_tot <- backs_water + backs_nap_spectra
  ## u
  u_spectra <- (backs_tot/(abs_tot + backs_tot))
  g1 = 0.949 # eqn 8.1
  g2 = 0.794 # eqn 8.1
  rrs_model_spectra <- (g1*u_spectra) + (g2 * u_spectra^2)
  # my_zeta = 0.52
  # my_gamma = 1.6 
  Rrs_model_spectra <- (0.52 * rrs_model_spectra)/(1-1.6*rrs_model_spectra) + beta[11]
  # sum(((y - Rrs_model_spectra)^2)/length(my_wls))
  return(Rrs_model_spectra)
  # return(Rrs_model_spectra)
}

y <- my_spectra_pt$reflectance
plot(y)
beta_init <- c(6, # cdom absorb
               0.008, # cdom slope
               440, # cdom ref wavelength
               10, # nap absorb
               0.005, # nap slope
               440, # nap ref wavelength
               30, # chl conc
               0.75, # backscatter
               555, # backscatter ref wavelength
               1, # backscatter gamma
               0.01) # error

data.frame(wl = my_wls, 
           rrs_meas = y,
           rrs_model = spectra_params_model(beta_init)) %>%
  ggplot(aes(x = wl, y = rrs_model)) +
  geom_vline(aes(xintercept = c(442)), col = 'green', lwd = 2) +
  geom_vline(aes(xintercept = c(665)), col = 'green', lwd = 2) +
  geom_vline(aes(xintercept = c(625)), col = 'darkgreen', lwd = 2, alpha = 0.5) +
  geom_line(aes(y = rrs_meas), col = 'blue') +
  geom_point(aes(y = rrs_meas), col = 'blue') +
  geom_point() +
  ylim(0, NA) +
  theme_bw()

# plot(Rrs_model_spectra)

spectra_error_model(beta, y)

# my_spectra_pt$reflectance %>% plot()

y <- my_spectra_pt$reflectance
plot(y)
beta_init <- c(10, # cdom absorb
               0.008, # cdom slope
               440, # cdom ref wavelength
               5, # nap absorb
               0.001, # nap slope
               440, # nap ref wavelength
               20, # chl conc
               0.5, # backscatter
               555, # backscatter ref wavelength
               1, # backscatter gamma
               0.01 # error
               )
my_optim <- optimx(
  beta_init, # initial guess! for betas
  spectra_error_model, # what function to optimize
  y = y, # arguments to your function
  # lower = c(0, 0, 100, 0, 0, 100, 0, 0, 100, 0, 0),
  # lower = c(),
  control = list(maxit = 10000,
                 all.methods = TRUE))  

my_optim %>% arrange(convcode)

beta_optim <- my_optim["BFGS",1:11] %>% unlist()
beta_optim <- my_optim["L-BFGS-B",1:11] %>% unlist()
beta_optim <- my_optim["Nelder-Mead",1:11] %>% unlist()

data.frame(wl = my_wls, 
           rrs_meas = y,
           rrs_model = spectra_params_model(beta_optim)) %>%
  ggplot(aes(x = wl, y = rrs_model)) +
  geom_line(col = 'black') +
  geom_line(aes(y = rrs_meas), col = 'blue') +
  geom_point(aes(y = rrs_meas), col = 'blue') +
  geom_point() +
  ylim(0, NA) +
  theme_bw()

my_optim <- optim(
  c(0, 1), # initial guess! for betas
  spectra_error_model, # what function to optimize
  y = y) # arguments to your function. 


# model_spectra <- spectra_model(5)

model_rsme <- sum(((my_spectra_pt$reflectance - model_spectra)^2)/length(my_wls))
model_rsme

data.frame(my_wls, model_spectra) %>%
  ggplot(aes(x = my_wls, y = model_spectra)) +
  geom_point() +
  geom_point(aes(y = my_spectra_pt$reflectance), col = 'red') +
  # geom_point(aes(y = model_spectra/pi), col = 'purple') +
  geom_point(aes(y = model_spectra + model_spectra^2), col = 'blue')

m1 <- purrr::map(1:100, ~spectra_model(.x))
names(m1) <- glue::glue('model{1:100}')

m1_rsme <- m1 %>% purrr::map_dbl(~sum((my_spectra_pt$reflectance - .x)^2)/length(my_wls))
names(which.min(m1_rsme))

gg <- m1 %>% 
  as.data.frame() %>%
  mutate(wl = my_wls) %>%
  tidyr::pivot_longer(cols = 1:100) %>%
  ggplot(aes(x = wl, y = value, group = name)) +
  geom_line(alpha = 0.8) +
  ylim(c(0,1)) +
  ggtitle('Modeled spectra')
gg

tm <- m1 %>% 
  as.data.frame() %>%
  mutate(wl = my_wls) %>%
  ggplot(aes(x = my_wls, y = model29)) +
  ylim(0, 0.1) +
  geom_line(col = 'green', lwd = 2)

rrs_meas <- my_spectra_pt %>%
  ggplot(aes(x = wavelength, y = reflectance/pi)) +
  geom_line(alpha = 0.5) +
  geom_smooth() +
  geom_point(pch = 21) +
  ylim(0,NA) +
  ggtitle(glue('Measured spectra at {my_aq_site} {my_aop_yr}'))
rrs_meas

plot_grid(gg, rrs_meas)
plot_grid(tm, rrs_meas)

# library(optimx)
# ?optimx
## modeled spectra

# optimx(par, fn, lower, upper)

# ggsave('figs/spectra-model.png', width = 5, height = 6)

