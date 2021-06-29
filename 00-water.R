# absorbance and backscattering of water
library(readxl)
library(ggplot2)
library(tidyverse)

popefry <- readxl::read_excel('../../data/water-spectra/water-spectra-popefry97.xlsx') %>%
  rename(wavelength_nm = 1,
         abs_m1 = 2,
         abs_sd_m1 = 3,
         percent = 4)

head(popefry)

popefry %>%
  ggplot(aes(x = wavelength_nm, y = abs_m1)) +
  geom_line()

water2 <- readxl::read_excel('../../data/water-spectra/water-spectra-bukata.xlsx')

water3 <- readxl::read_excel('../../data/water-spectra/water-spectra-buiteveld94.xlsx')

abs_water1 <- water2 %>%
  ggplot(aes(x = wavelength_nm, y = abs_SB)) +
  geom_vline(aes(xintercept = 700), col = 'green') +
  geom_vline(aes(xintercept = 760), col = 'green') +
  geom_line() +
  geom_line(data = popefry, aes(y = abs_m1), col = 'red') +
  geom_line(aes(y = abs_MA), col = 'blue') +
  geom_line(data = water3, aes(x = lambda_nm, y = absorption_m1), col = 'purple') +
  theme_bw() +
  # scale_x_continuous(breaks = seq(0,800, 50)) +
  # scale_y_continuous(breaks = seq(0, 3, 0.5)) +
  ylab(expression(absorbance~(m^-1)))
  # scale_y_log10()
  # ggtitle('Absorbance water')
cowplot::plot_grid(abs_water1, abs_water1 + scale_y_log10() + ylab(expression(log~abs~m^-1)))
ggsave('figs/water-absorbance.png', width = 6, height = 2.5)

water2 %>%
  ggplot(aes(x = wavelength_nm, y = backs_SB)) +
  geom_line() +
  geom_line(data = water3, aes(x = lambda_nm, y = scattering_m1), col = 'purple') +
  theme_bw() +
  ggtitle('Backscattering water')
