# absorbance and backscattering of water
library(readxl)
library(ggplot2)

popefry <- read_excel('../../data/water-spectra/water-spectra-popefry97.xlsx') %>%
  rename(wavelength_nm = 1,
         abs_m1 = 2,
         abs_sd_m1 = 3,
         percent = 4)

head(popefry)

popefry %>%
  ggplot(aes(x = wavelength, y = abs_m1)) +
  geom_line()

water2 <- read_excel('../../data/water-spectra/water-spectra-bukata.xlsx')

water2 %>%
  ggplot(aes(x = wavelength_nm, y = abs_SB)) +
  geom_line() +
  geom_line(data = popefry, aes(y = abs_m1), col = 'red') +
  geom_line(aes(y = abs_MA), col = 'blue') +
  theme_bw() +
  ggtitle('Absorbance water')

water2 %>%
  ggplot(aes(x = wavelength_nm, y = backs_SB)) +
  geom_line() +
  theme_bw() +
  ggtitle('Backscattering water')
