library(tidyverse)
library(ggplot2)
library(colorspace)

spectra_dir <- '/Volumes/hondula/DATA/spectra/'
my_coltypes <- 'cccccdd'
my_vals_df <- spectra_dir %>%
  fs::dir_ls() %>%
  map_dfr(~read_csv(.x, col_types = my_coltypes), 
          .id = 'filename') %>%
  mutate(filename = basename(filename)) %>%
  mutate(filename = tools::file_path_sans_ext(filename)) %>%
  tidyr::separate(filename, into = c('site', 'year'), sep = "_")
  
 
my_vals_df %>% 
  # dplyr::filter(reflectance > 0) %>%
  ggplot(aes(wavelength, reflectance)) +
  geom_line(lwd = 1, aes(col = site)) +
  theme_bw() +
  scale_color_discrete_diverging("Tropic")
