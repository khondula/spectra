library(tidyverse)
library(ggplot2)
library(colorspace)
library(glue)

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
  mutate(siteyear = glue('{site}-{year}')) %>%
  mutate(location_type = substr(nmdLctn, 10, nchar(nmdLctn))) %>%
  ggplot(aes(wavelength, reflectance)) +
  geom_line(lwd = 0.8, aes(col = location_type)) +
  # scale_color_discrete_diverging("") +
  coord_cartesian(xlim = c(350, 1000), ylim = c(0, 0.1)) +
  facet_wrap(vars(siteyear)) +
  theme_bw()

ggsave('figs/spectra.png')

