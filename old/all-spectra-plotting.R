library(tidyverse)
library(glue)

cell_info_df <- fs::dir_ls('/Volumes/hondula/DATA/L1-reflectance/meta') %>% 
  vroom::vroom() %>% dplyr::filter(sensor_zenith != '-9999') %>%
  mutate(flightdate = lubridate::as_date(str_sub(flightline, 1, 8)))
# 221 non-NA l1 reflectance spectra

rad_info_df <- fs::dir_ls('/Volumes/hondula/DATA/L1-radiance/meta', glob = '*BLUE*', invert = TRUE) %>% vroom::vroom(delim = ',')
rad_info_df %>% dplyr::filter(gpstime_hrs != '-9999') %>% write_csv('results/radiance-gpstime.csv')

rad_info_df <- read_csv('results/radiance-gpstime.csv') %>%
  dplyr::select(siteID, nmdLctn, cellx, celly, radfile, gpstime_hrs) %>%
  mutate(flightdate = lubridate::as_date(str_sub(radfile, 1, 8))) %>%
  dplyr::filter(gpstime_hrs != '-9999') %>%
  distinct()
rad_info_df

cell_info_join <- cell_info_df %>% left_join(rad_info_df, by = c('siteID', 'nmdLctn', 'cellx', 'celly', 'flightdate'))

cell_info_join %>% write_csv('results/l1-info.csv')

cell_info_df <- read_csv('results/l1-info.csv')

# cell_info_df %>% write_csv('results/spectra-ids.csv')

mysite <- 'PRLA'

spectra_files <- fs::dir_ls('results/L1-reflectance/spectra')
barc_files <- grep(mysite, spectra_files, value = TRUE)

barc_meta <- cell_info_df %>% dplyr::filter(siteID == mysite)
meta_files <- fs::dir_ls('results/L1-reflectance/meta')
barc_meta <- grep(mysite, meta_files, value = TRUE) %>% 
  vroom::vroom(id = 'filename') %>%
  mutate(filename = tools::file_path_sans_ext(basename(filename))) %>%
  rename(flightline2 = filename)

barc_spectra <- barc_files %>% 
  purrr::map(~read_tsv(.x)) %>%
  # purrr::map(~dplyr::select(.x, -band)) %>%
  purrr::map(~tidyr::pivot_longer(.x, cols = -c(wl, band), names_to = 'loctype')) %>%
  bind_rows(.id = 'filename') %>%
  mutate(filename = tools::file_path_sans_ext(basename(filename))) %>%
  rename(flightline2 = filename) %>%
  dplyr::filter(value != -9999) %>%
  left_join(barc_meta) %>%
  dplyr::filter(!loctype %in% c('buoy.c1', 'buoy.c2'))
  
barc_spectra %>% 
  mutate(aop_yr = str_sub(flightline, 1, 4)) %>%
  ggplot(aes(x = wl, y = value/10000)) +
  geom_rect(aes(xmin = 1260, xmax = 1560, ymin = -Inf, ymax = Inf), fill = 'gray') +
  geom_rect(aes(xmin = 1760, xmax = 1960, ymin = -Inf, ymax = Inf), fill = 'gray') +
  # geom_text(aes(label = round(solar_zenith,1), x = 0, y = 0.09, vjust = "inward", hjust = "inward")) +
  geom_line(aes(col = flightline), lwd = 0.1) +
  geom_point(aes(col = flightline), alpha = 0.5, pch = 21, size = 0.5) +
  # xlim(c(NA, 900)) +
  coord_cartesian(ylim = c(0, 0.1)) +
  theme_bw() +
  facet_grid(vars(aop_yr), vars(loctype)) +
  # facet_wrap(vars(loctype)) +
  theme(legend.position = 'bottom') +
  ggtitle(glue('L1 Reflectance at {mysite} AOS pts'))

ggsave(glue('figs/{mysite}-spectra2-wv-bands.png'), width = 10, height = 6)


gg <- barc_spectra %>% 
  mutate(aop_yr = str_sub(flightline, 1, 4)) %>%
  ggplot(aes(x = wl, y = value/10000)) +
  geom_text(aes(label = round(solar_zenith,1), x = 0, y = 0.09, vjust = "inward", hjust = "inward")) +
  geom_line(aes(lty = loctype, col = clouds)) +
  # geom_point(aes(fill = loctype), alpha = 0.5, pch = 21) +
  # xlim(c(380, 900)) + 
  coord_cartesian(ylim = c(0, 0.2)) +
  theme_bw() +
  # facet_grid(vars(aop_yr), vars(flightline)) +
  facet_wrap(vars(flightline)) +
  theme(legend.position = 'bottom')

gg
ggsave(glue('figs/{mysite}-spectra.png'), gg, width = 10, height = 6)

  # ggtitle(glue('{my_aq_site} {my_aop_yr}\n {flightline}, {cloud_conditions}'))


