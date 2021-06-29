library(tidyverse)
library(glue)

# copy over meta and spectra from Volumes to results for reading in
cell_info_df <- fs::dir_ls('results/L1-reflectance2/meta') %>% 
  vroom::vroom() %>% dplyr::filter(sensor_zenith != '-9999')

cell_info_df <- cell_info_df %>% 
  dplyr::mutate(flightdate = purrr::map_chr(reflectfile, ~str_split(.x, "_", simplify = TRUE)[5])) %>%
  dplyr::mutate(flightdate = lubridate::as_date(flightdate))

# cell_info_df$reflectfile %>% purrr::map_chr(~str_split(.x, "_", simplify = TRUE)[5])

nrow(cell_info_df)
cell_info_df %>% pull(reflectfile) %>% unique() %>% length()
# 282 non-NA l1 reflectance spectra

radnames <- fs::dir_ls('results/L1-radiance/meta', glob = '*BLUE*', invert = TRUE) %>%
  purrr::map(~names(read_csv(.x)))
# radnames[18]

rad_info_df <- fs::dir_ls('results/L1-radiance/meta', glob = '*BLUE*', invert = TRUE) %>%
  vroom::vroom(delim = ',')
rad_info_df %>%
  dplyr::filter(gpstime_hrs != '-9999') %>%
  write_csv('results/radiance-gpstime.csv')
# rad_info_df %>%
#   dplyr::filter(gpstime_hrs != '-9999') %>%
#   dplyr::select(-loctype, -nmdLctn) %>%
#   distinct() %>% nrow()
# rad_info_df %>%
#   dplyr::filter(gpstime_hrs != '-9999') %>%
#   pull(radfile) %>% unique() %>% length()
# 
# rad_info_df$radfile[1]
rad_info_df <- read_csv('results/radiance-gpstime.csv') %>%
  mutate(flightdate = lubridate::as_date(str_sub(radfile, 1, 8))) %>%
  dplyr::filter(gpstime_hrs != '-9999') %>%
  dplyr::select(siteID, cellid, cellrow, cellcol, radfile, gpstime_hrs, flightdate) %>%
  distinct()
rad_info_df
# 258 non-NA L1 radiance data (no duplicates for buoy depths)
names(cell_info_df)
names(rad_info_df)

cell_info_join <- cell_info_df %>% 
  left_join(rad_info_df, by = c('siteID', 'cellid', 'cellrow', 'cellcol', 'flightdate'))

cell_info_join %>% write_csv('results/l1-info.csv')

cell_info_join <- read_csv('results/l1-info.csv')

# cell_info_df %>% write_csv('results/spectra-ids.csv')

mysite <- 'TOMB'

spectra_files <- fs::dir_ls('results/L1-reflectance2/spectra')
mysite_spectra_files <- grep(mysite, spectra_files, value = TRUE)



mysite_spectra <- mysite_spectra_files %>% 
  purrr::map(~read_tsv(.x)) %>%
  # purrr::map(~dplyr::select(.x, -band)) %>%
  purrr::map(~tidyr::pivot_longer(.x, cols = -c(wl, band), names_to = 'loctype')) %>%
  bind_rows(.id = 'filename') %>%
  mutate(filename = tools::file_path_sans_ext(basename(filename))) %>%
  mutate(flightline = str_sub(filename, 24, 38)) %>%
  dplyr::filter(value != -9999) %>%
  # left_join(barc_meta) %>%
  dplyr::filter(!loctype %in% c('buoy.c1', 'buoy.c2'))
  
mysite_spectra %>% 
  mutate(aop_yr = str_sub(flightline, 1, 4)) %>%
  ggplot(aes(x = wl, y = value/10000)) +
  geom_rect(aes(xmin = 1260, xmax = 1560, ymin = -Inf, ymax = Inf), fill = 'gray') +
  geom_rect(aes(xmin = 1760, xmax = 1960, ymin = -Inf, ymax = Inf), fill = 'gray') +
  # geom_text(aes(label = round(solar_zenith,1), x = 0, y = 0.09, vjust = "inward", hjust = "inward")) +
  geom_line(aes(col = flightline), lwd = 0.5) +
  geom_point(aes(col = flightline), alpha = 0.5, pch = 21, size = 0.5) +
  xlim(c(NA, 900)) +
  coord_cartesian(ylim = c(0, 0.1)) +
  theme_bw() +
  facet_grid(vars(aop_yr), vars(loctype)) +
  # facet_wrap(vars(loctype)) +
  # theme(legend.position = 'bottom') +
  ggtitle(glue('L1 Reflectance at {mysite} AOS pts'))

ggsave(glue('figs/{mysite}-spectra2-bands.png'), width = 6, height = 7)
# ggsave(glue('figs/{mysite}-spectra2-wv-bands.png'), width = 10, height = 7)


library(ggiraph)

gg <- mysite_spectra %>% 
  mutate(aop_yr = str_sub(flightline, 1, 4)) %>%
  ggplot(aes(x = wl, y = value/10000)) +
  geom_rect(aes(xmin = 1260, xmax = 1560, ymin = -Inf, ymax = Inf), fill = 'gray') +
  geom_rect(aes(xmin = 1760, xmax = 1960, ymin = -Inf, ymax = Inf), fill = 'gray') +
  # geom_text(aes(label = round(solar_zenith,1), x = 0, y = 0.09, vjust = "inward", hjust = "inward")) +
  geom_line(aes(col = flightline), lwd = 0.5) +
  geom_point(aes(col = flightline), alpha = 0.5, pch = 21, size = 0.5) +
  xlim(c(NA, 900)) +
  coord_cartesian(ylim = c(0, 0.15)) +
  theme_bw() +
  facet_grid(vars(aop_yr), vars(loctype)) +
  # facet_wrap(vars(loctype)) +
  theme(legend.position = 'bottom') +
  ggtitle(glue('L1 Reflectance at {mysite} AOS pts'))

gg2 <- gg + geom_point_interactive(aes(tooltip = flightline, data_id = flightline, col = flightline), size = 0.001)
girafe(code = print(gg2))

site_meta <- cell_info_join %>% dplyr::filter(siteID == mysite) %>% 
  mutate(aop_yr = lubridate::year(flightdate)) %>%
  mutate(flightline = str_sub(reflectfile, 19, 33))
# get metadata for a spectra

site_meta %>% 
  dplyr::filter(str_detect(nmdLctn, "AOS.buoy.c0")) %>%
  dplyr::select(nmdLctn, flightline, sensor_zenith, solar_zenith, clouds, flightdate, gpstime_hrs, aop_yr)
