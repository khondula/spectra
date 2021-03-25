# extract spectra from pixels within
# a buffered distance of a point

library(neonhs)
library(viridis)
library(raster)
library(hdf5r)
library(sf)
library(tidyverse)
library(glue)
library(terra)
library(stars)
library(scico)
source('R/myfxns.R')

# set up
spectra_dir <- 'H:/DATA/spectra_buff5m'
shp_dir <- 'H:/DATA/spatial'

# inputs
my_aq_site <- 'BARC'
my_aop_yr <- '2018'
my_loc_type <- 'buoy.c0'
nbands <- 100


save_spectra_buff5m <- function(my_aq_site, my_aop_yr,
                                my_loc_type, nbands, max_y_refl = NA){
  
  # set up
  spectra_dir <- 'H:/DATA/spectra_buff5m'
  shp_dir <- 'H:/DATA/spatial'
  # to get AOP name for AQ site
my_aop_site <- 'results/sites_join_aop_dates.csv' %>%
  readr::read_csv(col_types = 'ccccccccddD') %>%
  dplyr::filter(siteID %in% my_aq_site) %>%
  dplyr::pull(aop_site_id) %>% unique()

# read in point for AQ site
my_aq_sf <- glue('{shp_dir}/swchem_sites/swchem_sites.shp') %>%
  st_read() %>% dplyr::filter(siteID == my_aq_site) %>%
  mutate(location_type = substr(nmdLctn, 10, nchar(nmdLctn))) %>% 
  dplyr::filter(location_type == my_loc_type)
my_domain <- my_aq_sf$domanID[1]

# folder with AOP data
my_site_dir <- glue('D:/{my_aop_yr}/FullSite/{my_domain}') %>%
  fs::dir_ls(glob = glue("*{my_aop_site}*"), type = 'directory')
my_site_files <- glue('{my_site_dir}/L3/Spectrometer/Reflectance') %>%
  fs::dir_ls(glob = '*.h5')

# get the projection info from the first file
my_h5_file <- my_site_files[1]
my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
epsg_path <- glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/EPSG Code')
my_epsg <- my_h5[[epsg_path]]$read() %>% as.integer()

# project and buffer point
my_aq_prj <- my_aq_sf %>% st_transform(crs = my_epsg)
my_aq_prj_buff5 <- st_buffer(my_aq_prj, 5)
# convert to sp objects
my_aq_prj_sp <- as(my_aq_prj, "Spatial")
my_buff5_sp <- as(my_aq_prj_buff5, "Spatial")
# get corners of extent
bbox_pts_sf <- my_aq_prj_buff5 %>% 
  st_bbox() %>% st_as_sfc() %>% 
  st_as_sf() %>% st_cast(to = "POINT")
my_xys <- sf::st_coordinates(bbox_pts_sf) %>% as.data.frame()

# find the images that span the entire buffered point
my_site_files_df <- my_site_files %>%
  as.data.frame() %>% 
  rename(fullname = 1) %>%
  mutate(filename = basename(fullname)) %>%
  mutate(easting = as.numeric(str_sub(filename, 19, 24)),
         easting2 = easting + 1000,
         northing = as.numeric(str_sub(filename, 26, 32)),
         northing2 = northing + 1000)

find_which_file <- function(my_id){
  my_x <- my_xys[['X']][my_id]
  my_y <- my_xys[['Y']][my_id]
  my_h5_file <- my_site_files_df %>%
    mutate(x_in_range = {{ my_x }} >= easting & {{ my_x }} <= easting2,
           y_in_range = {{ my_y }} >= northing & {{ my_y }} <= northing2) %>%
    dplyr::filter(x_in_range, y_in_range) %>% 
    pull(fullname) %>% unique()
  return(my_h5_file)
}

my_h5_files <- 1:nrow(my_xys) %>% 
  purrr::map(~find_which_file(.x)) %>% 
  unlist() %>% unique()

# read in and merge h5 files using hs_read
if(length(my_h5_files) == 1){
  my_h5 <- hs_read(my_h5_files, bands = 1:nbands)
  names(my_h5) <- hs_wavelength(my_h5_files, bands = 1:nbands)
}
if(length(my_h5_files) > 1){
  my_h5 <- hs_read(my_h5_files[1], bands = 1:nbands)
  names(my_h5) <- hs_wavelength(my_h5_files, bands = 1:nbands)
  my_hs_list <- purrr::map(my_h5_files, ~hs_read(.x, bands = 1:nbands))
  for(i in 2:length(my_h5_files)){
    my_h5 <- merge(my_h5, my_hs_list[[i]])
  }
}

# extract spectra from pixels within the 5m buffered point
my_pt_spatvec <- my_aq_prj_sp %>% as('SpatVector')
my_pts_spatvec <- my_buff5_sp %>% as('SpatVector')

my_hs_terra <- my_h5 %>% terra::rast()

# find which cell the exact coords of sampling
my_cellid <- cells(my_hs_terra, my_pt_spatvec)[2]

my_spectra <- terra::extract(my_hs_terra, my_pts_spatvec, cells = TRUE)
my_spectra_df <- my_spectra %>% 
  as_tibble() %>%
  mutate(x = terra::xFromCell(my_hs_terra, cell),
         y = terra::yFromCell(my_hs_terra, cell)) %>%
  group_by(cell) %>%
  pivot_longer(cols = starts_with('band'), 
               names_to = 'band', 
               values_to = 'reflectance') %>% 
  tidyr::separate(band, into = c('index', 'wavelength'), sep = "_") %>%
  dplyr::mutate(wavelength = parse_number(wavelength)) %>%
  dplyr::mutate(my_cellid = cell == my_cellid)

my_spectra_pt <- dplyr::filter(my_spectra_df, my_cellid)

my_spectra_df %>%
  ggplot(aes(wavelength, reflectance)) +
  geom_line(aes(group = cell)) +
  geom_line(data = my_spectra_pt, col = 'green', lwd = 1) +
  theme_bw() +
  coord_cartesian(ylim = c(0, max_y_refl)) +
  theme(legend.position = 'none') +
  ggtitle(glue('{my_aq_site}-{my_aop_yr}-5mbuff-{str_replace_all(my_loc_type, "[:punct:]", "")}'))

ggsave(glue('figs/buff5m-spectra/{my_aq_site}-{my_aop_yr}-5mbuff-{str_replace_all(my_loc_type, "[:punct:]", "")}.png'))

my_spectra_df %>% 
  write_csv(glue('{spectra_dir}/{my_aq_site}_{my_aop_yr}_{str_replace_all(my_loc_type, "[:punct:]", "")}-buff5m.csv'))


bands_to_plot <- seq(1, 100, length.out = 12)
my_h5 <- hs_read(my_h5_files, bands = bands_to_plot)
names(my_h5) <- hs_wavelength(my_h5_files, bands = bands_to_plot)

my_aq_prj_buff <- st_buffer(my_aq_prj, 5)
my_buff_sp <- as(my_aq_prj_buff, "Spatial")
my_aq_prj_buff30 <- st_buffer(my_aq_prj, 25)
my_buff30_sp <- as(my_aq_prj_buff30, "Spatial")

my_hs_crop <- raster::crop(my_h5, my_buff30_sp)
hs_crop_stars <- my_hs_crop %>% st_as_stars()

gg <- ggplot() +
  geom_stars(data = hs_crop_stars) +
  coord_equal() +
  facet_wrap(~band) +
  geom_sf(data = my_aq_prj, pch = 3, col = 'green') +
  geom_sf(data = my_aq_prj_buff, fill = NA, col = 'white') +
  theme_void() + 
  theme(legend.title = element_blank()) +
  scale_fill_scico(palette = 'davos') +
  ggtitle(glue('{my_aq_site} {my_aop_yr} {my_loc_type}'))

ggsave(glue('figs/buff5m-maps/{my_aq_site}_{my_aop_yr}_{str_replace_all(my_loc_type, "[:punct:]", "")}.pdf'), plot = gg, width = 10, height = 8)
}

save_spectra_buff5m('BARC', '2018', 'buoy.c0', 100)
save_spectra_buff5m('BARC', '2019', 'buoy.c0', 100)
#
save_spectra_buff5m('SUGG', '2017', 'buoy.c0', 100)
save_spectra_buff5m('SUGG', '2018', 'buoy.c0', 100)
save_spectra_buff5m('SUGG', '2019', 'buoy.c0', 100)
#
save_spectra_buff5m('CRAM', '2017', 'buoy.c0', 100)
save_spectra_buff5m('CRAM', '2019', 'buoy.c0', 100)
#
save_spectra_buff5m('PRLA', '2017', 'buoy.c0', 100)
save_spectra_buff5m('PRLA', '2019', 'buoy.c0', 100)
save_spectra_buff5m('PRLA', '2020', 'buoy.c0', 100)

save_spectra_buff5m('PRPO', '2017', 'buoy.c0', 100)
save_spectra_buff5m('PRPO', '2019', 'buoy.c0', 100)

save_spectra_buff5m('LIRO', '2017', 'buoy.c0', 100)
save_spectra_buff5m('LIRO', '2020', 'buoy.c0', 100)

save_spectra_buff5m('TOOK', '2017', 'buoy.c0', 100, 0.1)
save_spectra_buff5m('TOOK', '2018', 'buoy.c0', 100, 0.1)
#save_spectra_buff5m('TOOK', '2019', 'buoy.c0', 100)

save_spectra_buff5m('FLNT', '2017', 'buoy.c0', 100)

save_spectra_buff5m('BLWA', '2017', 'buoy.c0', 100)
save_spectra_buff5m('BLWA', '2018', 'buoy.c0', 100)
save_spectra_buff5m('BLWA', '2019', 'buoy.c0', 100, 0.1)

save_spectra_buff5m('TOMB', '2017', 'buoy.c0', 100)
save_spectra_buff5m('TOMB', '2018', 'buoy.c0', 100)
save_spectra_buff5m('TOMB', '2019', 'buoy.c0', 100, 0.1)

save_spectra_buff5m('BLUE', '2017', 'buoy.c0', 100)


save_spectra_buff5m('SYCA', '2019', 'S2', 100)

save_spectra_buff5m('ARIK', '2017', 'S2', 100)
save_spectra_buff5m('ARIK', '2020', 'S2', 100)

save_spectra_buff5m('CARI', '2017', 'S2', 100)
save_spectra_buff5m('CARI', '2018', 'S2', 100)
save_spectra_buff5m('CARI', '2019', 'S2', 100)

save_spectra_buff5m('TECR', '2019', 'S2', 100)

save_spectra_buff5m('MART', '2017', 'S2', 100)
save_spectra_buff5m('MART', '2018', 'S2', 100)
save_spectra_buff5m('MART', '2019', 'S2', 100)
