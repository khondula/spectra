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
library(vroom)
source('R/myfxns.R')

# set up
spectra_dir <- 'H:/DATA/spectra_site/'

my_water_sf <- st_read('H:/DATA/AOP/BARC_2018.shp')
nbands <- 100
my_aq_site <- 'BARC'
my_aop_yr <- '2017'
my_aop_site <- 'OSBS'
my_domain <- 'D03'


# save_spectra <- function(my_aq_site, my_aop_yr,
#                                 my_loc_type, nbands, max_y_refl = NA){
  
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

# project polygon to AOP data
my_aq_prj <- my_water_sf %>% st_transform(crs = my_epsg)
# convert to sp objects
my_aq_prj_sp <- as(my_aq_prj, "Spatial")
# get corners of extent
bbox_pts_sf <- my_aq_prj_sp %>% 
  st_bbox() %>% st_as_sfc() %>% 
  st_as_sf() %>% st_cast(to = "POINT")
my_xys <- sf::st_coordinates(bbox_pts_sf) %>% as.data.frame()

# find the images that span the entire polygon
# first get coords from filenames
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
  names(my_h5) <- hs_wavelength(my_h5_files[1], bands = 1:nbands)
  my_hs_list <- purrr::map(my_h5_files, ~hs_read(.x, bands = 1:nbands))
  for(i in 2:length(my_h5_files)){
    my_h5 <- merge(my_h5, my_hs_list[[i]])
  }
  names(my_h5) <- hs_wavelength(my_h5_files[1], bands = 1:nbands)
  
}

# extract spectra from pixels within the polygon
my_pt_spatvec <- my_aq_prj_sp %>% as('SpatVector')

my_hs_terra <- my_h5 %>% terra::rast()

my_spectra <- terra::extract(my_hs_terra, my_pt_spatvec, cells = TRUE)
my_spectra_df <- my_spectra %>% 
  as_tibble() %>%
  mutate(x = terra::xFromCell(my_hs_terra, cell),
         y = terra::yFromCell(my_hs_terra, cell)) %>%
  group_by(cell) %>%
  pivot_longer(cols = starts_with('band'), 
               names_to = 'band', 
               values_to = 'reflectance') %>% 
  tidyr::separate(band, into = c('index', 'wavelength'), sep = "_") %>%
  dplyr::mutate(wavelength = parse_number(wavelength)) 

unique(my_spectra_df$cell) %>% length()
mysamp <- sample(unique(my_spectra_df$cell), size = 100)

my_spectra_df %>%
  dplyr::filter(cell %in% mysamp) %>%
  # dplyr::filter(cell == 442414) %>%
  ggplot(aes(wavelength, reflectance)) +
  geom_line(aes(group = cell)) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme(legend.position = 'none') +
  ggtitle(glue('{my_aq_site}-{my_aop_yr}'))

# ggsave(glue('figs/buff5m-spectra/{my_aq_site}-{my_aop_yr}-5mbuff-{str_replace_all(my_loc_type, "[:punct:]", "")}.png'))


my_spectra_df %>% 
  vroom::vroom_write(glue('{spectra_dir}/{my_aq_site}_{my_aop_yr}.tsv'))

bands_to_plot <- seq(1, 100, length.out = 12)
# bands_to_plot <- 1:10
# my_h5 <- hs_read(my_h5_files, bands = bands_to_plot)
# names(my_h5) <- hs_wavelength(my_h5_files, bands = bands_to_plot)

# my_aq_prj_buff <- st_buffer(my_aq_prj, 1000)
my_h5_bands <- raster::subset(my_h5, subset = bands_to_plot)
# my_hs_crop <- raster::crop(my_h5_bands, my_aq_prj_buff)
my_hs_crop <- raster::crop(my_h5_bands, my_aq_prj)
hs_crop_stars <- my_hs_crop %>% st_as_stars()

gg <- ggplot() +
  geom_stars(data = hs_crop_stars) +
  coord_equal() +
  facet_wrap(~band) +
  geom_sf(data = my_aq_prj, fill = NA, col = 'white') +
  theme_void() + 
  theme(legend.title = element_blank()) +
  scale_fill_scico(palette = 'davos') +
  ggtitle(glue('{my_aq_site} {my_aop_yr}'))

# ggsave(glue('figs/site-maps/{my_aq_site}_{my_aop_yr}_firstbands.pdf'), plot = gg, width = 10, height = 8)
ggsave(glue('figs/site-maps/{my_aq_site}_{my_aop_yr}.pdf'), plot = gg, width = 10, height = 8)
# }

#
#

# STREAMS
