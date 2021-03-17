# extract spectra from pixels within
# a buffered distance of a point

library(neonhs)
library(viridis)
library(raster)
library(hdf5r)
library(sf)
library(tidyverse)
library(glue)
source('R/myfxns.R')

# set up
spectra_dir <- 'H:/DATA/spectra/'
shp_dir <- 'H:/DATA/spatial'

# inputs
my_aq_site <- 'PRLA'
my_aop_yr <- '2019'
my_loc_type <- 'buoy.c0'

get_aop_dates <- function(aq_siteids){
  aop_file <- 'results/sites_join_aop_dates.csv'
  aop_dates <- readr::read_csv(aop_file, col_types = 'ccccccccddD') %>%
    dplyr::filter(siteID %in% aq_siteids) %>%
    dplyr::select(siteID, aop_site_id, flightdate) %>%
    dplyr::arrange(flightdate) %>% distinct()
  return(aop_dates)
}

# AOP site name
my_aop_site <- get_aop_dates(my_aq_site) %>% 
  pull(aop_site_id) %>% unique()

# read in points
my_aq_sf <- glue('{shp_dir}/swchem_sites/swchem_sites.shp') %>%
  st_read() %>% dplyr::filter(siteID == my_aq_site) %>%
  mutate(location_type = substr(nmdLctn, 10, nchar(nmdLctn)))
my_loc_types <- my_aq_sf$location_type
my_domain <- my_aq_sf$domanID[1]
my_aq_pt_sf <- my_aq_sf %>% 
  dplyr::filter(location_type == my_loc_type)

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

# only for one point
my_aq_prj <- my_aq_pt_sf %>% st_transform(crs = my_epsg)
# buffer point
my_aq_prj_sp <- as(my_aq_prj, "Spatial")
my_aq_prj_buff5 <- st_buffer(my_aq_prj, 5)
my_buff5_sp <- as(my_aq_prj_buff5, "Spatial")

bbox_pts_sf <- my_aq_prj_buff5 %>% 
  st_bbox() %>% st_as_sfc() %>% 
  st_as_sf() %>% st_cast(to = "POINT")

# find the images that span the entire buffered point
my_site_files_df <- my_site_files %>%
  as.data.frame() %>% 
  rename(fullname = 1) %>%
  mutate(filename = basename(fullname)) %>%
  mutate(easting = as.numeric(str_sub(filename, 19, 24)),
         easting2 = easting + 1000,
         northing = as.numeric(str_sub(filename, 26, 32)),
         northing2 = northing + 1000)

# find which file one point is in
my_xys <- sf::st_coordinates(bbox_pts_sf) %>% as.data.frame()
# my_id <- 1
# need my_xys, my_site_files_df
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

if(length(my_h5_files) == 1){
  my_h5 <- hs_read(my_h5_files, bands = 1:100)
}
if(length(my_h5_files) > 1){
  my_h5 <- hs_read(my_h5_files[1], bands = 1:100)
  my_hs_list <- purrr::map(my_h5_files, ~hs_read(.x, bands = 1:100))
  for(i in 2:length(my_h5_files)){
    my_h5 <- merge(my_h5, my_hs_list[[i]])
  }
}

my_buff5_sp <- as(my_aq_prj_buff5, "Spatial")
my_aq_prj_sp <- as(my_aq_prj, "Spatial")
plot(my_h5[[100]], col = cividis(100), axes = FALSE, box = FALSE,
     main = glue('{my_aq_site} {my_aop_yr} \n {names(my_h5)[3]}'))
plot(my_buff5_sp, add = TRUE, border = 'green', col = NA)
plot(my_aq_prj_sp, add = TRUE, col = 'red', pch = 3)

##
my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
hdf5r::h5close(my_h5)

my_hs <- hs_read(my_h5_file, bands = c(1, 50, 100, 400))
plot(my_hs[[3]], col = cividis(100), axes = FALSE, box = FALSE,
     main = glue('{my_aq_site} {my_aop_yr} \n {names(my_hs)[3]}'))

# buffer point

my_aq_prj_buff30 <- st_buffer(my_aq_prj, 30)
my_buff30_sp <- as(my_aq_prj_buff30, "Spatial")

my_aq_prj_buff100 <- st_buffer(my_aq_prj, 100)
my_buff100_sp <- as(my_aq_prj_buff100, "Spatial")

my_hs_crop5 <- raster::crop(my_hs, my_buff5_sp)
my_hs_crop30 <- raster::crop(my_hs, my_buff30_sp)
my_hs_crop100 <- raster::crop(my_hs, my_buff100_sp)

plot(my_hs_crop5[[3]], col = cividis(100), axes = FALSE, box = FALSE)

plot(my_hs[[3]], col = cividis(100), axes = FALSE, box = FALSE)
plot(my_hs_crop100[[3]], col = cividis(100), axes = FALSE, box = FALSE)
plot(my_buff30_sp, add = TRUE, border = 'green', col = NA)
plot(my_buff5_sp, add = TRUE, border = 'green', col = NA)
plot(my_aq_prj_sp, add = TRUE, col = 'red', pch = 3)

# extract spectra from pixels within the 5m buffered point
my_hs
library(terra)
my_pts_spatvec1 <- my_aq_prj_sp %>% as('SpatVector')
my_pts_spatvec <- my_buff5_sp %>% as('SpatVector')

raster::extract(my_hs, my_aq_prj_sp)

my_hs <- hs_read(my_h5_file, bands = 1:200) # or 1:426
# this is slow. faster with terra?
# didnt work with all bands...
my_hs_terra <- my_hs %>% terra::rast()
my_spectra <- terra::extract(my_hs_terra, my_pts_spatvec)
my_spectra1 <- terra::extract(my_hs_terra, my_pts_spatvec1)
my_spectra_df <- my_spectra %>% 
  as_tibble() %>%
  rowid_to_column() %>%
  group_by(rowid) %>%
  pivot_longer(cols = starts_with('band'), 
               names_to = 'band', 
               values_to = 'reflectance') %>% 
  tidyr::separate(band, into = c('index', 'wavelength'), sep = "_") %>%
  dplyr::mutate(wavelength = parse_number(wavelength))

my_spectra_df1 <- my_spectra1 %>% 
  as_tibble() %>%
  rowid_to_column() %>%
  group_by(rowid) %>%
  pivot_longer(cols = starts_with('band'), 
               names_to = 'band', 
               values_to = 'reflectance') %>% 
  tidyr::separate(band, into = c('index', 'wavelength'), sep = "_") %>%
  dplyr::mutate(wavelength = parse_number(wavelength))

my_spectra_df %>% head()

my_spectra_df %>%
  ggplot(aes(wavelength, reflectance)) +
  geom_line(aes(group = rowid)) +
  geom_line(data = my_spectra_df1, col = 'green', lwd = 2) +
  # coord_cartesian(ylim = c(0, 0.1),
  #                 xlim = c(380, 1000)) +
  theme_bw() +
  theme(legend.position = 'none')

ggsave('prpo2019-5mbuff-buoyc0.png')

# make layers of cell IDs and or XYs
# extract all bands
# OR convert all pixel centroids to use hs_extract 

