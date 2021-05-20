# extract solar zenith angle from h5

library(hdf5r)
library(tidyverse)
library(glue)
library(neonhs)
library(sf)
library(terra)
library(raster)
library(dplyr)

# set up
source('R/myfxns.R')
setupdir <- '/Volumes/hondula/DATA' # mac
# setupdir <- 'H:' # DELL
spectra_dir <- glue('{setupdir}/spectra')
shp_dir <- glue('{setupdir}/spatial')
zenith_dir <- glue('{setupdir}/zenith-angles')

my_aq_site <- 'BARC'
my_aop_yr <- '2017'
my_loc_type <- 'buoy.c0'

# get domain and aop site
my_aop_site <- get_aop_dates(my_aq_site) %>% 
  pull(aop_site_id) %>% unique()
my_aq_sf <- glue('{shp_dir}/swchem_sites/swchem_sites.shp') %>%
  st_read() %>% dplyr::filter(siteID == my_aq_site) %>%
  mutate(location_type = substr(nmdLctn, 10, nchar(nmdLctn))) %>%
  dplyr::filter(location_type == my_loc_type)

my_loc_types <- my_aq_sf$location_type
my_domain <- my_aq_sf$domanID[1]

# get the site year directory
# my_site_dir <- glue('D:/{my_aop_yr}/FullSite/{my_domain}') %>%
#   fs::dir_ls(glob = glue("*{my_aop_site}*"), type = 'directory')
# my_site_files <- glue('{my_site_dir}/L3/Spectrometer/Reflectance') %>%
#   fs::dir_ls(glob = '*.h5')

my_site_files <- glue('{setupdir}/AOP/Reflectance/{my_aop_yr}/{my_aq_site}') %>%
  fs::dir_ls(glob = "*.h5")

# get the projection info from the first file
my_h5_file <- my_site_files[1]
my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
epsg_path <- glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/EPSG Code')
my_epsg <- my_h5[[epsg_path]]$read() %>% as.integer()

# only for one point
my_aq_prj <- my_aq_sf %>% 
  dplyr::filter(location_type %in% my_loc_type) %>%
  st_transform(crs = my_epsg)
my_xy <- sf::st_coordinates(my_aq_prj) %>% as.data.frame()
my_x <- my_xy[['X']]
my_y <- my_xy[['Y']]

# find the right h5 file - assume file names are correct UTMs

my_h5_file <- my_site_files %>%
  as.data.frame() %>% 
  rename(fullname = 1) %>%
  mutate(filename = basename(fullname)) %>%
  mutate(easting = as.numeric(str_sub(filename, 19, 24)),
         easting2 = easting + 1000,
         northing = as.numeric(str_sub(filename, 26, 32)),
         northing2 = northing + 1000) %>%
  mutate(x_in_range = {{ my_x }} >= easting & {{ my_x }} <= easting2,
         y_in_range = {{ my_y }} >= northing & {{ my_y }} <= northing2) %>%
  dplyr::filter(x_in_range, y_in_range) %>% 
  pull(fullname)

my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")

# try to get zenith angle dataset
my_h5$ls()
hdf5r::list.objects(my_h5) # everything? 59 things

best_date <- '164805'
zenith_path <- glue('{my_aop_site}/Reflectance/Metadata/Logs/{best_date}/Solar_Zenith_Angle')
my_h5[[zenith_path]]$read()
my_h5[['OSBS/Reflectance/Metadata/Spectral_Data/Wavelength']]$read()
my_h5[['OSBS/Reflectance/Metadata/Ancillary_Imagery/Data_Selection_Index']]$read()

zenith_path <- glue('{my_aop_site}/Reflectance/Metadata/Logs/{best_date}/Solar_Zenith_Angle')
zenith_angle <- my_h5[[zenith_path]]$read()
# 
# # zenith_path <- glue('{my_aop_site}/Reflectance/Metadata/to-sensor_zenith_angle')
# # zenith_raster <- my_h5[[zenith_path]]$read() %>% raster()
# 
# source('R/hs_mapinfo_zenith.R')
# 
# my_mapinfo <- hs_mapinfo_zenith(my_h5_file)
# extent(zenith_raster) <- raster::extent(my_mapinfo[["xmin"]],
#                                         my_mapinfo[["xmax"]],
#                                         my_mapinfo[["ymin"]],
#                                         my_mapinfo[["ymax"]])
# crs(zenith_raster) <- my_epsg
# 
# zenith_raster
# plot(zenith_raster)
# # plot(my_buff5_sp, add = TRUE)
# 
# # mask to buffered point
# # project and buffer point
# my_aq_prj <- my_aq_sf %>% st_transform(crs = my_epsg)
# my_aq_prj_buff5 <- st_buffer(my_aq_prj, 5)
# # convert to sp objects
# my_aq_prj_sp <- as(my_aq_prj, "Spatial")
# my_buff5_sp <- as(my_aq_prj_buff5, "Spatial")
# 
# my_pts_spatvec <- my_buff5_sp %>% as('SpatVector')
# my_pt_spatvec <- my_aq_prj_sp %>% as('SpatVector')
# my_zenith_terra <- zenith_raster %>% terra::rast()
# 
# # find which cell the exact coords of sampling
# my_cellid <- cells(my_zenith_terra, my_pt_spatvec)[2]
# 
# my_zeniths <- terra::extract(my_zenith_terra, my_pts_spatvec, cells = TRUE)
# my_zeniths_df <- my_zeniths %>% 
#   as_tibble() %>%
#   mutate(x = terra::xFromCell(my_zenith_terra, cell),
#          y = terra::yFromCell(my_zenith_terra, cell)) %>%
#   dplyr::mutate(my_cellid = cell == my_cellid) %>%
#   dplyr::rename(to_sensor_zenith_angle = layer) %>%
#   dplyr::mutate(aq_site = my_aq_site, aop_year = my_aop_yr)
# 
# zenith_path <- glue('{zenith_dir}/{my_aq_site}_{my_aop_yr}_{str_replace_all(my_loc_type, "[:punct:]", "")}-buff5m.csv')  
# my_zeniths_df %>% write_csv(zenith_path)

# weather_path <- glue('{my_aop_site}/Reflectance/Metadata/Ancillary_Imagery/Weather_Quality_Indicator')
# weather_raster <- my_h5[[weather_path]]$read()
# as.raster(weather_raster[3,,]) %>% str()
# str(weather_raster)
# weather_raster[1,,]
# ?as.raster.array
