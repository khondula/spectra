# install.packages('BiocManager')
# BiocManager::install('rhdf5')
library(neonhs)
# library(viridis)
library(raster)
library(hdf5r)
library(sf)
library(tidyverse)
library(glue)
library(lubridate)
library(rhdf5)
library(hms)
# library(terra)
# library(stars)
# library(scico)
# source('R/myfxns.R')
source('R/hs_mapinfo.R')

# set up
shp_dir <- 'H:/DATA/spatial'

# inputs
my_aq_site <- 'ARIK'
my_aop_yr <- '2020'
my_loc_type <- 'S2'



save_radiance_times(my_aq_site = 'ARIK', '2020', 'S2')



save_radiance_times <- function(my_aq_site, my_aop_yr, my_loc_type){
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

# folder with AOP data - RADIANCE
my_site_dir <- glue('D:/{my_aop_yr}/FullSite/{my_domain}') %>%
  fs::dir_ls(glob = glue("*{my_aop_site}*"), type = 'directory')
my_site_files <- glue('{my_site_dir}/L1/Spectrometer/RadianceH5') %>%
  fs::dir_ls(glob = '*.h5', recurse = TRUE)

# get the projection info from the first file
my_h5_file <- my_site_files[1]
my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
epsg_path <- glue('{my_aop_site}/Radiance/Metadata/Coordinate_System/EPSG Code')
my_epsg <- my_h5[[epsg_path]]$read() %>% as.integer()
my_h5$close_all()

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
my_xys <- sf::st_coordinates(my_aq_prj) %>% as.data.frame()

# find the images that span the entire buffered point
my_site_files_df <- my_site_files %>%
  as.data.frame() %>% 
  rename(fullname = 1) %>%
  mutate(xmin = purrr::map_dbl(fullname, ~hs_mapinfo(.x, 'Radiance')[['xmin']]),
         xmax = purrr::map_dbl(fullname, ~hs_mapinfo(.x, 'Radiance')[['xmax']]),
         ymin = purrr::map_dbl(fullname, ~hs_mapinfo(.x, 'Radiance')[['ymin']]),
         ymax = purrr::map_dbl(fullname, ~hs_mapinfo(.x, 'Radiance')[['ymax']]))


find_which_file <- function(my_id){
  my_x <- my_xys[['X']][my_id]
  my_y <- my_xys[['Y']][my_id]
  my_h5_file <- my_site_files_df %>%
    mutate(x_in_range = {{ my_x }} >= xmin & {{ my_x }} <= xmax,
           y_in_range = {{ my_y }} >= ymin & {{ my_y }} <= ymax) %>%
    dplyr::filter(x_in_range, y_in_range) %>% 
    pull(fullname) %>% unique()
  return(my_h5_file)
}

my_h5_files <- 1:nrow(my_xys) %>% 
  purrr::map(~find_which_file(.x)) %>% 
  unlist() %>% unique()

### extract time from OBS raster
### for each file overlapping point

# my_h5_file <- my_h5_files[1]

extract_gpstime <- function(my_h5_file, my_aq_prj_sp){
  file_h5 <- hdf5r::H5File$new(my_h5_file, mode = 'r')
  site <- file_h5$ls()$name
  obs_path <- glue::glue('{site}/Radiance/Metadata/Ancillary_Rasters/OBS_Data')
  my_info <- h5readAttributes(my_h5_file, obs_path)
  my_nrows <- my_info$Dimensions[2]
  my_ncols <- my_info$Dimensions[1]
  # gps time is the 10th band
  obs1 <- h5read(my_h5_file, obs_path, index = list(10, 1:my_nrows, 1:my_ncols))
  obs1_matrix <- matrix(obs1[1,,], nrow = my_nrows, ncol = my_ncols)
  
  # save this for getting all obs data layers
  # my_stack <- purrr::map(1:10, ~matrix(obs1[.x,,], nrow = 8264, ncol = 1017)) %>%
  #   purrr::map(~raster(.x)) %>%
  #   stack()
  # ids <- values(my_stack[[1]]) == my_info$Data_Ignore_Value
  # my_stack[ids] <- NA
  
  # make sure to transpose!! 
  my_r <- raster(t(obs1_matrix))
  # replace No data pixels
  ids <- raster::values(my_r[[1]]) == my_info$Data_Ignore_Value
  my_r[ids] <- NA
  
  # set extent and CRS
  extent(my_r) <- raster::extent(hs_mapinfo(my_h5_file)[["xmin"]],
                                 hs_mapinfo(my_h5_file)[["xmax"]],
                                 hs_mapinfo(my_h5_file)[["ymin"]],
                                 hs_mapinfo(my_h5_file)[["ymax"]])
  crs(my_r) <- my_epsg
  
  plot(my_r)
  plot(my_aq_prj_sp, add = TRUE, pch = 4, col = 'blue', cex = 1)
  
  # extract spectra from pixels within the 5m buffered point
  my_pt_spatvec <- my_aq_prj_sp %>% as('SpatVector')
  my_hs_terra <- my_r %>% terra::rast()
  
  my_gpstime <- terra::extract(my_hs_terra, my_pt_spatvec, cells = TRUE)
  gsp_df1 <- my_gpstime %>%
    rename(gps_time = layer) %>% 
    mutate(x = terra::xFromCell(my_hs_terra, cell),
           y = terra::yFromCell(my_hs_terra, cell)) %>%
    mutate(filename = basename(my_h5_file)) %>%
    mutate(date = str_sub(filename, 19, 26)) %>%
    mutate(datestring = glue('{date} {hms(hours = gps_time)}'))
  gsp_df <- gsp_df1 %>%
    mutate(datetime = lubridate::as_datetime(datestring))
  

  return(gsp_df)
}


times_df <- my_h5_files %>% purrr::map_dfr(~extract_gpstime(.x, my_aq_prj_sp))

rad_times_dir <- 'H:/DATA/radiance-time'
times_df %>% write_csv(glue('{rad_times_dir}/{my_aq_site}-{my_aop_yr}_radianceTime.csv'))
}

# my_h5_files
# my_h5_file <- my_h5_files[1]
# 
# extract_gpstime(my_h5_files[1], my_aq_prj_sp)
