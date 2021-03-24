# library(neonhs)
# library(viridis)
library(raster)
library(hdf5r)
library(sf)
library(tidyverse)
library(glue)
# library(terra)
# library(stars)
# library(scico)
source('R/myfxns.R')

# set up
shp_dir <- 'H:/DATA/spatial'

# inputs
my_aq_site <- 'BLWA'
my_aop_yr <- '2017'
my_loc_type <- 'buoy.c0'
# nbands <- 100

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
source('R/hs_mapinfo.R')

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

my_h5_file <- my_h5_files[1]



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

