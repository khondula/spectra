# extract solar zenith angle from h5

library(hdf5r)
library(tidyverse)
library(glue)
library(neonhs)

# set up
source('R/myfxns.R')
spectra_dir <- 'H:/DATA/spectra/'
shp_dir <- '../../DATA/spatial'

my_aq_site <- 'PRLA'
my_aop_yr <- '2019'
my_loc_type <- 'buoy.c0'

# get domain and aop site
my_aop_site <- get_aop_dates(my_aq_site) %>% 
  pull(aop_site_id) %>% unique()
my_aq_sf <- glue('{shp_dir}/swchem_sites/swchem_sites.shp') %>%
  st_read() %>% dplyr::filter(siteID == my_aq_site) %>%
  mutate(location_type = substr(nmdLctn, 10, nchar(nmdLctn)))
my_loc_types <- my_aq_sf$location_type
my_domain <- my_aq_sf$domanID[1]

# get the site year directory
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
my_aq_prj <- my_aq_sf %>% 
  dplyr::filter(location_type %in% my_loc_type) %>%
  st_transform(crs = my_epsg)
my_xy <- sf::st_coordinates(my_aq_prj) %>% as.data.frame()
my_x <- my_xy[['X']]
my_y <- my_xy[['Y']]

# find the right h5 file - assume file names are correct UTMs
# which h5 file has point..
# for each file name get yy range from northing + 1km
# for each filename get xx range from easting + 1km
my_site_files %>%
  basename() %>%
  as.data.frame() %>% 
  rename(filename = 1) %>%
  mutate(easting = as.numeric(str_sub(filename, 19, 24)),
         easting2 = easting + 1000,
         northing = as.numeric(str_sub(filename, 26, 32)),
         northing2 = northing + 1000) %>%
  mutate(x_in_range = {{ my_x }} >= easting & {{ my_x }} <= easting2,
         y_in_range = {{ my_y }} >= northing & {{ my_y }} <= northing2) %>%
  dplyr::filter(x_in_range, y_in_range) %>% 
  pull(filename)

# which filename has range that contains pt coords

