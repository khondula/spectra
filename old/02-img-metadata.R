# find exact timing of image
# for a flight year

library(hdf5r)
library(tidyverse)
library(glue)
library(neonhs)

my_aq_site <- 'PRLA'
my_aop_yr <- '2019'
shp_dir <- '../../DATA/spatial'
my_loc_type <- 'buoy.c0'
  
whichh5 <- read_csv('results/which-h5.csv', col_types = 'cccccdc') %>% 
  dplyr::filter(year %in% my_aop_yr, siteID %in% my_aq_site)

my_h5_files <- whichh5 %>% 
  dplyr::filter(location_type == 'buoy.c0') %>%
  pull(h5_files)

# read in points
my_aq_sf <- glue('{shp_dir}/swchem_sites/swchem_sites.shp') %>%
  st_read() %>% dplyr::filter(siteID == my_aq_site) %>%
  mutate(location_type = substr(nmdLctn, 10, nchar(nmdLctn)))
my_loc_types <- my_aq_sf$location_type
my_domain <- my_aq_sf$domanID[1]
# AOP site name
my_aop_site <- get_aop_dates(my_aq_site) %>% 
  pull(aop_site_id) %>% unique()
my_site_dir <- glue('D:/{my_aop_yr}/FullSite/{my_domain}') %>%
  fs::dir_ls(glob = glue("*{my_aop_site}*"), type = 'directory')


# find metadata with image acquisition!
my_h5_files_full <- glue('{my_site_dir}/L3/Spectrometer/Reflectance/{my_h5_files}')
# which one has THE point of interest

my_prj <- hs_epsg(my_h5_files_full[1]) %>% as.integer()
my_aq_prj <- my_aq_sf %>% 
  dplyr::filter(location_type %in% my_loc_type) %>%
  st_transform(crs = my_prj)
pt_in_img <- function(path_to_file){
  my_extent <- hs_extent(path_to_file) %>%
    as("SpatialPolygons") %>% 
    st_as_sf() %>% 
    st_set_crs(my_prj)
  mat1 <- st_intersects(my_extent, my_aq_prj, sparse = FALSE)
  pts_in_extent <- which(apply(mat1, 1, any))
  pt_in_img <- my_extent[pts_in_extent,]
  return(nrow(pt_in_img)>0)
}

site_files_list <- my_h5_files_full %>% 
  purrr::map_lgl(~pt_in_img(.x))
my_file <- my_h5_files_full[which(site_files_list)]

# open a file for reading
my_h5 <- hdf5r::H5File$new(my_file, mode = "r")
class(my_h5)

hdf5r::list.objects(my_h5) # everything? 59 things
hdf5r::list.datasets(my_h5) # 47 datasets
hdf5r::list.groups(my_h5) # 12 groups

my_h5$ls()$name

my_h5[[glue('{my_aop_site}/Reflectance/Metadata/to-sensor_azimuth_angle')]]
my_h5[[glue('{my_aop_site}/Reflectance/Metadata/Ancillary_Imagery/Weather_Quality_Indicator')]]

my_ft <- my_h5[[glue('{my_aop_site}/Reflectance/Metadata/Flight_Trajectory')]]
my_ft
  class(my_ft)
my_h5$ls()


hdf5r::h5close(my_h5)
# look at weather quality indicator
# get to sensor azimuth angle
# get to sensor zenith angle

my_aq_sp <- as(my_aq_prj, "Spatial")

my_vals <- hs_extract_pts(my_file, pts = my_aq_sp, bands = 1:426)

my_vals_df <- my_vals %>% st_as_sf() %>% 
  st_drop_geometry() %>% as_tibble() %>%
  pivot_longer(cols = starts_with('band'), 
               names_to = 'band', 
               values_to = 'reflectance') %>% 
  tidyr::separate(band, into = c('index', 'wavelength'), sep = "_") %>%
  dplyr::mutate(wavelength = parse_number(wavelength))

my_vals_df %>% write_csv(glue('{spectra_dir}/{my_aq_site}_{my_aop_yr}.csv'))

