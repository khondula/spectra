l1_refl_id <- 'DP1.30006.001'
library(tidyverse)
library(fs)
library(glue)
library(httr)
library(tidyverse)
library(jsonlite)
my_site <- 'OSBS'

base_url <- 'http://data.neonscience.org/api/v0/'
data_id <- 'DP1.30006.001' # l1 refl
req_avail <- GET(glue('{base_url}/products/{data_id}'))
avail_resp <- content(req_avail, as = 'text') %>% 
  fromJSON(simplifyDataFrame = TRUE, flatten = TRUE)

# List of products by site code with month
data_urls_list <- avail_resp$data$siteCodes$availableDataUrls
# make table of urls with site and months
avail_df <- data_urls_list %>%
  unlist() %>% as.data.frame() %>%
  dplyr::rename(url = 1) %>%
  mutate(siteid = str_sub(url, 56, 59)) %>%
  mutate(month = str_sub(url, 61, 67)) %>%
  dplyr::select(siteid, month, url)

my_site_to_get <- avail_df %>% 
  dplyr::filter(siteid == my_site) %>%
  dplyr::filter(as.numeric(substr(month, 1, 4)) > 2011)

my_site_urls <- my_site_to_get %>% pull(url)

# filter to just the waq_instantaneous basic files
my_url <- my_site_urls[4]

# find the KML file that overlaps the points of interest
# check out the RGB reflectance - GLINT? 
# 
myglob = 'kml'
get_pattern_files <- function(my_url, myglob = 'kml'){
  data_files_req <- GET(my_url)
  data_files <- content(data_files_req, as = "text") %>%
    fromJSON(simplifyDataFrame = TRUE, flatten = TRUE)
  data_files_df <- data_files$data$files %>% 
    filter(str_detect(name, glue('{myglob}')))
  # filter(str_detect(name, "(kml).*(kml)"))
  # future enhancement: check md5 sums for changes! 
  return_list <- NULL
  if(nrow(data_files_df) > 0){
    return_list <- list(files = data_files_df$name, urls = data_files_df$url)}
  return(return_list)
}

my_files_list <- my_site_urls %>% purrr::map(~get_pattern_files(.x, 'kml'))

new_files <- my_files_list %>% map_lgl(~!is.null(.x))
any_new <- any(new_files)
if(!any_new){message(glue('No new {myglob} data from {mysite}'))}
months_newfiles <- my_site_to_get[['month']][which(new_files)]

my_files_list[4]

kmls_dir <- '/Volumes/hondula/DATA/AOP/test-l1'
my_files <- my_files_list[[4]]

download_month <- function(my_files){
  my_files_local <- glue('{kmls_dir}/{my_files$files}')
  purrr::walk2(.x = my_files$urls, .y = my_files_local, ~download.file(.x, .y))
}
library(sf)

# find which kmls overlap aos_pts locations
my_aq_polygon <- 'BARC_AOSpts'
my_water_sf <- st_read(glue('/Volumes/hondula/DATA/AOP/site-polygons/{my_aq_polygon}.shp'))
my_kml_file <- '/Volumes/hondula/DATA/AOP/test-l1/20170927_170305_hsi_kml_0000.kml'
my_kml_file <- my_files_local[43]
test_in_kml <- function(my_kml_file, my_water_sf){
  kml1 <- st_read(my_kml_file) %>%
    st_zm() %>%
    filter(st_is(geometry, c("POLYGON", "MULTIPOLYGON"))) %>%
    st_make_valid()
  my_aq_prj <- my_water_sf %>% st_transform(st_crs(kml1)) # use EPSG for dell
  my_mat <- st_intersects(my_water_sf, kml1, sparse = FALSE)
  any_in_kml <- apply(my_mat, 1, any) %>% any()
  return(any_in_kml)
}

safe_test_in_kml <- possibly(test_in_kml, otherwise = FALSE)
safe_test_in_kml(my_files_local[1], my_water_sf)
safe_test_in_kml(my_kml_file, my_water_sf)

kml_tests <- my_files_local %>% purrr::map_lgl(~safe_test_in_kml(.x, my_water_sf))
# kml_tests <- my_files_local %>% purrr::map_lgl(~safe_test_in_kml(.x, my_water_sf))
# kml_tests %>% purrr::map(~.x$result) %>% flatten_lgl()
kml_tests
my_kmls <- my_files_local[which(kml_tests)]

my_flightlines <- my_kmls %>% basename() %>% str_sub(1, 15)


ggplot() +
  geom_sf(data = kml1) +
  geom_sf(data = my_water_sf)

## Once you know which flightlines (from kmz)

get_pattern_files <- function(my_url, myglob = 'reflectance.h5'){
  data_files_req <- GET(my_url)
  data_files <- content(data_files_req, as = "text") %>%
    fromJSON(simplifyDataFrame = TRUE, flatten = TRUE)
  data_files_df <- data_files$data$files %>% 
    filter(str_detect(name, glue('{myglob}'))) 
  # filter(str_detect(name, "(kml).*(kml)"))
  # future enhancement: check md5 sums for changes! 
  return_list <- NULL
  if(nrow(data_files_df) > 0){
    return_list <- list(files = data_files_df$name, urls = data_files_df$url)}
  return(return_list)
}

my_files_list <- my_site_urls %>% purrr::map(~get_pattern_files(.x, 'reflectance.h5'))
my_files_list[[4]] %>% str()
# first flightline
my_h5_file_id <- my_files_list[[4]]$files %>% purrr::map_lgl(~str_detect(.x, my_flightlines[2])) %>% which()
my_h5_file <- my_files_list[[4]]$files[my_h5_file_id]
my_h5_url <- my_files_list[[4]]$urls[my_h5_file_id]

h5_local_dir <- '/Volumes/hondula/DATA/AOP/ReflectanceL1'
my_files_local <- fs::dir_ls(h5_local_dir, recurse = TRUE, regexp = my_h5_file)
my_files_local <- glue('{kmls_dir}/{my_h5_file}')
## works up until here
getOption('timeout')
options(timeout=1200)
# approx 5 GB
library(curl)
curl_download(url = my_h5_url, destfile = my_files_local)
# download.file(url = my_h5_url, destfile = my_files_local)
# purrr::walk2(.x = my_files$urls, .y = my_files_local, ~download.file(.x, .y))

my_h5_file_local <- my_files_local

neonhs::hs_epsg(my_h5_file_local)
my_h5 <- hdf5r::H5File$new(my_h5_file_local, mode = "r")
hdf5r::h5close(my_h5)
