## function to search KMLs for flightlines overlapping points

library(tidyverse)
library(fs)
library(glue)
library(httr)
library(tidyverse)
library(jsonlite)
library(sf)

# set up
l1_refl_id <- 'DP1.30006.001'
my_aq_site <- 'BARC'
my_aop_site <- 'OSBS'
my_domain <- 'D03'
my_aq_polygon <- 'BARC_AOSpts'
my_water_sf <- st_read(glue('/Volumes/hondula/DATA/AOP/site-polygons/{my_aq_polygon}.shp'))

get_flightline_ids('BLUE', 'BLUE_AOSpts')
get_flightline_ids('DELA', 'BLWA_AOSpts')
get_flightline_ids('LENO', 'TOMB_AOSpts')
get_flightline_ids('JERC', 'FLNT_AOSpts')
get_flightline_ids('LIRO', 'LIRO_AOSpts')
get_flightline_ids('OSBS', 'SUGG_AOSpts')
get_flightline_ids('OSBS', 'BARC_AOSpts')
get_flightline_ids('WOOD', 'PRPO_AOSpts')
get_flightline_ids('WOOD', 'PRLA_AOSpts')
get_flightline_ids('UNDE', 'CRAM_AOSpts')
get_flightline_ids('TOOL', 'TOOK_AOSpts')

get_flightline_ids <- function(my_aop_site, 
                               # my_aop_yr, 
                               my_aq_polygon,
                               kmls_base_dir = '/Volumes/hondula/DATA/AOP/ReflectanceL1',
                               keep_kmls = FALSE){
  
base_url <- 'http://data.neonscience.org/api/v0/'
l1_refl_id <- 'DP1.30006.001'
data_id <- l1_refl_id # l1 refl DP1.30006.001
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
  dplyr::select(siteid, month, url) %>% 
  mutate(aop_yr = str_sub(month, 1, 4))

  my_site_aop_df <- avail_df %>% 
    dplyr::filter(siteid == my_aop_site)

# my_site_aop_df$month
# my_site_aop_df$aop_yr

# filter to year of interest 
# my_site_aop_df <- my_site_aop_df %>% 
#   dplyr::filter(aop_yr %in% my_aop_yr)

my_site_urls <- my_site_aop_df %>% pull(url)

# get URL
# my_url <- my_site_urls[1]

# find the KML file that overlaps the points of interest

get_pattern_files <- function(my_url, myglob = 'kml'){
  data_files_req <- GET(my_url)
  data_files <- content(data_files_req, as = "text") %>%
    fromJSON(simplifyDataFrame = TRUE, flatten = TRUE)
  data_files_df <- data_files$data$files %>% 
    filter(str_detect(name, glue('{myglob}')))
  # future enhancement: check md5 sums for changes! 
  return_list <- NULL
  if(nrow(data_files_df) > 0){
    return_list <- list(files = data_files_df$name, urls = data_files_df$url)}
  return(return_list)
}

my_files_list <- my_site_urls %>% purrr::map(~get_pattern_files(.x, 'kml'))
my_files <- my_files_list %>% map(~as.data.frame(.x)) %>% bind_rows()
# my_files <- my_files_list[[1]]
kmls_dir <- glue('{kmls_base_dir}/{my_aop_yr}/{my_aop_site}/kmls')
fs::dir_create(kmls_dir)
# my_files <- my_files_list[[4]]

n_kmls <- length(my_files$urls)
message(glue('downloading {n_kmls} kmls from {my_aop_site} {my_aop_yr}'))
# download_month <- function(my_files){
my_kmls_local <- glue('{kmls_dir}/{my_files$files}')
# TODO: check for files already downloaded
purrr::walk2(.x = my_files$urls, .y = my_kmls_local, ~download.file(.x, .y))
# }
n_kmls_local <- length(fs::dir_ls(kmls_dir))
message(glue('downloaded {n_kmls_local} kmls from {my_aop_site} {my_aop_yr}'))


# NOW find which kmls overlap aos_pts locations
my_water_sf <- st_read(glue('/Volumes/hondula/DATA/AOP/site-polygons/{my_aq_polygon}.shp'))
# my_kml_file <- my_kmls_local[1]
test_in_kml <- function(my_kml_file, my_water_sf){
  kml1 <- st_read(my_kml_file) %>%
    st_zm() %>% filter(st_is(geometry, c("POLYGON", "MULTIPOLYGON"))) %>%
    st_make_valid()
  my_aq_prj <- my_water_sf %>% st_transform(st_crs(kml1)) # use EPSG for dell
  my_mat <- st_intersects(my_water_sf, kml1, sparse = FALSE)
  any_in_kml <- apply(my_mat, 1, any) %>% any()
  return(any_in_kml)
}

safe_test_in_kml <- possibly(test_in_kml, otherwise = FALSE)
kml_tests <- my_kmls_local %>% purrr::map_lgl(~safe_test_in_kml(.x, my_water_sf))
my_kmls <- my_kmls_local[which(kml_tests)]

# get weather info from kmls of interest
get_kml_colorcode <- function(kml_file){
  my_xml <- kml_file %>% xml2::read_xml() %>% xml2::as_list()
  my_colorcode <- unlist(my_xml$kml$Document$Folder$Placemark$Style$LineStyle$color)
  return(my_colorcode)
}

my_colorcodes <- my_kmls %>% purrr::map_chr(~get_kml_colorcode(.x))

# put data frame together
my_flightlines <- my_kmls %>% basename() %>% str_sub(1, 15)
kml_cols <- read_csv('kml-colors.csv')

my_flightlines_df <- data.frame(aop_site = my_aop_site,
           shp = my_aq_polygon,
           flightlines = my_flightlines,
           kml_color = my_colorcodes) %>%
  left_join(kml_cols)
if(!keep_kmls){fs::file_delete(my_kmls_local)}
write_csv(my_flightlines_df, glue('results/l1-flightlines/{my_aop_site}-{my_aq_polygon}-flightlines.csv'))
return(my_flightlines_df)
}



## Once you know which flightlines (from kmz)

# get_pattern_files <- function(my_url, myglob = 'reflectance.h5'){
#   data_files_req <- GET(my_url)
#   data_files <- content(data_files_req, as = "text") %>%
#     fromJSON(simplifyDataFrame = TRUE, flatten = TRUE)
#   data_files_df <- data_files$data$files %>% 
#     filter(str_detect(name, glue('{myglob}'))) 
#   # future enhancement: check md5 sums for changes! 
#   return_list <- NULL
#   if(nrow(data_files_df) > 0){
#     return_list <- list(files = data_files_df$name, urls = data_files_df$url)}
#   return(return_list)
# }
# 
# my_files_list <- my_site_urls %>% purrr::map(~get_pattern_files(.x, 'reflectance.h5'))
# 
# my_h5s_df <- my_files_list[[1]] %>% 
#   as.data.frame() %>% 
#   mutate(flightline = str_sub(files, 19, 33)) %>%
#   dplyr::filter(flightline %in% my_flightlines)
# # first flightline
# # my_h5_file_id <- my_files_list[[1]]$files %>% 
# #   purrr::map_lgl(~str_detect(.x, my_flightlines)) %>% which()
# # my_h5_file <- my_files_list[[4]]$files[my_h5_file_id]
# # my_h5_url <- my_files_list[[4]]$urls[my_h5_file_id]
# 
# h5_local_dir <- '/Volumes/hondula/DATA/AOP/ReflectanceL1'
# my_files_local <- fs::dir_ls(h5_local_dir, recurse = TRUE, regexp = my_h5_file)
# my_files_local <- glue('{kmls_dir}/{my_h5_file}')
# ## works up until here
# getOption('timeout')
# options(timeout=1200)
# # approx 5 GB
# library(curl)
# curl_download(url = my_h5_url, destfile = my_files_local)
# # download.file(url = my_h5_url, destfile = my_files_local)
# # purrr::walk2(.x = my_files$urls, .y = my_files_local, ~download.file(.x, .y))
# 
# my_h5_file_local <- my_files_local
# 
# neonhs::hs_epsg(my_h5_file_local)
# my_h5 <- hdf5r::H5File$new(my_h5_file_local, mode = "r")
# hdf5r::h5close(my_h5)
