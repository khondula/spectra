# extract spectra from pixels within a polygon
# devtools::install_github('khondula/neonhs', force = TRUE)
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

# fs::dir_ls('H:/DATA/AOP/site-polygons')

# my_water_sf <- st_read('H:/DATA/AOP/site-polygons/PRLA_2016.shp')

nbands <- 426
my_aq_site <- 'BARC'
my_aop_yr <- '2017'
my_aop_site <- 'OSBS'
my_domain <- 'D03'
my_aq_polygon <- 'BARC_AOSpts'

# maybe working
# save_spectra('PRLA', '2017', 'WOOD', 'PRLA_2016', 'D09', 100)
# save_spectra('PRLA', '2019', 'WOOD', 'PRLA_2016', 'D09', 100)
# save_spectra('PRLA', '2020', 'WOOD', 'PRLA_2016', 'D09', 100)

# save_spectra('PRPO', '2017', 'WOOD', 'PRPO_2016', 'D09', 100)
# save_spectra('PRPO', '2019', 'WOOD', 'PRPO_2016', 'D09', 100)
# save_spectra('PRPO', '2020', 'WOOD', 'PRPO_2016', 'D09', 100)

#
save_spectra('CRAM', '2019', 'UNDE', 'CRAM_2015', 'D05', 100)
save_spectra('TOOK', '2018', 'TOOL', 'TOOK_2016', 'D18', 100)
save_spectra('FLNT', '2017', 'JERC', 'FLNT_2017', 'D03', 100)
save_spectra('FLNT', '2018', 'JERC', 'FLNT_2017', 'D03', 100)
save_spectra('FLNT', '2019', 'JERC', 'FLNT_2017', 'D03', 100)
save_spectra('BLWA', '2019', 'DELA', 'BLWA_2019', 'D08', 100)
save_spectra('TOMB', '2018', 'LENO', 'TOMB_2017_100m', 'D08', 100)


# save_spectra('BARC', '2018', 'OSBS', 'BARC_2018', 'D03', 100)
# save_spectra('BARC', '2019', 'OSBS', 'BARC_2018', 'D03', 100)

# save_spectra('SUGG', '2018', 'OSBS', 'SUGG_2018', 'D03', 100)
# save_spectra('SUGG', '2019', 'OSBS', 'SUGG_2018', 'D03', 100)

save_spectra <- function(my_aq_site, my_aop_yr, my_aop_site, 
                         my_aq_polygon,
                         my_domain, nbands){

  
  my_water_sf <- st_read(glue('H:/DATA/AOP/site-polygons/{my_aq_polygon}.shp'))
  
  spectra_dir <- 'H:/DATA/L1-refl-spectra/'
  
# folder with AOP data
  my_site_dir <- glue('D:/{my_aop_yr}/FullSite/{my_domain}') %>%
    fs::dir_ls(glob = glue("*{my_aop_site}*"), type = 'directory')
  my_site_files <- glue('{my_site_dir}/L1/Spectrometer/ReflectanceH5') %>%
    fs::dir_ls(glob = '*.h5', recurse = TRUE)

  # get the projection info from the first file
  my_h5_file <- my_site_files[1]
  my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
  epsg_path <- glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/EPSG Code')
  my_epsg <- my_h5[[epsg_path]]$read() %>% as.integer()
  hdf5r::h5close(my_h5)
  
  # wls_path <- glue('{my_aop_site}/Reflectance/Metadata/Spectral_Data/Wavelength')
  # my_wls <- my_h5[[wls_path]]$read()
  
  # project polygon to AOP data
  my_aq_prj <- my_water_sf %>% st_transform(crs = my_epsg)
  my_aq_prj_buff5 <- st_buffer(my_aq_prj, 5)
  # get corners of extent
  bbox_pts_sf <- my_aq_prj_buff5 %>% 
    st_bbox() %>% st_as_sfc() %>% 
    st_as_sf() %>% st_cast(to = "POINT")
  my_xys <- sf::st_coordinates(bbox_pts_sf) %>% 
    as.data.frame() %>% distinct()

# find the images that span the bbox of buffered points
# need to use map info not filenames
  my_site_files_df <- my_site_files %>%
    as.data.frame() %>% 
    rename(fullname = 1) %>%
    mutate(filename = basename(fullname))

  extent_list <- my_site_files %>% purrr::map(~hs_extent(.x))
  my_site_files_df <- my_site_files_df %>% 
    dplyr::mutate(xmin = purrr::map_dbl(extent_list, ~.x@xmin),
             xmax = purrr::map_dbl(extent_list, ~.x@xmax),
             ymin = purrr::map_dbl(extent_list, ~.x@ymin),
             ymax = purrr::map_dbl(extent_list, ~.x@ymax))
             
  
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

message(glue('{length(my_h5_files)} files'))
print(my_h5_files)

my_h5_files[3]

my_h5 <- hdf5r::H5File$new(my_h5_files[3], mode = "r")
solar_zenith_path <- glue('{my_aop_site}/Reflectance/Metadata/Logs/Solar_Zenith_Angle')
my_solar_zenith <- my_h5[[solar_zenith_path]]$read() 

epsg_path <- glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/EPSG Code')
my_epsg <- my_h5[[epsg_path]]$read() %>% as.integer()

# sensor zenith path
sensor_zenith_path <- glue('{my_aop_site}/Reflectance/Metadata/to-sensor_Zenith_Angle')
my_sensor_zenith <- my_h5[[sensor_zenith_path]]$read() %>% t() %>% terra::rast()
# replace No data pixels
my_sensor_zenith[my_sensor_zenith==-9999] <- NA


# read in and merge h5 files using hs_read
if(length(my_h5_files) == 1){
  my_h5 <- neonhs::hs_read(my_h5_files, bands = 1:nbands, crop = my_aq_prj)
  names(my_h5) <- hs_wavelength(my_h5_files, bands = 1:nbands)
}

if(length(my_h5_files) > 1){
  my_h5 <- neonhs::hs_read(my_h5_files[3], bands = 1:nbands)
  names(my_h5) <- hs_wavelength(my_h5_files[1], bands = 1:nbands)
  my_hs_list <- purrr::map(my_h5_files, ~hs_read(.x, bands = 1:nbands))
  for(i in 2:length(my_h5_files)){
    my_h5 <- merge(my_h5, my_hs_list[[i]])
  }
  names(my_h5) <- hs_wavelength(my_h5_files[1], bands = 1:nbands)
  
}


# extract spectra from pixels within the polygon
my_pt_spatvec <- my_aq_prj %>% as('SpatVector')
my_pts_spatvec <- my_aq_prj_buff5 %>% as('SpatVector')

my_hs_terra <- my_h5
my_cellid <- cells(my_hs_terra, my_pt_spatvec)[,2]
cellid_df <- data.frame(loctype = my_aq_prj$loctype,cellid = my_cellid)

terra::ext(my_sensor_zenith) <- terra::ext(my_hs_terra)
terra::crs(my_sensor_zenith) <- terra::crs(my_hs_terra)
names(my_sensor_zenith) <- 'sensorZenith'
# plot(my_sensor_zenith)
# plot(my_aq_prj$geometry, add = TRUE, pch = 4, col = 'blue', cex = 1)
my_hs_terra <- c(my_hs_terra, my_sensor_zenith)

message('extracting pixel values now')

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
  dplyr::mutate(wavelength = parse_number(wavelength))

my_spectra_df <- my_spectra_df %>% 
  left_join(cellid_df, by = c("cell" = "cellid")) %>%
  mutate(solar_zenith = my_solar_zenith)

ncells <- unique(my_spectra_df$cell) %>% length()
message(glue('{ncells} cells in your polygons'))

# mysamp <- sample(unique(my_spectra_df$cell), size = 100)
# 
# my_spectra_df %>%
#   dplyr::filter(cell %in% mysamp) %>%
#   # dplyr::filter(cell == 442414) %>%
#   ggplot(aes(wavelength, reflectance)) +
#   geom_line(aes(group = cell)) +
#   theme_bw() +
#   coord_cartesian(ylim = c(0, 0.2)) +
#   theme(legend.position = 'none') +
#   ggtitle(glue('{my_aq_site}-{my_aop_yr}'))

# ggsave(glue('figs/buff5m-spectra/{my_aq_site}-{my_aop_yr}-5mbuff-{str_replace_all(my_loc_type, "[:punct:]", "")}.png'))


my_spectra_df %>% 
  vroom::vroom_write(glue('{spectra_dir}/{my_aq_site}_{my_aop_yr}.tsv'))

# bands_to_plot <- seq(1, 100, length.out = 12)
# bands_to_plot <- 65:100
# my_h5 <- hs_read(my_h5_files, bands = bands_to_plot)
# names(my_h5) <- hs_wavelength(my_h5_files, bands = bands_to_plot)

# my_aq_prj_buff <- st_buffer(my_aq_prj, 50)
# my_hs_crop <- terra::crop(my_hs_terra, my_aq_prj_buff)
# hs_crop_stars <- my_hs_crop %>% st_as_stars()

# test_raster <- hs_crop_stars %>% as("Raster")
# library(leaflet)
# library(mapview)
# mapview(test_raster)

# gg <- ggplot() +
#   geom_stars(data = hs_crop_stars) +
#   coord_equal() +
#   facet_wrap(~band) +
#   geom_sf(data = my_aq_prj, fill = NA, col = 'white') +
#   theme_void() + 
#   theme(legend.title = element_blank()) +
#   scale_fill_scico(palette = 'romaO') +
#   # scale_fill_scico(palette = 'davos') +
#   ggtitle(glue('{my_aq_site} {my_aop_yr}'))
# gg
# ggsave(glue('figs/site-maps/{my_aq_site}_{my_aop_yr}-glint.pdf'), plot = gg, width = 10, height = 8)

# ggsave(glue('figs/site-maps/{my_aq_site}_{my_aop_yr}_firstbands.pdf'), plot = gg, width = 10, height = 8)
# ggsave(glue('figs/site-maps/{my_aq_site}_{my_aop_yr}.pdf'), plot = gg, width = 10, height = 8)

}

#
#

# STREAMS
