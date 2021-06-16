# extract spectra from pixels within a polygon
# devtools::install_github('khondula/neonhs', force = TRUE)

# library(neonhs)
# library(viridis)
# library(raster)
library(hdf5r)
library(sf)
# library(tidyverse)
library(glue)
# library(terra)
# library(stars)
# library(scico)
# library(vroom)
# source('R/myfxns.R')

# set up
# spectra_dir <- 'H:/DATA/spectra_site/'

# fs::dir_ls('H:/DATA/AOP/site-polygons')

# my_water_sf <- st_read('H:/DATA/AOP/site-polygons/PRLA_2016.shp')

data_dir <- '/Volumes/hondula/DATA'
my_aq_site <- 'OSBS'
my_aop_yr <- '2014'
my_aop_site <- 'OSBS'
my_aq_polygon <- 'SUGG_AOSpts'
my_domain <- 'D03'
flightline <- '20140507_140657'

aq_sites <- read_csv('../neon-sites/results/swchem_sites_df.csv') %>%
  dplyr::select(domainID, siteID) %>% distinct()

flightlines_df <- read_csv('results/l1-flightlines.csv') %>%
  mutate(aq_site = str_sub(shp, 1, 4), 
         aop_yr = str_sub(flightlines, 1, 4)) %>%
  left_join(aq_sites, by = c('aq_site' = 'siteID'))

# flightlines_df %>% arrange(shp, flightlines) %>% write_csv('results/l1-flightlines-wDomain.csv')
  
1:83 %>% 
  purrr::walk(~save_spectra(flightlines_df$aq_site[.x], 
                            flightlines_df$aop_yr[.x], 
                            flightlines_df$aop_site[.x], 
                            flightlines_df$shp[.x], 
                            flightlines_df$domainID[.x], 
                            flightlines_df$flightlines[.x]))

# save_spectra(flightlines_df$aq_site[1], 
#              flightlines_df$aop_yr[1], 
#              flightlines_df$aop_site[1], 
#              flightlines_df$shp[1], 
#              flightlines_df$domainID[1], 
#              flightlines_df$flightlines[1])

save_spectra('SUGG', '2014', 'OSBS', 'SUGG_AOSpts', 'D03', '20140507_154756')


save_spectra('CRAM', '2019', 'UNDE', 'CRAM_AOSpts', 'D05', '20190606_192344') # no cells?
save_spectra('PRPO', '2020', 'WOOD', 'PRPO_AOSpts', 'D09', '20200624_161319') # no cells
save_spectra('PRPO', '2020', 'WOOD', 'PRPO_AOSpts', 'D09', '20200624_162051') # no cells
save_spectra('PRPO', '2020', 'WOOD', 'PRPO_AOSpts', 'D09', '20200624_162804') # no cells

my_aq_site <- 'SUGG'
my_aop_yr <- '2014'
my_aop_site <- 'OSBS'
my_aq_polygon <- 'SUGG_AOSpts'
my_domain <- 'D03'
flightline <- '20140507_154756'

save_spectra <- function(my_aq_site, my_aop_yr, my_aop_site, 
                         my_aq_polygon, my_domain, flightline){
  
  data_dir <- '/Volumes/hondula/DATA'
  polygon_file <- glue::glue('{data_dir}/AOP/site-polygons/{my_aq_polygon}.shp')
  my_water_sf <- sf::st_read(polygon_file)
  spectra_out_dir <- glue::glue('{data_dir}/L1-refl-spectra')
  
  # check for duplicate named points in shp file... add new identifier?
  
  # get h5 file for flightline
  my_site_dir <- glue::glue('{data_dir}/AOP/ReflectanceL1/{my_aop_yr}/{my_aop_site}')
  my_h5_file <- fs::dir_ls(my_site_dir, glob = glue::glue('*{flightline}*.h5'))
  
  # file metadata
  my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
  epsg_path <- glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/EPSG Code')
  my_epsg <- my_h5[[epsg_path]]$read()
  
  solar_zenith_path <- glue('{my_aop_site}/Reflectance/Metadata/Logs/Solar_Zenith_Angle')
  my_solar_zenith <- my_h5[[solar_zenith_path]]$read() 
  
  wls_path <- glue('{my_aop_site}/Reflectance/Metadata/Spectral_Data/Wavelength')
  my_wls <- my_h5[[wls_path]]$read()
  
  # reflectance and metadata
  my_refl <- my_h5[[glue('{my_aop_site}/Reflectance/Reflectance_Data')]]
  scale_factor <- my_refl$attr_open('Scale_Factor')$read()
  na_value <- my_refl$attr_open('Data_Ignore_Value')$read()
  cloud_conditions <- my_refl$attr_open('Cloud conditions')$read()
  
  # extent
  
  map_info_path <- glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/Map_Info')
  map_info <- my_h5[[map_info_path]]$read() %>% strsplit(',') %>% unlist()
  my_dims <- my_refl$dims
  xy_resolution <- as.numeric(c(map_info[2], map_info[3]))
  xmin <- as.numeric(map_info[4])
  xmax <- xmin + my_dims[2] * xy_resolution[1]
  ymax <- as.numeric(map_info[5])
  ymin <- ymax - my_dims[3] * xy_resolution[2]
  my_extent <- terra::ext(xmin, xmax, ymin, ymax)
  
  # blank raster with map info
  my_epsg2 <- glue('EPSG:{my_epsg}')
  my_rast <- terra::rast(nrow = my_dims[2],
                         ncol = my_dims[3],
                         crs = my_epsg2,
                         extent = my_extent)
  
  # convert points to h5 projection
  my_aq_prj <- my_water_sf %>% sf::st_transform(my_epsg2) # use EPSG for dell?
  my_pts_spatvec <- my_aq_prj %>% as('SpatVector')
  
  # get spatial location of pts in raster dimensions
  my_cellids <- terra::cells(my_rast, my_pts_spatvec)[,2] # modify this for polygons
  my_cellrowcols <- terra::rowColFromCell(my_rast, my_cellids)
  my_cell_xys <- terra::xyFromCell(my_rast, my_cellids)

  cellinfo_df <- my_aq_prj %>% st_drop_geometry() %>%
    dplyr::mutate(cellid = my_cellids,
                  cellrow = my_cellrowcols[,1],
                  cellcol = my_cellrowcols[,2],
                  cellx = my_cell_xys[, 1],
                  celly = my_cell_xys[,2]) %>%
    dplyr::filter(!is.nan(cellid))
  
  if(nrow(cellinfo_df)<1){stop('cell ids are all NaN')}
  
  # should first check whether pixels are in NA region of raster!! 
  # otherwise getting lots of -9999s for no reason
  # get spectra 
  my_spectra_list <- purrr::map(1:nrow(cellinfo_df), ~my_refl[1:426, cellinfo_df[['cellrow']][.x], cellinfo_df[['cellcol']][.x]])
  names(my_spectra_list) <- cellinfo_df[['loctype']]
  
  spectra_df <- my_spectra_list %>% 
    as.data.frame() %>% 
    mutate(wl = my_wls) %>%
    mutate(band = glue('band_{str_pad(1:426, 3, "left", "0")}'))
  
  # divide by scale factor and take out no data value here??
  
  # spectra_df %>% 
  #   dplyr::select(-band) %>%
  #   tidyr::pivot_longer(cols = -wl, names_to = 'loctype') %>%
  #   ggplot(aes(x = wl, y = value/scale_factor)) +
  #   geom_line(aes(col = loctype)) +
  #   geom_point(aes(fill = loctype), alpha = 0.5, pch = 21) +
  #   xlim(c(NA, 750)) + ylim(c(0, 0.2)) +
  #   theme_bw() + 
  #   ggtitle(glue('{my_aq_site} {my_aop_yr}\n {flightline}, {cloud_conditions}'))
  # get pixel metadata - sensor angles

  # sensor zenith path
  sensor_zenith_path <- glue('{my_aop_site}/Reflectance/Metadata/to-sensor_Zenith_Angle')
  # same as reflectance
  my_sensor_zenith <- my_h5[[sensor_zenith_path]]
  my_sensor_zenith_list <- purrr::map(1:nrow(cellinfo_df), ~my_sensor_zenith[cellinfo_df[['cellrow']][.x], cellinfo_df[['cellcol']][.x]])
  names(my_sensor_zenith_list) <- cellinfo_df[['loctype']]
  
  cellinfo_df <- cellinfo_df %>%
    mutate(sensor_zenith = unlist(my_sensor_zenith_list),
           solar_zenith = my_solar_zenith,
           clouds = cloud_conditions,
           flightline = flightline,
           aop_site = my_aop_site,
           crs = my_epsg2)


  spectra_path <- glue('{data_dir}/L1-reflectance/spectra/{my_aq_site}_{flightline}.tsv')
  cellinfo_path <- glue('{data_dir}/L1-reflectance/meta/{my_aq_site}_{flightline}.tsv')
  spectra_df %>% vroom::vroom_write(spectra_path)
  cellinfo_df %>% vroom::vroom_write(cellinfo_path)
  
  message(glue('saved spectra for {nrow(cellinfo_df)} pts at {my_aq_site} from {flightline}'))
  hdf5r::h5close(my_h5)
  
}

#
#

