# extract spectra from pixels within a polygon
# devtools::install_github('khondula/neonhs', force = TRUE)

# library(neonhs)
# library(viridis)
# library(raster)
library(hdf5r)
library(sf)
library(tidyverse)
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

# data_dir <- '/Volumes/hondula/DATA'
data_dir <- 'H:/DATA/'

my_aq_site <- 'BARC'
my_aop_yr <- '2016'
my_aop_site <- 'OSBS'
my_aq_polygon <- 'BARC_AOSpts'
my_domain <- 'D03'


# save_radiance_gpstime('SUGG', '2014', 'OSBS', 'SUGG_AOSpts', 'D03')
# save_radiance_gpstime('BARC', '2014', 'OSBS', 'BARC_AOSpts', 'D03')
save_L1reflectance_spectra('PRLA', '2016', 'WOOD', 'PRLA_AOSpts', 'D09')
save_L1reflectance_spectra('PRPO', '2016', 'WOOD', 'PRPO_AOSpts', 'D09')

save_L1reflectance_spectra('BARC', '2016', 'OSBS', 'BARC_AOSpts', 'D03')
save_L1reflectance_spectra('BARC', '2019', 'OSBS', 'BARC_AOSpts', 'D03')
save_L1reflectance_spectra('SUGG', '2016', 'OSBS', 'SUGG_AOSpts', 'D03')
save_L1reflectance_spectra('SUGG', '2019', 'OSBS', 'SUGG_AOSpts', 'D03')

save_L1reflectance_spectra('FLNT', '2014', 'JERC', 'FLNT_AOSpts', 'D03')
save_L1reflectance_spectra('FLNT', '2016', 'JERC', 'FLNT_AOSpts', 'D03')

save_L1reflectance_spectra('BLWA', '2015', 'DELA', 'BLWA_AOSpts', 'D08')
save_L1reflectance_spectra('BLWA', '2016', 'DELA', 'BLWA_AOSpts', 'D08')

save_L1reflectance_spectra('TOMB', '2015', 'LENO', 'TOMB_AOSpts', 'D08')
save_L1reflectance_spectra('TOMB', '2016', 'LENO', 'TOMB_AOSpts', 'D08')


flightlines_df <- read_csv('results/l1-flightlines-wDomain.csv') %>%
  dplyr::select(aq_site, aop_yr, aop_site, shp, domainID) %>% dplyr::distinct()

1:nrow(flightlines_df) %>% 
  purrr::walk(~save_radiance_gpstime(flightlines_df$aq_site[.x], 
                                    flightlines_df$aop_yr[.x], 
                                    flightlines_df$aop_site[.x], 
                                    flightlines_df$shp[.x], 
                                    flightlines_df$domainID[.x]))

# note this assumes that all flightlines are in the same CRS for a siteyear
# need to redo BLUE to work with s2 loctype

save_L1reflectance_spectra <- function(my_aq_site, my_aop_yr, my_aop_site, 
                         my_aq_polygon, my_domain){
  
  data_dir <- 'H:/DATA/'
  polygon_file <- glue::glue('{data_dir}/AOP/site-polygons/{my_aq_polygon}.shp')
  my_water_sf <- sf::st_read(polygon_file)
  out_dir <- glue::glue('{data_dir}/L1-reflectance2')
  
  # first need to find which h5 files have points
  my_site_dir <- glue('D:/{my_aop_yr}/FullSite/{my_domain}') %>%
    fs::dir_ls(glob = glue('*{my_aop_site}*'), type = 'directory')
  my_site_files <- glue('{my_site_dir}/L1/Spectrometer/ReflectanceH5') %>%
    fs::dir_ls(glob = '*.h5', recurse = TRUE)
  
  # get the projection info from the first file
  my_h5_file <- my_site_files[1]
  
  my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
  epsg_path <- glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/EPSG Code')
  my_epsg <- my_h5[[epsg_path]]$read() %>% as.integer()
  my_epsg2 <- glue('EPSG:{my_epsg}')
  hdf5r::h5close(my_h5)
  
  # convert points to h5 projection
  my_aq_prj <- my_water_sf %>% sf::st_transform(my_epsg2)
  my_pts_spatvec <- my_aq_prj %>% as('SpatVector')

  # FOR each of the flightline files
  
  for(i in 1:length(my_site_files)){
    my_h5_file <- my_site_files[i]
    my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
    my_base_filename <- basename(my_h5_file) %>% tools::file_path_sans_ext()
    my_flightline_id <- my_h5_file %>% stringr::str_split('/') %>% unlist()
    my_reflectance_id <- glue('{my_flightline_id[9]}-{my_base_filename}')
    
    # reflectance and metadata
    my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
    my_refl <- my_h5[[glue('{my_aop_site}/Reflectance/Reflectance_Data')]]
    
    epsg_path <- glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/EPSG Code')
    my_epsg <- my_h5[[epsg_path]]$read()
    my_epsg2 <- glue('EPSG:{my_epsg}')
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
    my_rast <- terra::rast(ncols = my_dims[2],
                           nrows = my_dims[3],
                           crs = my_epsg2,
                           extent = my_extent)
  
    
    # get spatial location of pts in raster dimensions
    my_cellids <- terra::cells(my_rast, my_pts_spatvec)[,2] # modify this for polygons
    my_cellrowcols <- terra::rowColFromCell(my_rast, my_cellids)
    my_cell_xys <- terra::xyFromCell(my_rast, my_cellids)
    
    cellinfo_df <- my_water_sf %>% 
      st_drop_geometry() %>%
      dplyr::mutate(cellid = my_cellids,
                    cellrow = my_cellrowcols[,1],
                    cellcol = my_cellrowcols[,2],
                    cellx = my_cell_xys[, 1],
                    celly = my_cell_xys[,2]) %>%
      dplyr::filter(!is.nan(cellid)) %>%
      dplyr::mutate(reflectfile = my_reflectance_id)
    
    if(nrow(cellinfo_df)<1){message(glue('cell ids are all NaN for {my_reflectance_id}'))}
    if(nrow(cellinfo_df)>0){
    
    # should first check whether pixels are in NA region of raster!! 
    # otherwise getting lots of -9999s for no reason
    # get reflectance data and metadata
      
      my_epsg <- my_h5[[epsg_path]]$read()
      
      solar_zenith_path <- glue('{my_aop_site}/Reflectance/Metadata/Logs/Solar_Zenith_Angle')
      my_solar_zenith <- my_h5[[solar_zenith_path]]$read() 

      wls_path <- glue('{my_aop_site}/Reflectance/Metadata/Spectral_Data/Wavelength')
      my_wls <- my_h5[[wls_path]]$read()
      
      scale_factor <- my_refl$attr_open('Scale_Factor')$read()
      na_value <- my_refl$attr_open('Data_Ignore_Value')$read()
      cloud_conditions <- my_refl$attr_open('Cloud conditions')$read()
      
      my_spectra_list <- purrr::map(1:nrow(cellinfo_df), ~my_refl[1:426, cellinfo_df[['cellcol']][.x], cellinfo_df[['cellrow']][.x]])
      names(my_spectra_list) <- cellinfo_df[['loctype']]
      
      spectra_df <- my_spectra_list %>% 
        as.data.frame() %>% 
        mutate(wl = my_wls) %>%
        mutate(band = glue('band_{str_pad(1:426, 3, "left", "0")}'))
      
      # sensor zenith path
      sensor_zenith_path <- glue('{my_aop_site}/Reflectance/Metadata/to-sensor_Zenith_Angle')
      # same as reflectance
      my_sensor_zenith <- my_h5[[sensor_zenith_path]]
      my_sensor_zenith_list <- purrr::map(1:nrow(cellinfo_df), ~my_sensor_zenith[cellinfo_df[['cellcol']][.x], cellinfo_df[['cellrow']][.x]])
      names(my_sensor_zenith_list) <- cellinfo_df[['loctype']]
      
      cellinfo_df <- cellinfo_df %>%
        mutate(sensor_zenith = unlist(my_sensor_zenith_list),
               solar_zenith = my_solar_zenith,
               clouds = cloud_conditions,
               reflectfile = my_reflectance_id,
               aop_site = my_aop_site,
               crs = my_epsg2)
      
      spectra_path <- glue('{data_dir}/L1-reflectance2/spectra/{my_aq_site}_{my_reflectance_id}.tsv')
      cellinfo_path <- glue('{data_dir}/L1-reflectance2/meta/{my_aq_site}_{my_reflectance_id}.tsv')
      spectra_df %>% vroom::vroom_write(spectra_path)
      cellinfo_df %>% vroom::vroom_write(cellinfo_path)
      

    message(glue('saved spectra for {nrow(cellinfo_df)} pts at {my_aq_site} from {my_reflectance_id}'))}
    hdf5r::h5close(my_h5)    
  }
}

#
#

