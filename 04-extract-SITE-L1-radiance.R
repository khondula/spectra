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

# data_dir <- '/Volumes/hondula/DATA'
data_dir <- 'H:/DATA/'

my_aq_site <- 'SUGG'
my_aop_yr <- '2014'
my_aop_site <- 'OSBS'
my_aq_polygon <- 'SUGG_AOSpts'
my_domain <- 'D03'


save_radiance_gpstime('SUGG', '2014', 'OSBS', 'SUGG_AOSpts', 'D03')
save_radiance_gpstime('BARC', '2014', 'OSBS', 'BARC_AOSpts', 'D03')
save_radiance_gpstime('CRAM', '2016', 'UNDE', 'CRAM_AOSpts', 'D05')


flightlines_df <- read_csv('results/l1-flightlines-wDomain.csv') %>%
  dplyr::select(aq_site, aop_yr, aop_site, shp, domainID) %>% distinct()

1:nrow(flightlines_df) %>% 
  purrr::walk(~save_radiance_gpstime(flightlines_df$aq_site[.x], 
                                    flightlines_df$aop_yr[.x], 
                                    flightlines_df$aop_site[.x], 
                                    flightlines_df$shp[.x], 
                                    flightlines_df$domainID[.x]))


save_radiance_gpstime <- function(my_aq_site, my_aop_yr, my_aop_site, 
                         my_aq_polygon, my_domain){
  
  data_dir <- 'H:/DATA/'
  polygon_file <- glue::glue('{data_dir}/AOP/site-polygons/{my_aq_polygon}.shp')
  my_water_sf <- sf::st_read(polygon_file)
  out_dir <- glue::glue('{data_dir}/L1-radiance')
  
  # first need to find which h5 files have points
  my_site_dir <- glue('D:/{my_aop_yr}/FullSite/{my_domain}') %>%
    fs::dir_ls(glob = glue('*{my_aop_site}*'), type = 'directory')
  my_site_files <- glue('{my_site_dir}/L1/Spectrometer/RadianceH5') %>%
    fs::dir_ls(glob = '*.h5', recurse = TRUE)
  
  # get the projection info from the first file
  my_h5_file <- my_site_files[1]
  my_base_filename <- basename(my_h5_file) %>% tools::file_path_sans_ext()
  my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
  epsg_path <- glue('{my_aop_site}/Radiance/Metadata/Coordinate_System/EPSG Code')
  my_epsg <- my_h5[[epsg_path]]$read() %>% as.integer()
  my_epsg2 <- glue('EPSG:{my_epsg}')
  hdf5r::h5close(my_h5)
  
  # convert points to h5 projection
  my_aq_prj <- my_water_sf %>% sf::st_transform(my_epsg2)
  my_pts_spatvec <- my_aq_prj %>% as('SpatVector')

  # FOR each of the flightline files
  
  for(i in 1:length(my_site_files)){
    my_h5_file <- my_site_files[i]
    my_base_filename <- basename(my_h5_file) %>% tools::file_path_sans_ext()
    my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
    
    # radiance and metadata
    my_obs <- my_h5[[glue('{my_aop_site}/Radiance/Metadata/Ancillary_Rasters/OBS_Data')]]
    # my_rad <- my_h5[[glue('{my_aop_site}/Radiance/RadianceIntegerPart')]]
    # my_rad$attr_open('Dimension_Labels')$read() # line, sample, wavelength
    # my_rad$attr_open('Dimensions')$read()
    # my_obs$attr_open('Dimension_Labels')$read() # line, sample, band
    my_dims <- my_obs$attr_open('Dimensions')$read() # rows, cols, bands
    na_value <- my_obs$attr_open('Data_Ignore_Value')$read()
    
    # extent
    map_info_path <- glue('{my_aop_site}/Radiance/Metadata/Coordinate_System/Map_Info')
    map_info <- my_h5[[map_info_path]]$read() %>% strsplit(',') %>% unlist()
    # map origin = mapinfo 4 and 5, xMin and yMax
    xy_resolution <- as.numeric(c(map_info[2], map_info[3]))
    xmin <- as.numeric(map_info[4])
    xmax <- xmin + my_dims[2] * xy_resolution[1]
    ymax <- as.numeric(map_info[5])
    ymin <- ymax - my_dims[1] * xy_resolution[2]
    my_extent <- terra::ext(xmin, xmax, ymin, ymax)
    
    # blank raster with map info
    my_rast <- terra::rast(ncols = my_dims[2],
                           nrows = my_dims[1],
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
      dplyr::mutate(radfile = my_base_filename)
    
    if(nrow(cellinfo_df)<1){message(glue('cell ids are all NaN for {my_base_filename}'))}
    if(nrow(cellinfo_df)>0){
    
    # should first check whether pixels are in NA region of raster!! 
    # otherwise getting lots of -9999s for no reason
    # get gps time 
    my_gps_list <- purrr::map(1:nrow(cellinfo_df), ~my_obs[10, cellinfo_df[['cellcol']][.x], cellinfo_df[['cellrow']][.x]])
    names(my_gps_list) <- cellinfo_df[['loctype']]
    
    gsp_df <- my_gps_list %>% 
      as.data.frame() %>%
      tidyr::pivot_longer(cols = names(my_gps_list), names_to = 'loctype', values_to = 'gpstime_hrs')
    
    cellinfo_df <- cellinfo_df %>%
      mutate(aop_site = my_aop_site,
             crs = my_epsg2) %>%
      left_join(gsp_df, by = 'loctype')
    
    cellinfo_path <- glue('{data_dir}/L1-radiance/meta/{my_aq_site}_{my_base_filename}.csv')
    cellinfo_df %>% vroom::vroom_write(cellinfo_path, delim = ',')
    
    message(glue('saved GPS time for {nrow(cellinfo_df)} pts at {my_aq_site} from {my_base_filename}'))}
    hdf5r::h5close(my_h5)    
  }
}

#
#

