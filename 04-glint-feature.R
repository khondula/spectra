# extract spectra from pixels within a polygon
# devtools::install_github('khondula/neonhs', force = TRUE)

# library(neonhs)
# library(viridis)
# library(raster)
library(colorspace)
library(hdf5r)
library(sf)
library(tidyverse)
library(glue)
library(terra)
library(stars)
# library(scico)
# library(vroom)
# source('R/myfxns.R')

# set up
# spectra_dir <- 'H:/DATA/spectra_site/'

data_dir <- '/Volumes/hondula/DATA'
fs::dir_ls(glue('{data_dir}/AOP/site-polygons'))
my_aq_polygon <- 'PRLA_2016'
polygon_file <- glue::glue('{data_dir}/AOP/site-polygons/{my_aq_polygon}.shp')
my_water_sf <- sf::st_read(polygon_file)

my_aq_site <- 'PRLA'
my_aop_yr <- '2017'
my_aop_site <- 'WOOD'
# my_aq_polygon <- 'BARC_2018'
my_domain <- 'D09'

# note this assumes that all flightlines are in the same CRS for a siteyear
# need to redo BLUE to work with s2 loctype

# save_L1reflectance_spectra <- function(my_aq_site, my_aop_yr, my_aop_site, 
#                          my_aq_polygon, my_domain){
  
  # data_dir <- /Volumes/hondula/DATA'
  polygon_file <- glue::glue('{data_dir}/AOP/site-polygons/{my_aq_polygon}.shp')
  my_water_sf <- sf::st_read(polygon_file)
  # out_dir <- glue::glue('{data_dir}/L1-reflectance2')
  
  if(data_dir == '/Volumes/hondula/DATA'){
    my_site_dir <- glue('{data_dir}/AOP/ReflectanceL1/{my_aop_yr}/{my_aop_site}')
    my_site_files <- glue('{my_site_dir}') %>%
      fs::dir_ls(glob = '*.h5', recurse = TRUE)
  }
  
  # get the projection info from the first file
  my_h5_file <- my_site_files[1]
  my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
  epsg_path <- glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/EPSG Code')
  my_epsg <- my_h5[[epsg_path]]$read() %>% as.integer()
  my_epsg2 <- glue('EPSG:{my_epsg}')
  hdf5r::h5close(my_h5)
  
  # convert polygon to h5 projection
  my_aq_prj <- my_water_sf %>% sf::st_transform(my_epsg2)
  my_polygon_spatvec <- my_aq_prj %>% as('SpatVector')

  # FOR each of the flightline files
  
  # for(i in 1:length(my_site_files)){
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
    my_cellids <- terra::cells(my_rast, my_polygon_spatvec)[,2] # modify this for polygons
    
    if(!any(!is.nan(my_cellids))){message(glue('cell ids are all NaN for {my_reflectance_id}'))}
    if(any(!is.nan(my_cellids))){
      
    my_cellrowcols <- terra::rowColFromCell(my_rast, my_cellids)
    my_cell_xys <- terra::xyFromCell(my_rast, my_cellids)
    
    cellinfo_df <- data.frame(cellid = my_cellids,
                    cellrow = my_cellrowcols[,1],
                    cellcol = my_cellrowcols[,2],
                    cellx = my_cell_xys[, 1],
                    celly = my_cell_xys[,2]) %>%
      dplyr::filter(!is.nan(cellid)) %>%
      dplyr::mutate(polygon_id = my_aq_polygon,
                     reflectfile = my_reflectance_id)
    
    my_rows <- range(cellinfo_df$cellrow)
    my_cols <- range(cellinfo_df$cellcol)
    my_cellx_range <- range(cellinfo_df$cellx)
    my_celly_range <- range(cellinfo_df$celly)
    
    # get reflectance data and metadata
      
      my_epsg <- my_h5[[epsg_path]]$read()
      
      solar_zenith_path <- glue('{my_aop_site}/Reflectance/Metadata/Logs/Solar_Zenith_Angle')
      my_solar_zenith <- my_h5[[solar_zenith_path]]$read() 

      wls_path <- glue('{my_aop_site}/Reflectance/Metadata/Spectral_Data/Wavelength')
      my_wls <- my_h5[[wls_path]]$read()
      
      scale_factor <- my_refl$attr_open('Scale_Factor')$read()
      na_value <- my_refl$attr_open('Data_Ignore_Value')$read()
      cloud_conditions <- my_refl$attr_open('Cloud conditions')$read()
      
      # wl 76 is 759, closest to oxygen absorption feature
      # 75 = 754, 77 = 764
      # get glint_wl,  + 1 and 
      glint_wl <- which.min(abs(760-my_wls))
      glint_bands <- c(glint_wl-1, glint_wl, glint_wl+1)
      
      # make raster with 3 bands for glint
      wl760_rast <- terra::rast(my_refl[76, my_cols[1]:my_cols[2], my_rows[1]:my_rows[2]])
      plot(wl760_rast)
      my_spectra_list <- purrr::map(glint_bands,
                                    ~terra::rast(my_refl[.x, my_cols[1]:my_cols[2], my_rows[1]:my_rows[2]]))
      my_spectra_list <- purrr::map(my_spectra_list, ~t(.))
      my_rast3 <- rast(my_spectra_list)
      my_rast3[my_rast3 == -9999] <- NA
      names(my_rast3) <- c('wl1', 'wl2', 'wl3')
      plot(my_rast3)
      my_extent_crop <- terra::ext(my_cellx_range[1], my_cellx_range[2], my_celly_range[1], my_celly_range[2])
      terra::ext(my_rast3) <- my_extent_crop
      terra::crs(my_rast3) <- my_epsg2
      plot(my_rast3$wl1)
      
      # my_aq_prj
      my_aq_prj_buff <- st_buffer(my_aq_prj, 10)
      ggplot() +
        geom_stars(data = st_as_stars(my_rast3$wl1)) +
        geom_sf(data = my_aq_prj, fill = NA, col = 'green') +
        # geom_sf(data = my_aq_prj_buff, fill = NA, col = 'green') +
        ggtitle('testing glint stuff d760')
      
      # ggsave('figs/barc-wl75.png')  
      

      d760 <- function(img, ii, jj, kk){
        bi <- img[[ii]]
        bj <- img[[jj]]
        bk <- img[[kk]]
        d760 <- (bi + bk)/2 -bj
        return(d760)
      }
      my_d760 <- d760(my_rast3, 1, 2, 3)
      my_aq_prj_buff_spat <- my_aq_prj_buff %>% as('SpatVector')
      my_d760_mask <- terra::mask(my_d760, my_aq_prj_buff_spat, inverse = FALSE)
      # my_d760_mask[my_d760_mask<0] <- NA
      ggplot() +
        geom_stars(data = st_as_stars(my_d760_mask/1000)) +
        geom_sf(data = my_aq_prj, fill = NA, col = 'blue') +
        scale_fill_viridis_c() +
        # scale_fill_continuous_diverging(palette = 'Purple-Green', mid = 0, limits = c(-0.05, NA)) +
        ggtitle('testing glint stuff d760 ')
      
      ggsave(glue('figs/{my_aq_site}-{my_aop_yr}-d760.pdf'), height = 20, width = 20)
      # my_spectra_list <- purrr::map(1:nrow(cellinfo_df), ~my_refl[glint_bands, cellinfo_df[['cellcol']][.x], cellinfo_df[['cellrow']][.x]])
      my_spectra_list[1]
      # make into a rast for plotting...
      r <- raster::raster(my_d760_mask)
      r <- raster::raster(my_d760)
      library(leaflet)
      library(leafem)

      leaflet(r) %>%
        addProviderTiles(providers$Esri.WorldImagery) %>%
        addRasterImage(r, group = 'd760', layerId = "wl1") %>%
        addImageQuery(r, layerId = "wl1") %>%
        addLayersControl(overlayGroups = c('d760'))
      
      leaflet(r) %>%
        addProviderTiles(providers$Esri.WorldImagery) %>%
        addRasterImage(r$wl1) %>%
        addImageQuery(r$wl1)

      #
      #
      #
      #
      
      names(my_spectra_list) <- cellinfo_df$cellid
      
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

