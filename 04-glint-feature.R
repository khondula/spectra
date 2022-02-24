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
# my_aq_polygon <- 'PRLA_2016'
# polygon_file <- glue::glue('{data_dir}/AOP/site-polygons/{my_aq_polygon}.shp')
# my_water_sf <- sf::st_read(polygon_file)

my_aq_site <- 'CRAM'
my_aop_yr <- '2016'
my_aop_site <- 'UNDE'
my_aq_polygon <- 'CRAM_2015'
my_domain <- 'D05'

# note this assumes that all flightlines are in the same CRS for a siteyear
# need to redo BLUE to work with s2 loctype

save_d760_map('BARC', '2014', 'OSBS', 'BARC_2018', 'D03')
save_d760_map('BARC', '2016', 'OSBS', 'BARC_2018', 'D03')
save_d760_map('BARC', '2017', 'OSBS', 'BARC_2018', 'D03')
save_d760_map('BARC', '2018', 'OSBS', 'BARC_2018', 'D03')
save_d760_map('BARC', '2019', 'OSBS', 'BARC_2018', 'D03')

save_d760_map('SUGG', '2014', 'OSBS', 'SUGG_2018', 'D03')
save_d760_map('SUGG', '2016', 'OSBS', 'SUGG_2018', 'D03')
save_d760_map('SUGG', '2017', 'OSBS', 'SUGG_2018', 'D03')
save_d760_map('SUGG', '2018', 'OSBS', 'SUGG_2018', 'D03')
save_d760_map('SUGG', '2019', 'OSBS', 'SUGG_2018', 'D03')

save_d760_map('PRPO', '2016', 'WOOD', 'PRPO_2016', 'D09')
save_d760_map('PRPO', '2017', 'WOOD', 'PRPO_2016', 'D09')
save_d760_map('PRPO', '2019', 'WOOD', 'PRPO_2016', 'D09')
save_d760_map('PRPO', '2020', 'WOOD', 'PRPO_2016', 'D09')

save_d760_map('PRLA', '2016', 'WOOD', 'PRLA_2016', 'D09')
save_d760_map('PRLA', '2017', 'WOOD', 'PRLA_2016', 'D09')
save_d760_map('PRLA', '2019', 'WOOD', 'PRLA_2016', 'D09')
save_d760_map('PRLA', '2020', 'WOOD', 'PRLA_2016', 'D09')

save_d760_map('CRAM', '2016', 'UNDE', 'CRAM_2015', 'D05')
save_d760_map('CRAM', '2017', 'UNDE', 'CRAM_2015', 'D05')
save_d760_map('CRAM', '2019', 'UNDE', 'CRAM_2015', 'D05')
save_d760_map('CRAM', '2020', 'UNDE', 'CRAM_2015', 'D05')

save_d760_map('LIRO', '2017', 'LIRO', 'LIRO_2016', 'D05')
save_d760_map('LIRO', '2020', 'LIRO', 'LIRO_2016', 'D05')

save_d760_map('TOOK', '2017', 'TOOL', 'TOOK_2016', 'D18')
save_d760_map('TOOK', '2018', 'TOOL', 'TOOK_2016', 'D18')
save_d760_map('TOOK', '2019', 'TOOL', 'TOOK_2016', 'D18')

save_d760_map('BLWA', '2015', 'DELA', 'BLWA_2019', 'D08')
save_d760_map('BLWA', '2016', 'DELA', 'BLWA_2019', 'D08')
save_d760_map('BLWA', '2017', 'DELA', 'BLWA_2019', 'D08')
save_d760_map('BLWA', '2018', 'DELA', 'BLWA_2019', 'D08')
save_d760_map('BLWA', '2019', 'DELA', 'BLWA_2019', 'D08')

save_d760_map('FLNT', '2014', 'JERC', 'FLNT_2017', 'D03')
save_d760_map('FLNT', '2016', 'JERC', 'FLNT_2017', 'D03')
save_d760_map('FLNT', '2017', 'JERC', 'FLNT_2017', 'D03')
save_d760_map('FLNT', '2018', 'JERC', 'FLNT_2017', 'D03')
save_d760_map('FLNT', '2019', 'JERC', 'FLNT_2017', 'D03')

save_d760_map('TOMB', '2015', 'LENO', 'TOMB_2017', 'D08')
save_d760_map('TOMB', '2016', 'LENO', 'TOMB_2017', 'D08')
save_d760_map('TOMB', '2017', 'LENO', 'TOMB_2017', 'D08')
save_d760_map('TOMB', '2018', 'LENO', 'TOMB_2017', 'D08')
save_d760_map('TOMB', '2019', 'LENO', 'TOMB_2017', 'D08')


save_d760_map <- function(my_aq_site, my_aop_yr, my_aop_site, my_aq_polygon, my_domain){
  
  data_dir <- '/Volumes/hondula/DATA'
  polygon_file <- glue::glue('{data_dir}/AOP/site-polygons/{my_aq_polygon}.shp')
  my_water_sf <- sf::st_read(polygon_file)
  aos_pts_sf <- sf::st_read(glue::glue('{data_dir}/AOP/site-polygons/{my_aq_site}_AOSpts.shp'))

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
  aos_pts_sf_prj <- aos_pts_sf %>% sf::st_transform(my_epsg2)
  my_aos_spatvec <- aos_pts_sf_prj %>% as('SpatVector')
  my_polygon_spatvec <- my_aq_prj %>% as('SpatVector')

  # FOR each of the flightline files
  
  for(i in 1:length(my_site_files)){
    my_h5_file <- my_site_files[i]
    my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
    my_base_filename <- basename(my_h5_file) %>% tools::file_path_sans_ext()
    my_flightline_id <- my_h5_file %>% stringr::str_split('/') %>% unlist()
    my_reflectance_id <- glue('{my_base_filename}')
    
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
    my_aos_cellids <- terra::cells(my_rast, my_aos_spatvec)[,2]
    if(!any(!is.nan(my_cellids))){message(glue('cell ids in {my_aq_polygon} are all NaN for {my_reflectance_id}'))}
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
    
    aos_cellinfo_df <- cellinfo_df %>% filter(cellid %in% my_aos_cellids)
    my_rows <- range(cellinfo_df$cellrow)
    my_cols <- range(cellinfo_df$cellcol)
    my_cellx_range <- range(cellinfo_df$cellx)
    my_celly_range <- range(cellinfo_df$celly)
    
    # get reflectance data and metadata
      
      my_epsg <- my_h5[[epsg_path]]$read()
      
      solar_zenith_path <- glue('{my_aop_site}/Reflectance/Metadata/Logs/Solar_Zenith_Angle')
      my_solar_zenith <- my_h5[[solar_zenith_path]]$read() %>% round()
      
      wls_path <- glue('{my_aop_site}/Reflectance/Metadata/Spectral_Data/Wavelength')
      my_wls <- my_h5[[wls_path]]$read()
      
      scale_factor <- my_refl$attr_open('Scale_Factor')$read()
      na_value <- my_refl$attr_open('Data_Ignore_Value')$read()
      cloud_conditions <- my_refl$attr_open('Cloud conditions')$read()
      
      solar_azimuth_path <- glue('{my_aop_site}/Reflectance/Metadata/Logs/Solar_Azimuth_Angle')
      solar_azimuth <- my_h5[[solar_azimuth_path]]$read() %>% round()
      flight_heading_path <- glue('{my_aop_site}/Reflectance/Metadata/Flight_Trajectory/Flight_Heading')
      my_fh <- my_h5[[flight_heading_path]]
      fh_range <- range(my_fh[1,]) %>% round()

      # hdf5r::list.objects(my_h5) # everything? 59 things
      
      # wl 76 is 759, closest to oxygen absorption feature
      # 75 = 754, 77 = 764
      # get glint_wl,  + 1 and 
      glint_wl <- which.min(abs(760-my_wls))
      glint_bands <- c(glint_wl-1, glint_wl, glint_wl+1)
      # make raster with 3 bands for glint
      # wl760_rast <- terra::rast(my_refl[76, my_cols[1]:my_cols[2], my_rows[1]:my_rows[2]])
      # plot(wl760_rast)
      my_spectra_list <- purrr::map(glint_bands,
                                    ~terra::rast(my_refl[.x, my_cols[1]:my_cols[2], my_rows[1]:my_rows[2]]))
      my_spectra_list <- purrr::map(my_spectra_list, ~t(.))
      my_rast3 <- rast(my_spectra_list)
      my_rast3[my_rast3 == na_value] <- NA
      names(my_rast3) <- c('wl1', 'wl2', 'wl3')
      # plot(my_rast3)
      my_extent_crop <- terra::ext(my_cellx_range[1], my_cellx_range[2], my_celly_range[1], my_celly_range[2])
      terra::ext(my_rast3) <- my_extent_crop
      terra::crs(my_rast3) <- my_epsg2
      # plot(my_rast3$wl1)
      
      # my_aq_prj
      my_aq_prj_buff <- st_buffer(my_aq_prj, 0)
      # ggplot() +
      #   geom_stars(data = st_as_stars(my_rast3$wl1)) +
      #   geom_sf(data = my_aq_prj, fill = NA, col = 'green') +
      #   # geom_sf(data = my_aq_prj_buff, fill = NA, col = 'green') +
      #   ggtitle('testing glint stuff d760')
      
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
      my_d760_mask[my_d760_mask<0] <- 0
      my_minmax <- terra::minmax(my_d760_mask)
      my_range <- my_minmax[2,] - my_minmax[1,]
      
      gg <- ggplot() +
        # geom_stars(data = st_as_stars(my_d760_mask/1000)) +
        geom_stars(data = st_as_stars((my_d760_mask - my_minmax[1,])/my_range)) +
        geom_sf(data = my_aq_prj, fill = NA, col = 'cyan') +
        geom_sf(data = aos_pts_sf_prj, col = 'green', pch = 4, size = 5) +
        scale_fill_viridis_c(option = 'magma', 'd760') +
        # scale_fill_continuous_diverging(palette = 'Tropic', mid = 0.25) +
        # scale_fill_continuous_diverging(palette = 'Purple-Green', mid = 0, limits = c(-0.05, NA)) +
        ggtitle(glue('d760 in {my_aq_polygon} \n{my_base_filename}\n{my_wls[glint_wl]} nm\nHeading: {fh_range[1]} to {fh_range[2]}\nphi sun: {solar_azimuth}\ntheta sun: {my_solar_zenith}')) +
        theme_void()
      # gg
      gg2 <- gg + geom_sf_label(data = aos_pts_sf_prj, aes(label = loctype))
      ggsave(glue('figs/d760/{my_aq_site}-{my_aop_yr}-{my_base_filename}-d760.pdf'), gg, height = 20, width = 20)
      ggsave(glue('figs/d760/{my_aq_site}-{my_aop_yr}-{my_base_filename}-d760-labs.pdf'), gg2, height = 20, width = 20)
      # my_spectra_list <- purrr::map(1:nrow(cellinfo_df), ~my_refl[glint_bands, cellinfo_df[['cellcol']][.x], cellinfo_df[['cellrow']][.x]])
      # my_spectra_list[1]
      # make into a rast for plotting...
      
     
    hdf5r::h5close(my_h5)    
  }
  }
}

#
#

