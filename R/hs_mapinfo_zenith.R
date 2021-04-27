# function to get hs_extent from radiance data

# filename <- my_site_files[2]
# producttype <- 'Radiance'
# filename = my_h5_file
# producttype = 'Reflectance'

hs_mapinfo_zenith <- function(filename, producttype = 'Reflectance'){
  file_h5 <- hdf5r::H5File$new(filename, mode = 'r')
  site <- file_h5$ls()$name
  mapinfo_path <- glue::glue('{site}/{producttype}/Metadata/Coordinate_System/Map_Info')
  map_info <- file_h5[[mapinfo_path]]$read()
  map_info <- stringr::str_split(map_info, ',') %>% unlist()
  
  # add in line here to check for which name is used for radiance data
  # hdf5r::list.datasets(my_h5) 
  
  # dims_path <- glue::glue('{site}/{producttype}/Radiance_Data')
  # dims_path <- glue::glue('{site}/{producttype}/RadianceIntegerPart')
  dims_path <- glue::glue('{site}/{producttype}/Metadata/to-sensor_zenith_angle')
  
  my_dims <- file_h5[[dims_path]]$dims
  xy_res <- c(map_info[2] %>% as.numeric(),
              map_info[3] %>% as.numeric())
  file_h5$close_all()
  my_xmin = map_info[4] %>% as.numeric()
  my_xmax = my_xmin + my_dims[1] * xy_res[1] # only 2 dimensions for angles
  my_ymax = map_info[5] %>% as.numeric()
  my_ymin = my_ymax - my_dims[2] * xy_res[2] # only 2 dimensions for angles
  map_info_list <- list(xmin = my_xmin, 
                        xmax = my_xmax,
                        ymin = my_ymin,
                        ymax = my_ymax)
  return(map_info_list)
}

# my_site_files %>% purrr::map(~hs_mapinfo(.x, 'Radiance'))
