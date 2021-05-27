# https://gis.stackexchange.com/questions/364667/r-4-0-1-not-sure-i-understand-this-message-warning-message-in-proj4stringx

filename <- "D:/2017/FullSite/D09/2017_WOOD_2/L3/Spectrometer/Reflectance/NEON_D09_WOOD_DP3_490000_5222000_reflectance.h5"
bands <- 1:100

use_h5_extent <- TRUE
h5_extent <- hs_extent(filename)
out_extent <- h5_extent
dims <- hs_dims(filename)
index <- list(bands, seq_len(dims[2]), seq_len(dims[3]))

hs_clean <- function(vals, reflectance) {
    na_value <- reflectance$attr_open('Data_Ignore_Value')$read()
    vals[vals == na_value] <- NA
    scale_factor <- reflectance$attr_open('Scale_Factor')$read()
    vals <- vals / scale_factor
    vals
}
  
read_hs_values <- function(filename, index){
    bands <- index[[1]]
    stopifnot(all(bands > 0))
    n_band <- length(bands)
    
    file_h5 <- hdf5r::H5File$new(filename, mode = 'r+')
    site <- file_h5$ls()$name
    
    reflectance <- file_h5[[paste0(site, '/Reflectance/Reflectance_Data')]]
    
    arrays <- lapply(bands, function(x) {
      raster::raster(reflectance$read(args = list(x, index[[2]], index[[3]])))
    })
    
    r <- raster::t(raster::stack(arrays))
    r <- hs_clean(r, reflectance)
    names(r) <- hs_wavelength(filename, bands)
  
    r2 <- r %>% terra::rast()
    terra::crs(r2) <- st_crs(my_epsg) %>% as.character()
  }
  
r <- read_hs_values(my_h5_files[1], index)

raster::extent(r) <- out_extent
r