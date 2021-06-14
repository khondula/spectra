# library(stars)
library(hdf5r)
library(glue)
library(terra)
library(tidyverse)
library(sf)
library(neonhs)

# nbands <- 10
my_aq_site <- 'BARC'
my_aop_yr <- '2017'
my_aop_site <- 'OSBS'
my_domain <- 'D03'
my_aq_polygon <- 'BARC_AOSpts'

# l1 file path
my_h5_file <- 'D:/2017/FullSite/D03/2017_OSBS_3/L1/Spectrometer/ReflectanceH5/2017092713/NEON_D03_OSBS_DP1_20170927_170305_reflectance.h5'

my_dims <- neonhs::hs_dims(my_h5_file)
my_extent <- neonhs::hs_extent(my_h5_file)
my_epsg <- neonhs::hs_epsg(my_h5_file)
my_epsg <- glue('epsg:{my_epsg}')
my_rast <- terra::rast(nrow = my_dims[2],
                       ncol = my_dims[3],
                       crs = my_epsg,
                       extent = my_extent)

my_water_sf
my_aq_prj <- my_water_sf %>% st_transform(crs = my_epsg)
my_pts_spatvec <- my_aq_prj %>% as('SpatVector')

# get spatial location of pts in raster dimensions
my_cellids <- cells(my_rast, my_pts_spatvec)[,2] # modify this for polygons
my_cellrowcols <- rowColFromCell(my_rast, my_cellids)
my_cellrowcols
my_cell_xys <- xyFromCell(my_rast, my_cellids)
my_cellrowcols

my_index <- list(1:426, my_cellrowcols[,1], my_cellrowcols[,2])

my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
my_refl <- my_h5[[glue('{my_aop_site}/Reflectance/Reflectance_Data')]]
scale_factor <- my_refl$attr_open('Scale_Factor')$read()
na_value <- my_refl$attr_open('Data_Ignore_Value')$read()
cloud_conditions <- my_refl$attr_open('Cloud conditions')$read()
wls_path <- glue('{my_aop_site}/Reflectance/Metadata/Spectral_Data/Wavelength')
my_wls <- my_h5[[wls_path]]$read()

my_spectra_list <- purrr::map(1:length(my_cellids), ~my_refl[1:426,
                                          my_cellrowcols[.x, 1],
                                          my_cellrowcols[.x, 2]])

names(my_spectra_list) <- my_aq_prj$loctype
my_spectra_list %>% as.data.frame() %>% 
  mutate(wl = my_wls) %>% 
  tidyr::pivot_longer(cols = 1:3, names_to = 'loctype') %>%
  ggplot(aes(x = wl, y = value/scale_factor)) +
  geom_line(aes(col = loctype)) +
  geom_point(aes(col = loctype)) +
  xlim(c(380, 800)) + ylim(c(0, 0.1))

ggsave('barc-2017-zoom.png')  



### old below here


my_spectra <- my_refl[1:426, my_cellrowcols[1,1], my_cellrowcols[1,2]]



my_spectra <- terra::extract(my_hs_terra, my_pts_spatvec, cells = TRUE)
my_spectra_df <- my_spectra %>% 
  as_tibble() %>%
  mutate(x = terra::xFromCell(my_hs_terra, cell),
         y = terra::yFromCell(my_hs_terra, cell))

# make empty rast with extent, dimensions, epsg 
# find xy index of cells in pts
# get x and y of those points
# read in data values of spectra
# read in data values of sensor zenith angle

epsg_path <- glue::glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/EPSG Code')
my_epsg <- file_h5[[epsg_path]]$read() %>% as.integer()

terra::crs(b10_rast) <- sf::st_crs(paste0("epsg:" ,my_epsg)) %>% as.character()
terra::ext(b10_rast) <- neonhs::hs_extent(my_h5_file)

index <- list(1:426, seq_len(dims[2]), seq_len(dims[3]))

wls_path <- glue('{my_aop_site}/Reflectance/Metadata/Spectral_Data/Wavelength')
my_wls <- my_h5[[wls_path]]$read()
hdf5r::h5close(my_h5)

my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")

# my_h5_file <- 'D:/2017/FullSite/D03/2017_OSBS_3/L3/Spectrometer/Reflectance/NEON_D03_OSBS_DP3_402000_3283000_reflectance.h5'
# gdal_subdatasets(my_h5_file)
my_refl_file <- "HDF5:\"D:/2017/FullSite/D03/2017_OSBS_3/L1/Spectrometer/ReflectanceH5/2017092713/NEON_D03_OSBS_DP1_20170927_170305_reflectance.h5\"://OSBS/Reflectance/Reflectance_Data"
my_reflectance <- stars::read_stars(my_refl_file, 
                                    NA_value = -9999, 
                                    # RasterIO = list(x = 1:1046, y = 1:14058, bands = 1:426), 
                                    proxy = TRUE)

my_reflectance

my_refl <- st_set_dimensions(my_reflectance, xy = c('y', 'band')) %>% 
  st_set_dimensions(names = c('band', 'x', 'y'))

my_refl

my_refl <- st_set_dimensions(my_refl, "band", delta = NA)
my_refl <- st_set_dimensions(my_refl, "x", delta = 1)
my_refl <- st_set_dimensions(my_refl, "y", delta = -1)
# my_refl <- st_set_dimensions(my_refl, "x", offset = 1046) # ??


refl_path <- glue('{my_aop_site}/Reflectance/Reflectance_Data')

##### read in reflectance
my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
file_h5 <- hdf5r::H5File$new(my_h5_file, mode = 'r')
site <- file_h5$ls()$name
dims <- neonhs::hs_dims(my_h5_file)
index <- list(1:426, seq_len(dims[2]), seq_len(dims[3]))

reflectance <- file_h5[[paste0(site, '/Reflectance/Reflectance_Data')]]


read_refl <- function(x, reflectance){
  reflectance$read(args = list(x, 300:600, 1:10000))}
b10 <- purrr::map(1:100, ~read_refl(.x, reflectance))
  

b10_rast <- terra::rast(t(b10[[100]]))
b10_rast[b10_rast == -9999] <- NA
b10_rast <- b10_rast/10000
b10_rast %>% plot()

epsg_path <- glue::glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/EPSG Code')
my_epsg <- file_h5[[epsg_path]]$read() %>% as.integer()

terra::crs(b10_rast) <- sf::st_crs(paste0("epsg:" ,my_epsg)) %>% as.character()
terra::ext(b10_rast) <- neonhs::hs_extent(my_h5_file)
#########


#


arrays <- lapply(bands, function(x) {
  terra::rast(reflectance$read(args = list(x, index[[2]], index[[3]])))
})


my_refl2 <- my_h5[[refl_path]]
b10 <- my_refl2$read()
b10[b10==-9999] <- NA

plot(b10$lyr.1)

epsg_path <- glue('{my_aop_site}/Reflectance/Metadata/Coordinate_System/EPSG Code')
my_epsg <- my_h5[[epsg_path]]$read() %>% as.integer()


sf::gdal_read(my_h5_file)

# my_h5 <- hdf5r::H5File$new(my_h5_file, mode = "r")
# refl_data_path <- glue('{my_aop_site}/Reflectance/Reflectance_Data')
# my_reflectance <- my_h5[[refl_data_path]]

my_reflectance <- stars::read_stars(my_h5_file, proxy = TRUE)
gdal_subdatasets(my_h5_file)
hdf5r::h5attributes(my_reflectance)
my_h5$close_all()
neonhs::hs_dims(my_h5_file)

my_reflectance[1,1,426]

# instead of reading in whole reflectance data,
# use dims and map info to create a raster (terra?)
# and find location of points of interest
# read all bands of interest in from THOSE cells
# then don't forget to clean with NA value and scale factor

# my_stars <- my_reflectance[1:426,1:1000,1:1000] %>% st_as_stars()
my_stars

plot(my_stars)

reflectance$read(args = list(1, 1, 1))

# my_stars <- my_h5[[refl_data_path]]$read() %>% stars::read_stars()
