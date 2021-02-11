library(neonhs)
library(viridis)
library(raster)
library(hdf5r)
library(sf)
library(tidyverse)
library(glue)

shp_dir <- '../../DATA/spatial'
swchem_sf <- glue('{shp_dir}/swchem_sites/swchem_sites.shp') %>%
  st_read()

head(swchem_sf)
como_s2 <- dplyr::filter(swchem_sf, siteID == 'COMO')

ggplot(swchem_sf) + geom_sf()

niwo_2020_dir <- 'D:/2020/FullSite/D13/2020_NIWO_4/L3/Spectrometer/Reflectance'
niwo_2020_files <- fs::dir_ls(niwo_2020_dir, recurse = TRUE, glob = '*.h5')
path_to_file <- niwo_2020_files[1]

como_s2_prj <- como_s2 %>% st_transform(st_crs(my_extent))

pt_in_img <- function(path_to_file){
  my_hs <- hs_read(path_to_file, bands = c(1))
  my_extent <- raster::extent(my_hs) %>% 
    as("SpatialPolygons") %>% 
    st_as_sf() %>% st_set_crs(proj4string(my_hs))
  mat1 <- st_intersects(my_extent, como_s2_prj, sparse = FALSE)
  pts_in_extent <- which(apply(mat1, 1, any))
  pt_in_img <- my_extent[pts_in_extent,]
  message(glue('point in image extent is {nrow(pt_in_img)>0}'))
  return(nrow(pt_in_img)>0)
}

niwo_test_list <- niwo_2020_files[101:182] %>% purrr::map_lgl(~pt_in_img(.x))
any(niwo_test_list)

plot(my_hs$band1_382nm, col = cividis(100), axes = FALSE, box = FALSE)

mypts <- as(como_s2, "Spatial")

my_vals <- hs_extract_pts(path_to_file, pts = como_s2_prj, bands = 1:426)

hs_extract_pts
hs_wavelength
