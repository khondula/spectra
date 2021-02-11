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
my_file <- niwo_2020_files[c(101:182)[which(niwo_test_list)]]

my_hs <- hs_read(my_file, bands = c(1, 50, 100, 400))
plot(my_hs, col = cividis(100), axes = FALSE, box = FALSE)

plot(my_hs$band100_877nm, col = cividis(100), axes = FALSE, box = FALSE)
plot(as(como_s2_prj, "Spatial"), add = TRUE, col = "green", cex = 3)

como_s2_sp <- as(como_s2_prj, "Spatial")
my_vals <- hs_extract_pts(my_file, pts = como_s2_sp, bands = 1:426)
my_vals_df <- my_vals %>% st_as_sf() %>% 
  st_drop_geometry() %>% as_tibble() %>%
  pivot_longer(cols = starts_with('band'), 
               names_to = 'band', 
               values_to = 'reflectance') %>% 
  tidyr::separate(band, into = c('index', 'wavelength'), sep = "_") %>%
  dplyr::mutate(wavelength = parse_number(wavelength))

my_vals_df %>%
  ggplot(aes(wavelength, reflectance)) +
  geom_line(lwd = 1) +
  theme_bw() +
  scale_color_viridis() +
  ggtitle('COMO.AOS.S2 (NIWO) 2020')

ggsave('figs/como2020.png')
# hs_extract_pts
# hs_wavelength
