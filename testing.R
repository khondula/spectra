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

my_aq_site <- 'BLDE'
my_aop_yr <- '2020'
my_aop_site <- 'YELL'
my_domain <- 'D12'



# save_spectra_df <- function(my_aq_site, my_aop_yr, my_aop_site, my_domain){
  
my_aq_sf <- dplyr::filter(swchem_sf, siteID == my_aq_site)

my_site_dir <- glue('D:/{my_aop_yr}/FullSite/{my_domain}') %>%
  fs::dir_ls(glob = glue("*{my_aop_site}*"), type = 'directory')

my_site_files <- glue('{my_site_dir}/L3/Spectrometer/Reflectance') %>%
  fs::dir_ls(glob = '*.h5')

# read in h5 file
# read in just one to project sampling point
my_hs <- hs_read(my_site_files[1], bands = c(1))
my_extent <- raster::extent(my_hs) %>% 
  as("SpatialPolygons") %>% 
  st_as_sf() %>% 
  st_set_crs(proj4string(my_hs))

# project sampling point
my_aq_prj <- my_aq_sf %>% st_transform(proj4string(my_hs))
path_to_file <- my_site_files[1]
my_prj <- proj4string(my_hs)

pt_in_img <- function(path_to_file){
  # my_hs <- hs_read(path_to_file, bands = c(1))
  # my_extent <- raster::extent(my_hs) %>% 
  my_extent <- hs_extent(path_to_file) %>%
    as("SpatialPolygons") %>% 
    st_as_sf() %>% 
    st_set_crs(my_prj)
  mat1 <- st_intersects(my_extent, my_aq_prj, sparse = FALSE)
  pts_in_extent <- which(apply(mat1, 1, any))
  pt_in_img <- my_extent[pts_in_extent,]
  return(nrow(pt_in_img)>0)
}

site_files_list <- my_site_files %>% purrr::map_lgl(~pt_in_img(.x))
my_file <- my_site_files[which(site_files_list)]
my_file <- my_file[1]
# multiple files if there are multiple points... 


# add in saving an image of where the point is? #
my_hs <- hs_read(my_file, bands = c(1, 50, 100, 400))
plot(my_hs[[3]], col = cividis(100), axes = FALSE, box = FALSE)
plot(as(my_aq_prj, "Spatial"), add = TRUE, col = "green", cex = 3)
my_aq_prj
# click_point <- raster::click(my_hs[[3]], n = 1, xy = TRUE)
# shen_river_pt <- st_point(c(click_point[['x']], click_point[['y']])) %>% 
#   st_sfc(crs = st_crs(my_aq_prj)) %>% st_as_sf()

# plot(my_hs[[3]], col = cividis(100), axes = FALSE, box = FALSE)
# plot(as(shen_river_pt, "Spatial"), add = TRUE, col = "green", cex = 3)

my_aq_sp <- as(my_aq_prj, "Spatial")
my_vals <- hs_extract_pts(my_file, pts = my_aq_sp, bands = 1:426)

spectra_dir <- 'H:/DATA/spectra/'

my_vals_df <- my_vals %>% st_as_sf() %>% 
  st_drop_geometry() %>% as_tibble() %>%
  pivot_longer(cols = starts_with('band'), 
               names_to = 'band', 
               values_to = 'reflectance') %>% 
  tidyr::separate(band, into = c('index', 'wavelength'), sep = "_") %>%
  dplyr::mutate(wavelength = parse_number(wavelength))

my_vals_df %>% write_csv(glue('{spectra_dir}/{my_aq_site}_{my_aop_yr}.csv'))

my_vals_df %>%
  ggplot(aes(wavelength, reflectance)) +
  geom_line(lwd = 1, aes(col = nmdLctn)) +
  coord_cartesian(ylim = c(0, 0.1),
                  xlim = c(380, 1000)) +
  theme_bw()
# }

# save_spectra_df(my_aq_site = 'PRLA', my_aop_yr = '2019', my_aop_site = 'WOOD', my_domain = 'D09')
# save_spectra_df(my_aq_site = 'ARIK', my_aop_yr = '2020', my_aop_site = 'ARIK', my_domain = 'D10')
# save_spectra_df(my_aq_site = 'BLDE', my_aop_yr = '2019', my_aop_site = 'YELL', my_domain = 'D12')
# save_spectra_df(my_aq_site = 'BLDE', my_aop_yr = '2020', my_aop_site = 'YELL', my_domain = 'D12')

ggsave('figs/shen-river2019.png')
# hs_extract_pts
# hs_wavelength
