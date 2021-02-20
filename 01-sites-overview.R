## extracting spectra
library(neonhs)
library(viridis)
library(raster)
library(hdf5r)
library(sf)
library(tidyverse)
library(glue)
library(stars)
source('R/myfxns.R')

# set up
spectra_dir <- 'H:/DATA/spectra/'
shp_dir <- '../../DATA/spatial'

# inputs
my_aq_site <- 'BLDE'
my_aop_yr <- '2020'

save_my_plots <- function(my_aq_site, my_aop_yr){
  spectra_dir <- 'H:/DATA/spectra/'
  shp_dir <- '../../DATA/spatial'
  
   get_aop_dates <- function(aq_siteids){
    aop_file <- 'results/sites_join_aop_dates.csv'
    aop_dates <- readr::read_csv(aop_file, col_types = 'ccccccccddD') %>%
      dplyr::filter(siteID %in% aq_siteids) %>%
      dplyr::select(siteID, aop_site_id, flightdate) %>%
      dplyr::arrange(flightdate) %>% distinct()
    return(aop_dates)
  }
  
# AOP site name
my_aop_site <- get_aop_dates(my_aq_site) %>% 
  pull(aop_site_id) %>% unique()

# read in points
my_aq_sf <- glue('{shp_dir}/swchem_sites/swchem_sites.shp') %>%
  st_read() %>% dplyr::filter(siteID == my_aq_site) %>%
  mutate(location_type = substr(nmdLctn, 10, nchar(nmdLctn)))
my_loc_types <- my_aq_sf$location_type
my_domain <- my_aq_sf$domanID[1]

# folder with AOP data
my_site_dir <- glue('D:/{my_aop_yr}/FullSite/{my_domain}') %>%
  fs::dir_ls(glob = glue("*{my_aop_site}*"), type = 'directory')
my_site_files <- glue('{my_site_dir}/L3/Spectrometer/Reflectance') %>%
  fs::dir_ls(glob = '*.h5')

# read in first file to project sampling points
# read in h5 file
# read in just one to project sampling point
my_prj <- hs_epsg(my_site_files[1]) %>% as.integer()
my_aq_prj <- my_aq_sf %>% st_transform(crs = my_prj)

# identify which h5 file contains points SLOW
# will need to update to look for file in table
message('finding the right h5 file...')
pt_in_img <- function(path_to_file){
  my_extent <- hs_extent(path_to_file) %>%
    as("SpatialPolygons") %>% 
    st_as_sf() %>% 
    st_set_crs(my_prj)
  mat1 <- st_intersects(my_extent, my_aq_prj, sparse = FALSE)
  pts_in_extent <- which(apply(mat1, 1, any))
  pt_in_img <- my_extent[pts_in_extent,]
  return(nrow(pt_in_img)>0)
}

site_files_list <- my_site_files %>% 
  purrr::map_lgl(~pt_in_img(.x))
my_file <- my_site_files[which(site_files_list)]

file_df <- my_aq_prj %>% st_drop_geometry() %>%
  mutate(year = my_aop_yr) %>%
  left_join(data.frame(siteID = my_aq_site, h5_files = basename(my_file)),
            by = 'siteID')

# SAVE file name and site, year
file_df %>% write_csv('results/which-h5.csv', append = TRUE)

if(length(my_file) == 1){
  my_hs <- hs_read(my_file, bands = c(1, 50, 100, 400))
}
if(length(my_file) > 1){
  my_hs <- hs_read(my_file[1], bands = c(1, 50, 100, 400))
  my_hs_list <- purrr::map(my_file, ~hs_read(.x, bands = c(1, 50, 100, 400)))
  for(i in 2:length(my_file)){
    my_hs <- merge(my_hs, my_hs_list[[i]])
  }
}


# save showing where point is
my_aq_prj_buff100 <- st_buffer(my_aq_prj, 100) %>% st_union()
my_buff100_sp <- as(my_aq_prj_buff100, "Spatial")
my_aq_prj_buff30 <- st_buffer(my_aq_prj, 30) %>% st_union()
my_buff30_sp <- as(my_aq_prj_buff30, "Spatial")

# which bands to read in?
# whats a good COLOR PALETTE?

my_hs_crop100 <- raster::crop(my_hs, my_buff100_sp)
my_hs_crop30 <- raster::crop(my_hs, my_buff30_sp)

# stars plotting... maybe when I need more bands
# my_hs_stars <- stars::st_as_stars(my_hs)
# z <- my_hs_stars %>% dplyr::select(band1_382nm) %>% adrop()
# 
# ggplot() +
#   geom_stars(data = z[1]) +
#   geom_sf(data = my_aq_prj) +
#   facet_wrap("band") +
#   scale_fill_viridis(option = 'cividis')
# 
# plot(my_hs, col = cividis(100), axes = FALSE, box = FALSE)
my_sp_list <- my_loc_types %>% 
  purrr::map(~filter(my_aq_prj, location_type == .x)) %>%
  purrr::map(~as(.x, 'Spatial'))
names(my_sp_list) <- my_loc_types

message('making zoomed out image')

pdf(glue('figs/maps/{my_aq_site}-{my_aop_yr}.pdf'))
plot(my_hs[[3]], col = cividis(100), axes = FALSE, box = FALSE,
     main = glue('{my_aq_site} {my_aop_yr} \n {names(my_hs)[3]}'))
if('S2' %in% my_loc_types){
  plot(my_sp_list[['S2']], add = TRUE, col = "green", cex = 3, pch = 3)
}
if('buoy.c0' %in% my_loc_types){
  plot(my_sp_list[['buoy.c0']], add = TRUE, col = "green", cex = 2, pch = 3)
}
if('reach' %in% my_loc_types){
  plot(my_sp_list[['reach']], add = TRUE, col = "red", cex = 1, pch = 4)
}
if('inlet' %in% my_loc_types){
  plot(my_sp_list[['inlet']], add = TRUE, col = "hot pink", cex = 1, pch = 4)
}
if('outlet' %in% my_loc_types){
  plot(my_sp_list[['outlet']], add = TRUE, col = "cyan", cex = 1, pch = 4)
}
dev.off()

message('making 100m image')

pdf(glue('figs/maps/{my_aq_site}-{my_aop_yr}_100m.pdf'))
plot(my_hs_crop100[[3]], col = cividis(100), axes = FALSE, box = FALSE,
     main = glue('{my_aq_site} {my_aop_yr} \n {names(my_hs_crop100)[3]}'))
if('S2' %in% my_loc_types){
  plot(my_sp_list[['S2']], add = TRUE, col = "green", cex = 3, pch = 3)
}
if('buoy.c0' %in% my_loc_types){
  plot(my_sp_list[['buoy.c0']], add = TRUE, col = "green", cex = 2, pch = 3)
}
if('reach' %in% my_loc_types){
  plot(my_sp_list[['reach']], add = TRUE, col = "red", cex = 1, pch = 4)
}
if('inlet' %in% my_loc_types){
  plot(my_sp_list[['inlet']], add = TRUE, col = "hot pink", cex = 1, pch = 4)
}
if('outlet' %in% my_loc_types){
  plot(my_sp_list[['outlet']], add = TRUE, col = "cyan", cex = 1, pch = 4)
}
dev.off()

message('making 30m image')
pdf(glue('figs/maps/{my_aq_site}-{my_aop_yr}_30m.pdf'))
plot(my_hs_crop30[[3]], col = cividis(100), axes = FALSE, box = FALSE,
     main = glue('{my_aq_site} {my_aop_yr} \n {names(my_hs_crop30)[3]}'))
if('S2' %in% my_loc_types){
  plot(my_sp_list[['S2']], add = TRUE, col = "green", cex = 3, pch = 3)
}
if('buoy.c0' %in% my_loc_types){
  plot(my_sp_list[['buoy.c0']], add = TRUE, col = "green", cex = 2, pch = 3)
}
if('reach' %in% my_loc_types){
  plot(my_sp_list[['reach']], add = TRUE, col = "red", cex = 1, pch = 4)
}
if('inlet' %in% my_loc_types){
  plot(my_sp_list[['inlet']], add = TRUE, col = "hot pink", cex = 1, pch = 4)
}
if('outlet' %in% my_loc_types){
  plot(my_sp_list[['outlet']], add = TRUE, col = "cyan", cex = 1, pch = 4)
}
dev.off()

}

save_my_plots('LEWI', '2019')
save_my_plots('POSE', '2017')

save_my_plots('FLNT', '2018')
save_my_plots('SUGG', '2019')
save_my_plots('CUPE', '2018')
save_my_plots('GUIL', '2019')
save_my_plots('CRAM', '2020')
save_my_plots('LIRO', '2020')
save_my_plots('KING', '2020')
save_my_plots('MCDI', '2020')
save_my_plots('LECO', '2018')
save_my_plots('BWLA', '2019')
save_my_plots('MAYF', '2019')
save_my_plots('TOMB', '2019')
save_my_plots('MAYF', '2019')
save_my_plots('PRLA', '2020')
save_my_plots('PRPO', '2020')
save_my_plots('ARIK', '2020')
save_my_plots('PRIN', '2017')
save_my_plots('BLDE', '2020')
save_my_plots('WLOU', '2020')
save_my_plots('REDB', '2019')
save_my_plots('MART', '2019')
save_my_plots('MCRA', '2018')
save_my_plots('BIGC', '2019')
save_my_plots('TECR', '2019')
save_my_plots('CARI', '2019')






# extract spectra for point
# identify adjacent cells? 
# water mask... 