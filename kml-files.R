glue('{shp_dir}/swchem_sites/swchem_sites.shp') %>%
  st_read() %>% 
  rename(Name = siteID, Description = nmdLctn) %>%
  st_write('my_aq_sf.kml')

my_kml <- '20190730_164923_hsi_kml_0000.kml'
my_kml_full <- fs::dir_ls(my_site_dir, 
                          glob = glue('*{my_kml}'), 
                          recurse = TRUE)

my_kmls <- fs::dir_ls(my_site_dir, glob = glue('*FlightlineBoundary*.kml'), recurse = TRUE)

library(xml2)
my_kml <- my_kmls[1]

get_flightline_color <- function(my_kml){
  mylist <- my_kml %>% xml2::read_xml() %>% xml2::as_list()
  return(mylist$kml$Document$Folder$Placemark$Style$LineStyle$color[[1]])
  }

my_cols <- my_kmls %>% purrr::map_chr(~get_flightline_color(.x))
kml_colors <- read_csv('kml-colors.csv', col_types = 'ccc')

my_cols %>% as.data.frame() %>% 
  tibble::rownames_to_column() %>%
  dplyr::mutate(filename = basename(rowname)) %>%
  dplyr::rename(kml_color = 2) %>%
  left_join(kml_colors, by = 'kml_color') %>%
  dplyr::select(filename, kml_color, color, condition) %>%
  head()
unique(my_cols)

get_flightline_color(my_kmls[68])
my_kmls[37] %>% get_flightline_color()
my_kmls[1] %>% get_flightline_color()


# red ff0000ff
# yellow ff00ffff
# green ff00ff00


flightline_sf <- st_read(my_kmls[68])

ggplot() +
  geom_sf(data = flightline_sf)
