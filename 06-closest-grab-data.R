# get AOS data based on date and site
library(glue)

rad_times_dir <- '/Volumes/hondula/DATA/radiance-time'

my_aq_site <- 'BLWA'
my_aop_yr <- '2017'
my_loc_type <- 'buoy.c0'

times_df <- fs::dir_ls(rad_times_dir, glob = glue('*{my_aq_site}*{my_aop_yr}*')) %>%
  read_csv()

times_df

# will need to figure out which flightline
# has the best weather? 

my_time <- times_df[['datetime']][1]

# get grab sample data for site closest to my_time
suva254_path <- 'suva254_all.csv'
my_analyte <- 'SUVA254'

grab_path <- suva254_path
grab_path <- 'suva280_all.csv'

grab_paths <- list(SUVA254 = 'suva254_all.csv',
                   SUVA280 = 'suva280_all.csv',
                   doc = 'doc_all.csv',
                   tss = 'tss_all.csv')


get_closest_sample <- function(my_analyte){
  
grab_path <-grab_paths[[my_analyte]]

grab_df <- glue('../neon-sites/results/{grab_path}') %>% 
  read_csv() %>%
  dplyr::mutate(loc_type = substr(namedLocation, 10, nchar(namedLocation))) %>%
  dplyr::filter(siteID %in% my_aq_site, loc_type %in% my_loc_type) %>%
  filter(!is.na(analyteConcentration))

# ggplot(grab_df, aes(x = collectDate, y = analyteConcentration)) +
#   geom_vline(aes(xintercept = my_time), col = 'red') +
#   geom_line() +
#   geom_point() +
#   theme_bw() +
#   ylim(0, NA) +
#   ggtitle(glue('{my_aq_site} {my_aop_yr} {my_analyte}'))

# time difference is the number of days before (-) 
# or after (+) the flight time that the grab sampling occurred
time_diff_df <- grab_df %>%
  dplyr::mutate(time_diff_days = difftime(collectDate, my_time, units = 'day')) %>% 
  mutate(time_diff_days = round(time_diff_days, 1)) %>% 
  mutate(min_diff = time_diff_days == min(abs(time_diff_days))) %>%
  dplyr::select(-uid)

closest_sample <- time_diff_df %>% filter(min_diff)
return(closest_sample)
}

my_analytes <- c('SUVA254', 'SUVA280', 'doc', 'tss')
my_samps <- my_analytes %>% purrr::map_dfr(~get_closest_sample(.x))
