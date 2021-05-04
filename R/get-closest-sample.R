library(glue)
library(tidyverse)
# my_analyte <- 'chla'
# my_time <- my_radtime
get_closest_sample <- function(my_analyte, my_time){
  
  my_loc_type <- 'buoy.c0'
  grab_paths <- list(SUVA254 = 'suva254_all.csv',
                     SUVA280 = 'suva280_all.csv',
                     doc = 'doc_all.csv',
                     tss = 'tss_all.csv',
                     chla = 'all-phyto-data.csv')
  grab_path <- grab_paths[[my_analyte]]
  
  grab_df <- glue('~/Documents/projects/neon-sites/results/{grab_path}') %>% 
    read_csv() %>%
    dplyr::mutate(loc_type = substr(namedLocation, 10, nchar(namedLocation))) %>%
    dplyr::filter(siteID %in% my_aq_site, loc_type %in% my_loc_type) %>%
    filter(!is.na(analyteConcentration))
  
  # if(my_analyte %in% c('chla')){
  #   grab_df <- grab_df %>%
  #     dplyr::select(-collectDate) %>%
  #     dplyr::mutate(collectDate = lubridate::as_datetime(glue('{collect_date} {time_hms}')))
  # }
  
  # ggplot(grab_df, aes(x = collectDate, y = analyteConcentration)) +
  #   geom_vline(aes(xintercept = my_time), col = 'red') +
  #   geom_line() +
  #   geom_point() +
  #   theme_bw() +
  #   ylim(0, NA) +
  #   facet_wrap(vars(analyte), ncol = 1) +
  #   ggtitle(glue('{my_aq_site} {my_aop_yr} {my_analyte}'))
  
  # time difference is the number of days before (-) 
  # or after (+) the flight time that the grab sampling occurred
  time_diff_df <- grab_df %>%
    dplyr::mutate(time_diff_days = difftime(collectDate, my_time, units = 'day')) %>% 
    mutate(time_diff_days = round(time_diff_days, 1)) %>%
    mutate(time_diff_abs = abs(time_diff_days)) %>%
    mutate(min_diff = time_diff_abs == min(time_diff_abs))
  
  closest_sample <- time_diff_df %>% filter(min_diff)
  return(closest_sample)
}