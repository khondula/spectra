# plot some spectra
library(tidyverse)
library(vroom)
library(glue)

spectra_dir <- '/Volumes/hondula/DATA/spectra_site'

refl_dir <- '/Volumes/hondula/DATA/AOP/Reflectance'
tomb_h5 <- fs::dir_ls(glue('{refl_dir}/2018/TOMB'))
sugg_h5 <- fs::dir_ls(glue('{refl_dir}/2018/SUGG'))

library(neonhs)

wls <- hs_wavelength(tomb_h5[1], bands = 1:426)
wls[100]
wls[195]
wls[200]

fs::dir_ls(spectra_dir)

barc2019 <- vroom::vroom(glue('{spectra_dir}/BARC_2019.tsv'))
head(barc2019)

save(barc2019, file = "barc2019.RData")
