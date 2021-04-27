# spectra in WASI

wasi_data_dir <- 'C:/Users/hondula/Downloads/WASI4/DATA'
fs::dir_ls(wasi_data_dir)

# Mean spectrum typical for Lake Constance (Bodensee)
# - Weighted sum of CRY-L.A, CRY-H.A, DIA.A, DINO.A, GREEN.A
# - Calculated by Heege using 
# (a) 213 in-situ phytoplankton absorption spectra of Gege from 32 days in 1990 and 1991
# Gege derived these by inverse modeling of his above-water remote sensing reflectance measurements
# (b) 32 measurements of chlorophyll-a and phaeophytin-a from 32 days in 1990 and 1991
# Measured at the University of Constance by Beese, Richter and Kenter
# - Validated by Heege using 139 irradiance reflectance and 278 attenuation measurements from 1990 to 1996


phyto_a <- glue('{wasi_data_dir}/PHYTO.A') %>%
  read_delim(delim = ' ', skip = 12,
             col_types = 'nn',
             col_names = c('nm', 'm2_mg')) %>%
  mutate(phyto_class = "phyto_mixA")

# Species:   CRYPTOMONAS OVATA
# Source:    Absorption measurement of Gabi Hartmann in the laboratory (13 October 1993) using a pure culture
# Reference: G. M. Hartmann (1995), Diplomarbeit, Albert-Ludwigs-Universität Freiburg, 91 pp.
# Modifications by Gege:
#   - Fit with a sum of Lorentzian functions
# - Normalisation to a*(440) = 0.0409 m^2/mg 

phyto_b <- glue('{wasi_data_dir}/CRY-LO.A') %>%
  read_delim(delim = ' ', skip = 10,
             col_types = 'nn',
             col_names = c('nm', 'm2_mg')) %>%
  mutate(phyto_class = "Cryptophyta_L")

# Species:   RHODOMONAS LENS
# Source:    Optical Density measurement of Haxo and Fork (1959)
# Reference: Haxo and Fork (1959), Photosynthetically active accessory pigments of cryptomonads, Nature 184, 1051-1052
# Modifications by Gege:
#   - Digitalisation of the published Optical Density measurement 
# - Fit with a sum of Lorentzian functions to get 1 nm intervals
# - Normalisation to a*(440) = 0.0228 m^2/mg


phyto_c <- glue('{wasi_data_dir}/CRY-HI.A') %>%
  read_delim(delim = ' ', skip = 11,
             col_types = 'nn',
             col_names = c('nm', 'm2_mg')) %>%
  mutate(phyto_class = "Cryptophyta_H")

# Species:   FRAGILLARIA CROTONENSIS
# Source:    Absorption measurement of Gabi Hartmann in the laboratory (25 October 1993) using a pure culture
# Reference: G. M. Hartmann (1995), Diplomarbeit, Albert-Ludwigs-Universität Freiburg, 91 pp.
# Modifications by Gege:
#   - Fit with a sum of Lorentzian functions
# - Normalisation to a*(440) = 0.036356 m^2/mg 

phyto_d <- glue('{wasi_data_dir}/DIA.A') %>%
  read_delim(delim = ' ', skip = 10,
             col_types = 'nn',
             col_names = c('nm', 'm2_mg')) %>%
  mutate(phyto_class = "Diatoms")

# Specific absorption of Dinoflagellates (m^2/mg)
# 
# Species:   mixture of five marine dinoflagellates
# Source:    Absorption measurement of von Smekot-Wensierski et al.
# Reference: W. von Smekot-Wensierski et al (1992). Die Absorptionseigenschaften
#           des marinen Phytoplanktons. GKSS Bericht Nr. 92/E/105. 104 pp.
# Modifications by Gege:
#   - Fit with a sum of Lorentzian functions
# - Normalisation to a*(440) = 0.0315 m^2/mg 
# - Rescaled with factor 1.6 according to Gege (2012)

phyto_e <- glue('{wasi_data_dir}/Dinoflagellates.A') %>%
  read_delim(delim = ' ', skip = 11,
             col_types = 'nn',
             col_names = c('nm', 'm2_mg')) %>%
  mutate(phyto_class = "Dinoflagellates")

# Specific absorption of Green algae (m^2/mg)
# 
# Species:   MOUGEOTIA SP.
# Source:    Absorption measurement of Gabi Hartmann in the laboratory (13 October 1993) using a pure culture
# Reference: G. M. Hartmann (1995). Diplomarbeit. Albert-Ludwigs-Universität Freiburg. 91 pp.
# Modifications by Gege:
#   - Fit with a sum of Lorentzian functions
# - Normalisation to a*(440) = 0.0498 m^2/mg  
# - Rescaled with factor 0.7 according to Gege (2012)

phyto_f <- glue('{wasi_data_dir}/GREEN_ALGAE.A') %>%
  read_delim(delim = ' ', skip = 11,
             col_types = 'nn',
             col_names = c('nm', 'm2_mg')) %>%
  mutate(phyto_class = "Green_algae")

astar_chl <- bind_rows(phyto_a, phyto_b, phyto_c, phyto_d, phyto_e, phyto_f) %>%
  tidyr::pivot_wider(names_from = phyto_class, values_from = m2_mg)

astar_chl %>% write_csv('data/gege-wasi-spectra.csv')
