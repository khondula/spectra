library(hdf5r)
library(tidyverse)
library(glue)

# exploring h5 file structure

# open a file for reading
my_h5 <- hdf5r::H5File$new(my_file, mode = "r")
class(my_h5)

# get all names of objects in a group or in root directory 
names(my_h5)

hdf5r::list.objects(my_h5) # everything? 59 things
hdf5r::list.datasets(my_h5) # 47 datasets
hdf5r::list.groups(my_h5) # 12 groups

# "BLAN/Reflectance/Metadata/Logs/185613/Solar_Azimuth_Angle"
# "BLAN/Reflectance/Metadata/Logs/184816/Solar_Zenith_Angle" 
# "BLAN/Reflectance/Metadata/Spectral_Data"                                         
# "BLAN/Reflectance/Metadata/Spectral_Data/FWHM"                                    
# "BLAN/Reflectance/Metadata/Spectral_Data/Wavelength"                              
# "BLAN/Reflectance/Metadata/to-sensor_azimuth_angle"                               
# "BLAN/Reflectance/Metadata/to-sensor_zenith_angle"                                
# "BLAN/Reflectance/Reflectance_Data" 


