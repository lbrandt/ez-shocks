#rm(list=ls())

require(tidyverse)
require(readxl)

# Set working directory to script file location

# Retrieve script location when code is run within RStudio
location.thisfile = dirname(rstudioapi::getActiveDocumentContext()$path)

# Retrieve script location if file is called via source()
location.thisfile = dirname(sys.frame(1)$ofile)

# Setting working directory to file location
if(getwd() != location.thisfile){setwd(location.thisfile)}else{"yes"}

# Set location of data relative to working directory which contains this script
location.data = normalizePath(file.path("..", "..", "..", "Data"), winslash = "/")




# # Read raw data from datastream request file
# datastream_de = read_excel(file.path(location.data, "DATASTREAM_REQUEST_MONTHLY.xls"), 
#                            sheet = "de", col_names = TRUE, na = "NA", skip = 1)
# 
# datastream_gs = read_excel(file.path(location.data, "DATASTREAM_REQUEST_MONTHLY.xls"), 
#                            sheet = "de2", col_names = TRUE, na = "NA", skip = 1)