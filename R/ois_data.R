# EUR OIS daily data

#rm(list=ls())

require(tidyverse)
require(readxl)
require(h5)

# Set working directory to script file location

# Retrieve script location when code is run within RStudio
location.thisfile = dirname(rstudioapi::getActiveDocumentContext()$path)
# Retrieve script location if file is called via source()
location.thisfile = dirname(sys.frame(1)$ofile)
# Setting working directory to file location
setwd(location.thisfile)

# Set location of data relative to working directory which contains this script
location.data = normalizePath(file.path("..", "..", "..", "Data"), winslash = "/")




datastream_ois = read_excel(file.path(location.data, "DATASTREAM_REQUEST_DAILY.xls"), 
                            sheet = "ois", col_names = TRUE, na = "NA", skip = 1, guess_max = 6000)

# Extract dates
dates = datastream_ois[[1]]
#ta = first(dates)
#te = last(dates)

# Restrict sample
ta = as.POSIXct("1999-01-01", tz = "UTC")
te = as.POSIXct("2018-01-01", tz = "UTC")

ta.index = which(dates == ta)
te.index = which(dates == te)

sample = ta.index:te.index

# Build dataset
dates = dates[sample]

ois_data = datastream_ois %>%
  
  rename(EUEONIA = "EUEONIA(IO)") %>%
  
  select(-starts_with("X__")) %>% # Remove empty columns
  select(-starts_with("Code")) %>% # Remove date columns
  
  slice(sample)


# Save data to workspace
save(dates, ois_data, file = "ois_data.RData")

# Prepare data for export
out.varnames  = colnames(ois_data)
out.dates     = as.character.Date(dates)
out.data      = as.matrix(ois_data) # level data

# Save results to HDF5
out.file = h5file("ois_data.h5", mode = "w")
out.file["/varnames"] = out.varnames
out.file["/dates"] = out.dates
out.file["/data"] = out.data
h5close(out.file)
