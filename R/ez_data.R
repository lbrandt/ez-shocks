# -------------------------------------------------------------------------
# Loads time series on Eurozone economy from Datastream request file and
# transforms them according to accompanying meta data file. Saves results
# to HDF5 container for further use.
#
# Lennart Brandt, Feb 2018
# -------------------------------------------------------------------------

# Initialize CRAN packages

#install.packages("tidyverse")
#install.packages("readxl")
#install.packages("h5")
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




# EZ data monthly --------------------------------------------------------
datastream_ez = read_excel(file.path(location.data, "DATASTREAM_REQUEST_MONTHLY.xls"), 
                           sheet = "ez", col_names = TRUE, na = "NA", skip = 1, guess_max = 6000)


# Extract dates
source("dateShift.R")
dates = dateShift(datastream_ez[[1]], unit = 'month', rule = 'start')


# Build dataset
ez_data = datastream_ez %>%
  
  select(-starts_with("X__")) %>% # Remove empty columns
  select(-starts_with("Code")) # Remove date columns

colnames(ez_data) = gsub("[&, %, $]", "", colnames(ez_data)) # Remove R operators from names

# Read metadata
ez_meta = read_excel(file.path(location.data, "ez_data_meta.xlsx"), col_names = TRUE, na = "NA")
ez_meta$code = gsub("[&, %, $]", "", ez_meta$code) # Remove R operators from names

# Write transform column to ez_data attributes
for(i in 1:dim(ez_data)[2]){
  attr(ez_data, "transformation")[i] = ez_meta$transform[which(ez_meta$code == names(ez_data)[i])]
}

ez_tdata = ez_data %>%
  mutate_if(.predicate = attributes(ez_data)$transformation == 1, funs(. - lag(.))) %>% # Diffs
  mutate_if(.predicate = attributes(ez_data)$transformation == 2, funs(log(.) - lag(log(.)))) # Logdiffs




# Inquire position of the first and last row which contains no empty cells in array
ta.index = first(which( rowSums(is.na(ez_data)) == 0 ))
te.index = last(which( rowSums(is.na(ez_data)) == 0 ))

# Construct time range of largest complete dataset
ta = dates[ta.index]
te = dates[te.index]
sample = ta.index:te.index

# Save data to workspace
save(dates, ez_data, ez_tdata, ez_meta, file = "ez_data.RData")

# Prepare data for export
out.names = colnames(ez_data)
out.dates = as.character.Date(dates[sample])
out.data  = as.matrix(ez_data[sample, ])
out.tdata = as.matrix(ez_tdata[sample, ])

# Save results to HDF5
out.file = h5file("ez_data.h5", mode = "w")
out.file["/names"] = out.names
out.file["/dates"] = out.dates
out.file["/data"]  = out.data
out.file["/tdata"] = out.tdata
h5close(out.file)
