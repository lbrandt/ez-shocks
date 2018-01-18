# EZ data

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




# Read Fagan AWM data (quarterly)
awmdata_ez = read_csv(file.path(location.data, "ez_awmdata.csv"))

# Read raw data from datastream request file
#datastream_de = read_excel(file.path(location.data, "DATASTREAM_REQUEST_MONTHLY.xls"), 
#                           sheet = "de", col_names = TRUE, na = "NA", skip = 1)


# Extract dates
dates = awmdata_ez[[1]]
ta = first(dates)
te = last(dates)

# Inquire position of the first and last row which contains no empty cells in array
ta.index = first(which( rowSums(is.na(awmdata_ez)) == 0 ))
te.index = last(which( rowSums(is.na(awmdata_ez)) == 0 ))

# Preliminary cleanup
ez_data = awmdata_ez %>%

  select(-starts_with("X1")) # Remove data column

# Restrict sample to 1990Q1
index90 = which(dates == "1990Q1")
ez_data90 = ez_data[index90:te.index, ]
dates = dates[index90:te.index]

# Delete short series
indexNA = which( colSums(is.na(ez_data90[])) != 0 )
nameNA = colnames(ez_data90)[indexNA]
ez_data90 = select(ez_data90, -one_of(nameNA))

# Reset ta te
ta.index = first(which( rowSums(is.na(ez_data90)) == 0 ))
te.index = last(which( rowSums(is.na(ez_data90)) == 0 ))

ta = dates[ta.index]
te = dates[te.index]


# Do not log series which contain negative numbers
indexNEG = which( colSums((ez_data90<=0)) != 0)
nameNEG = colnames(ez_data90)[indexNEG]

# Apply transformations for stationarity
dlndata = ez_data90 %>%
  
  mutate_at(vars(-one_of(nameNEG)),
            funs(log)) %>%
  
  mutate_at(vars(), # Ifo Inventories
            funs(. - lag(.))) %>%
  
  slice(2:te.index) # delete first row b/c diff

# Save data to workspace
save(dates, ez_data90, dlndata, file = "ez_data90_q.RData")


# Prepare data for export
out.varnames  = colnames(dlndata)
out.dates     = as.character.Date(dates[ta.index:te.index])
out.data      = as.matrix(ez_data90[ta.index:te.index, ]) # level data
out.dlndata   = as.matrix(dlndata)[,] # diffed series are one observation shorter


# Save results to HDF5
out.file = h5file("ez_data90_q.h5", mode = "w")
out.file["/varnames"] = out.varnames
out.file["/dates"] = out.dates
out.file["/data"] = out.data
out.file["/dlndata"] = out.dlndata
h5close(out.file)
