# Babecka Kucharcukova Replication data

#rm(list=ls())

require(tidyverse)
require(readxl)
require(haven)
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




datastream_bab = read_excel(file.path(location.data, "DATASTREAM_REQUEST_MONTHLY.xls"), 
                            sheet = "bab", col_names = TRUE, na = "NA", skip = 1, guess_max = 6000)

# Extract dates
source("monthStart.R")
dates = monthStart(datastream_bab[[1]])
#ta = first(dates)
#te = last(dates)

# Build dataset
bab_data = datastream_bab %>%
  
  select(-starts_with("X__")) %>% # Remove empty columns
  select(-starts_with("Code")) %>% # Remove date columns
  
  mutate_at(vars(EMECASM), funs(replace(., is.na(.), 0))) # Set NA values in EMECASM to zero
  
  mutate_at(vars(-c(EMPRATE., EMIBOR3., EMIBOR1Y, EMGBOND.,
                  OIEURSW, OIEUR2W, OIEUR1M, OIEUR3M, OIEUR10, OIEUR1Y)), # Exclude interest rates
          funs(. - lag(., n = 12))) # Year on year differences
  

plot(bab_data$EUDOLLR, type = "l")
lines(1/bab_data$USEURSP, col = "blue")



# Inquire position of the first and last row which contains no empty cells in array
ta.index = first(which( rowSums(is.na(bab_data)) == 0 ))
te.index = last(which( rowSums(is.na(bab_data)) == 0 ))

# Construct time range of largest complete dataset
ta = dates[ta.index]
te = dates[te.index]
sample = ta.index:te.index


# Save data to workspace
save(dates, bab_data, file = "bab_data.RData")

# Prepare data for export
out.varnames  = colnames(bab_data)
out.dates     = as.character.Date(dates)
out.data      = as.matrix(bab_data[sample, ]) # level data

# Save results to HDF5
out.file = h5file("bab_data.h5", mode = "w")
out.file["/varnames"] = out.varnames
out.file["/dates"] = out.dates
out.file["/data"] = out.data
h5close(out.file)
