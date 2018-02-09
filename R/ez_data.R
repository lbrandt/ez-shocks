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




# EZ data monthly --------------------------------------------------------
datastream_ez = read_excel(file.path(location.data, "DATASTREAM_REQUEST_MONTHLY.xls"), 
                           sheet = "ez", col_names = TRUE, na = "NA", skip = 1, guess_max = 6000)


# Extract dates
source("dateShift.R")
dates = dateShift(datastream_ez[[1]], unit = 'month', rule = 'start')
#ta = first(dates)
#te = last(dates)


# Build dataset
ez_data = datastream_ez %>%
  
  rename(EKESCDER = "EKESCDE&R") %>%
  rename(USUNTOTQ = "USUN%TOTQ") %>%
  rename(SWHWWAEF = "SWHWWAE$F") %>%
  rename(SPCOMP = "S&PCOMP") %>%
  
  select(-starts_with("X__")) %>% # Remove empty columns
  select(-starts_with("Code")) %>% # Remove date columns
  
  select(-c(EKESIENGG__1)) # Remove duplicate series


# Read metadata
ez_meta = read_csv(file.path(location.data, "ez_data_varlist.csv"), col_names = TRUE, col_types = "icccciiii", na = "NA")
  

ez_data = ez_data %>%
  
  mutate_at(vars(-c(EMPRATE., EMIBOR3., EMIBOR1Y, EMGBOND., OIEURSW, OIEUR2W, OIEUR1M, OIEUR3M, OIEUR10, OIEUR1Y, USEURSP, # Interest rates
                    EKEBUN..O, EMUNPTOTO, EKEBUN..O, EMTOTUNQ, EKUNTOTQ, # Unemployment rates
                    EMEBCPGS, EKCAFBALA)), # CPI rate & Balance 
            funs(log(.))) %>% # Log
  
  #mutate_at(vars(c(EMM1....B, EMM2....B, EMM3....B, EMASTOT, EMECASM, EMECLBC, EMECLEM,
  #                 EKEBCARRO, EMACECARP, EMECBALEA, EMSHRPRCF, EMCRDCONA)), funs(log(.))) %>% # Log
  
  mutate_at(vars(EMECASM), funs(replace(., is.na(.), 0))) # Set NA values in EMECASM to zero


# Inquire position of the first and last row which contains no empty cells in array
ta.index = first(which( rowSums(is.na(ez_data)) == 0 ))
te.index = last(which( rowSums(is.na(ez_data)) == 0 ))

# Construct time range of largest complete dataset
ta = dates[ta.index]
te = dates[te.index]
sample = ta.index:te.index

# Save data to workspace
save(dates, ez_vardata, file = "ez_data.RData")

# Prepare data for export
out.varnames  = colnames(ez_data)
out.dates     = as.character.Date(dates[sample])
out.data      = as.matrix(ez_data[sample, ])

# Save results to HDF5
out.file = h5file("ez_data.h5", mode = "w")
out.file["/varnames"] = out.varnames
out.file["/dates"] = out.dates
out.file["/data"] = out.data
h5close(out.file)






# VAR data monthly --------------------------------------------------------
datastream_var = read_excel(file.path(location.data, "DATASTREAM_REQUEST_MONTHLY.xls"), 
                            sheet = "var", col_names = TRUE, na = "NA", skip = 1, guess_max = 6000)


# Extract dates
source("dateShift.R")
dates = dateShift(datastream_var[[1]], unit = 'month', rule = 'start')
#ta = first(dates)
#te = last(dates)


# Build dataset
ez_vardata = datastream_var %>%
  
  rename(EMEBCPGS = "EMEBCPGS%") %>%
  rename(EMTOTUNQ = "EMTOTUN%Q") %>%
  rename(EKUNTOTQ = "EKUN%TOTQ") %>%
  
  select(-starts_with("X__")) %>% # Remove empty columns
  select(-starts_with("Code")) %>% # Remove date columns
  select(-c(EKIMPPRCF)) %>% # Remove short series
  
  mutate_at(vars(-c(EMPRATE., EMIBOR3., EMIBOR1Y, EMGBOND., OIEURSW, OIEUR2W, OIEUR1M, OIEUR3M, OIEUR10, OIEUR1Y, USEURSP, # Interest rates
                    EKEBUN..O, EMUNPTOTO, EKEBUN..O, EMTOTUNQ, EKUNTOTQ, # Unemployment rates
                    EMEBCPGS, EKCAFBALA)), # CPI rate & Balance 
            funs(log(.))) %>% # Log

  #mutate_at(vars(c(EMM1....B, EMM2....B, EMM3....B, EMASTOT, EMECASM, EMECLBC, EMECLEM,
  #                 EKEBCARRO, EMACECARP, EMECBALEA, EMSHRPRCF, EMCRDCONA)), funs(log(.))) %>% # Log

  mutate_at(vars(EMECASM), funs(replace(., is.na(.), 0))) # Set NA values in EMECASM to zero


# Inquire position of the first and last row which contains no empty cells in array
ta.index = first(which( rowSums(is.na(ez_vardata)) == 0 ))
te.index = last(which( rowSums(is.na(ez_vardata)) == 0 ))

# Construct time range of largest complete dataset
ta = dates[ta.index]
te = dates[te.index]
sample = ta.index:te.index

# Save data to workspace
save(dates, ez_vardata, file = "ez_vardata.RData")

# Prepare data for export
out.varnames  = colnames(ez_vardata)
out.dates     = as.character.Date(dates[sample])
out.data      = as.matrix(ez_vardata[sample, ])

# Save results to HDF5
out.file = h5file("ez_vardata.h5", mode = "w")
out.file["/varnames"] = out.varnames
out.file["/dates"] = out.dates
out.file["/data"] = out.data
h5close(out.file)




# AWM data quarterly ------------------------------------------------------
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
