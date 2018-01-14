# Grimme-Stoeckli data DE

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




datastream_gs = read_excel(file.path(location.data, "DATASTREAM_REQUEST_MONTHLY.xls"), 
                           sheet = "de2", col_names = TRUE, na = "NA", skip = 1)

# Extract dates
dates = datastream_gs[[1]]
ta = first(dates)
te = last(dates)


# Preliminary cleanup
gs_data = datastream_gs %>%
  
  # Replace R operator % in variable names
  rename(BDUN_TOTQ = "BDUN%TOTQ") %>%
  
  select(-starts_with("X__")) %>% # Remove empty columns
  select(-starts_with("Code")) # Remove date columns


# Extending time series via growth rates of related variables
source("fun.chain.R")

# Chain link labour market data
gs_data$constrempl = fun.chain(gs_data$BDUSMC01B, gs_data$BDUSMB01B, 0, 0)
gs_data$constrwage = fun.chain(gs_data$BDUSMC08B, gs_data$BDUSMB08B, 0, 0)
gs_data$constrhour = fun.chain(gs_data$BDUSMC11B, gs_data$BDUSMB11B, 0, 0)
gs_data$unemplmanu = fun.chain(gs_data$BDEMPMFTP, gs_data$BDEMPMF.P, 0, 0)

# Chain link world trade indices
gs_data$wtradeprod = fun.chain(gs_data$WDCPBPWWG, gs_data$WDCPBPW5G, 0, 0)
gs_data$wtradebala = fun.chain(gs_data$WDCPBTBWG, gs_data$WDCPBTB5G, 0, 0)

# Chain link short term interest rates with corresponding German rates
gs_data$eonia = fun.chain(gs_data$BDSU0304R, gs_data$BDSU0101, 0, 0)
gs_data$is1mo = fun.chain(gs_data$BDSU0310R, gs_data$BDSU0104, 0, 0)
gs_data$is3mo = fun.chain(gs_data$BDSU0316R, gs_data$BDSU0107, 0, 0)


# Remove unnecessary or discontinued series from data
gs_data = gs_data %>%
  
  select(-c(BDTOTEMPP, BDWAGES.F, BDWAGMANF, BDCPCF..F, BDCAR...P)) %>% # Seasonal
  
  select(-c(BDUSMC01B, BDUSMC08B, BDUSMC11B, BDEMPMFTP)) %>% # New labour market data
  select(-c(BDUSMB01B, BDUSMB08B, BDUSMB11B, BDEMPMF.P)) %>% # Discontinued labour market data
  
  select(-c(BDSU0304R, BDSU0310R, BDSU0316R)) %>% # New short term interest rates
  select(-c(BDSU0101, BDSU0104, BDSU0107, BDINTER3)) %>% # Old short term interest rates
  
  select(-c(WDCPBPWWG, WDCPBTBWG)) %>% # New world trade indices
  select(-c(WDCPBPW5G, WDCPBTB5G)) %>% # Discontinued world trade indices
  
  select(-c(BDES7IVPG, BDES77FBG)) %>% # Building permits
  select(-c(BDI..NELF, BDI..RELF)) %>% # Alternative FX computation
  select(-c(BDESPPIEF, BDESPPITF, BDESPPIDF, BDESPPINF, BDCPSERVF)) # Prices




# Build dataset for factor analysis

# Inquire position of the first and last row which contains no empty cells in array
ta.index = first(which( rowSums(is.na(gs_data)) == 0 ))
te.index = last(which( rowSums(is.na(gs_data)) == 0 ))

# Construct time range of largest complete dataset
ta = dates[ta.index]
te = dates[te.index]


# Apply transformations as in GS
dlndata = gs_data %>%
  
  slice(ta.index:te.index) %>%
  
  mutate_at(vars(-c(BDIFDCTNQ, BDIFOBUSQ, BDIFOMTAQ, BDIFORTAQ, BDIFOBDOQ, BDIFOWHAQ, # Ifo Business Conditions
                    BDIFDCTIQ, BDIFDCBIQ, BDIFDCPIQ, BDIFDCCIQ, BDIFDCDIQ, BDIFDCEIQ, # Ifo Construction
                    BDIFDMTFQ, BDIFDMCFQ, BDIFDMPFQ, BDIFDMIFQ, BDIFDMDFQ, BDIFDMNFQ, # Ifo Consumption
                    BDIFDMTCQ, BDIFDMCCQ, BDIFDMPCQ, BDIFDMICQ, BDIFDMDCQ, BDIFDMNCQ, 
                    BDIFDFBCQ, BDIFRS.CQ, BDIFWS.CQ, 
                    BDEUSCMPQ, BDEUSCSAQ, BDEUSCFHQ, # DG ECFIN
                    
                    BDUN_TOTQ, BDESUNEMO, BDESUNUPQ, # UNP rates
                    
                    BDUUCG05P, # Short time workers
                    
                    eonia, is1mo, is3mo, BDWU0913, BDWU0915, BDWU8612, # Interest rates and yields
                    USTRCN3., USTRCN5., USTRCN10, BDWU0022, BDWU0004R
  )), funs(log)) %>%
  
  mutate_at(vars(-c(BDIFDCTNQ, BDIFOBUSQ, BDIFOMTAQ, BDIFORTAQ, BDIFOBDOQ, BDIFOWHAQ, # Ifo Business Conditions
                    BDIFDMCCQ, BDIFDMPCQ, BDIFDMICQ, BDIFDMDCQ, BDIFDMNCQ, BDIFDFBCQ, BDIFRS.CQ, BDIFWS.CQ, # Ifo Inventories
                    
                    BDUUCG05P, # Short time workers
                    
                    USTRCN3., USTRCN5., USTRCN10, BDWU0022 # Bond Yield Spreads
  )), funs(. - lag(.)))



# Save data to workspace
save(dates, gs_data, dlndata, file = "gs_data.RData")


# Prepare data for export
out.varnames  = colnames(dlndata)
out.dates     = as.character.Date(dates[ta.index:te.index])
out.data      = as.matrix(gs_data[ta.index:te.index, ]) # level data
out.dlndata   = as.matrix(dlndata[(ta.index+1):te.index, ]) # diffed series are one observation shorter


# Save results to HDF5
out.file = h5file("gs_data.h5", mode = "w")
out.file["/varnames"] = out.varnames
out.file["/dates"] = out.dates
out.file["/data"] = out.data
out.file["/dlndata"] = out.dlndata
h5close(out.file)
