# German data

#rm(list=ls())

require(tidyverse)
require(readxl)

# Set working directory to script file location

# Retrieve script location when code is run within RStudio
location.thisfile = dirname(rstudioapi::getActiveDocumentContext()$path)

# Retrieve script location if file is called via source()
location.thisfile = dirname(sys.frame(1)$ofile)

# Setting working directory to file location
setwd(location.thisfile)

# Set location of data relative to working directory which contains this script
location.data = normalizePath(file.path("..", "..", "..", "Data"), winslash = "/")




# Read raw data from datastream request file
datastream_de = read_excel(file.path(location.data, "DATASTREAM_REQUEST_MONTHLY.xls"), 
                           sheet = "de", col_names = TRUE, na = "NA", skip = 1)

datastream_gs = read_excel(file.path(location.data, "DATASTREAM_REQUEST_MONTHLY.xls"), 
                           sheet = "de2", col_names = TRUE, na = "NA", skip = 1)


# datastream_de = read_excel("C:/Dateien/My Dropbox/Lennart/Thesis/Data/DATASTREAM_REQUEST_MONTHLY.xls", 
#                           sheet = "de", col_names = TRUE, na = "NA", skip = 1)
# 
# datastream_gs = read_excel("C:/Dateien/My Dropbox/Lennart/Thesis/Data/DATASTREAM_REQUEST_MONTHLY.xls", 
#                           sheet = "de2", col_names = TRUE, na = "NA", skip = 1)



# Extract dates
dates = datastream_de[[1]]
ta = first(dates)
te = last(dates)


# Clean data and transform data
data = datastream_de %>%
  
  # Replace R operator % in variable names
  rename(WGUN_TOTQ = "WGUN%TOTQ") %>%
  
  select(-starts_with("X__")) %>% # Remove empty columns
  select(-starts_with("Code")) # Remove date columns

  
data2 = datastream_gs %>%
  
  # Replace R operator % in variable names
  rename(BDUN_TOTQ = "BDUN%TOTQ") %>%
  
  select(-starts_with("X__")) %>% # Remove empty columns
  select(-starts_with("Code")) # Remove date columns


# Which series are in both datesets?
var1 = colnames(data)
var2 = colnames(data2)

test2 = cbind(c(var1, var2))

test1 = matrix(NA,100,105)

for(i in 1:100){
  
  for(j in 1:105){
    
    test1[i,j] = var1[i] == var2[j]
  }
}

sum(test1)

tindex = which(test1, arr.ind = TRUE)

# There are 8 series which are in both data sets. Merging should result in a 323x197 df
data3 = merge(data, data2)



ta.index = first(which( rowSums(is.na(data)) == 0 ))
te.index = last(which( rowSums(is.na(data)) == 0 ))




  





# Extending time series via growth rates of related variables
source("fun.chain.R")

# Chain link labour market data
data2$constrempl = fun.chain(data2$BDUSMC01B, data2$BDUSMB01B, 0, 0)
data2$constrwage = fun.chain(data2$BDUSMC08B, data2$BDUSMB08B, 0, 0)
data2$constrhour = fun.chain(data2$BDUSMC11B, data2$BDUSMB11B, 0, 0)
data2$unemplmanu = fun.chain(data2$BDEMPMFTP, data2$BDEMPMF.P, 0, 0)

# Chain link world trade indices
data2$wtradeprod = fun.chain(data2$WDCPBPWWG, data2$WDCPBPW5G, 0, 0)
data2$wtradebala = fun.chain(data2$WDCPBTBWG, data2$WDCPBTB5G, 0, 0)

# Chain link short term interest rates with corresponding German rates
data2$eonia = fun.chain(data2$BDSU0304R, data2$BDSU0101, 0, 0)
data2$is1mo = fun.chain(data2$BDSU0310R, data2$BDSU0104, 0, 0)
data2$is3mo = fun.chain(data2$BDSU0316R, data2$BDSU0107, 0, 0)


View(colnames(data2))

# Remove unnecessary or discontinued series from data
data2 = data2 %>%
  
  select(-c(BDUSMC01B, BDUSMC08B, BDUSMC11B, BDEMPMFTP)) %>% # New labour market data
  select(-c(BDUSMB01B, BDUSMB08B, BDUSMB11B, BDEMPMF.P)) %>% # Discontinued labour market data
  
  select(-c(BDSU0304R, BDSU0310R, BDSU0316R)) %>% # New short term interest rates
  select(-c(BDSU0101, BDSU0104, BDSU0107, BDINTER3)) %>% # Old short term interest rates
  
  select(-c(WDCPBPWWG, WDCPBTBWG)) %>% # New world trade indices
  select(-c(WDCPBPW5G, WDCPBTB5G)) %>% # Discontinued world trade indices
  
  
  select(-c(BDES7IVPG, BDES77FBG)) %>% # Building permits
  select(-c(BDI..NELF, BDI..RELF)) %>% # Alternative FX computation
  select(-c(BDESPPIEF, BDESPPITF, BDESPPIDF, BDESPPINF, BDCPSERVF)) # Prices



#### Add more finance from pw data


# Build dataset for factor analysis 

# Apply natural log to variables excluding interest rates and survey data
lndata2 = mutate_at(data2, vars(-c(BDIFDCTNQ, BDIFOBUSQ, BDIFOMTAQ, BDIFORTAQ, BDIFOBDOQ, BDIFOWHAQ, # Ifo Business Conditions
                                   BDIFDCTIQ, BDIFDCBIQ, BDIFDCPIQ, BDIFDCCIQ, BDIFDCDIQ, BDIFDCEIQ, # Ifo Construction
                                   BDIFDMTFQ, BDIFDMCFQ, BDIFDMPFQ, BDIFDMIFQ, BDIFDMDFQ, BDIFDMNFQ, # Ifo Consumption
                                   BDIFDMTCQ, BDIFDMCCQ, BDIFDMPCQ, BDIFDMICQ, BDIFDMDCQ, BDIFDMNCQ, 
                                   BDIFDFBCQ, BDIFRS.CQ, BDIFWS.CQ, 
                                   BDEUSCMPQ, BDEUSCSAQ, BDEUSCFHQ, # DG ECFIN
                                  
                                   BDUN_TOTQ, BDESUNEMO, BDESUNUPQ, # UNP rates
                                  
                                   BDUUCG05P, # Short time workers
                                  
                                   eonia, is1mo, is3mo, BDWU0913, BDWU0915, BDWU8612, # Interest rates and yields
                                   USTRCN3., USTRCN5., USTRCN10, BDWU0022, BDWU0004R
                                   )), funs(log))


# Diff for stationarity
dlndata2  = data.frame(sapply(lndata2, FUN = diff))


# # Assign descriptive column names
# varlist = read.csv("de_varlist.csv")
# varnames = varlist[2]
# 
# colnames(dlndata) = varnames[[1]]


# Inquire position of the first and last row which contains no empty cells in array
ta.index = first(which( rowSums(is.na(data2)) == 0 ))
te.index = last(which( rowSums(is.na(data2)) == 0 ))

# Construct time range of largest complete dataset
ta = dates[ta.index]
te = dates[te.index]



# Save data to workspace
save(dates, data2, lndata2, dlndata2, file = "de_gsdata.RData")


# # Save data and varnames to csv
# out.dsrequest = colnames(datastream)
# out.varcodes  = colnames(data)
# out.varnames  = colnames(dlndata)
# out.dates     = dates[ta.index:te.index]
# 
# out.data      = data[ta.index:te.index, ] # level data
out.dlndata2   = dlndata2[ta.index:te.index-1, ] # diffed series are one observation shorter
# 
# out.surveys   = surveys[ta.index:te.index, ] # surveys still contain NA elements!
# out.dsurveys  = dsurveys[ta.index:te.index-1, ]
# 
# # Write .csv
# write.table(out.dates, file = 'de_dates.csv', row.names = FALSE, col.names = FALSE, sep = ',')
# write.table(out.dsrequest, file = 'de_request.csv', row.names = FALSE, col.names = FALSE, sep = ',')
# write.table(out.varnames, file = 'de_varnames.csv', row.names = FALSE, col.names = FALSE, sep = ',')
# 
# write.table(out.data, file = 'de_data.csv', row.names = FALSE, col.names = FALSE, sep = ',')
 write.table(out.dlndata2, file = 'de_gsdata.csv', row.names = FALSE, col.names = FALSE, sep = ',')
# 
# write.table(out.surveys, file = 'de_surveys.csv', row.names = FALSE, col.names = FALSE, sep = ',')
# write.table(out.dsurveys, file = 'de_surveys2.csv', row.names = FALSE, col.names = FALSE, sep = ',')




#####################################################################################################
#### OLD

# Chain link UNP with UNP in Western Germany
data$unrate = fun.chain(data$BDUSCC02Q, data$WGUN_TOTQ, 0, 0)

# Chain link short term interest rates with corresponding German rates
data$eonia = fun.chain(data$BDSU0304R, data$BDSU0101, 0, 0)
data$is1mo = fun.chain(data$BDSU0310R, data$BDSU0104, 0, 0)
data$is3mo = fun.chain(data$BDSU0316R, data$BDSU0107, 0, 0)

# Chain link Turnover in construction
data$toconstind = fun.chain(data$BDUSMC31B, data$BDUSMB31B, 0, 0)
data$toconstres = fun.chain(data$BDUSMC36B, data$BDUSMB36B, 0, 0)
data$toconstpub = fun.chain(data$BDUSMC37B, data$BDUSMB37B, 0, 0)




# Extract survey indicators
surveys = select(data, starts_with("BDIFD"))


View(colnames(data))

# Remove unnecessary or discontinued series from data
data = data %>%
  
  select(-c(BDIPMAN.G__1)) %>% # Duplicate series
  select(-c(BDGOVBALA, BDBQ9059A, BDGGOVBLA)) %>% # Government balance
  
  select(-c(BDUSCC02Q, WGUN_TOTQ)) %>% # Unemployment rate
  select(-c(BDSU0304R, BDSU0310R, BDSU0316R)) %>% # New short term interest rates
  select(-c(BDSU0101, BDSU0104, BDSU0107)) %>% # Old short term interest rates
  select(-c(BDUSMC36B, BDUSMC31B, BDUSMC37B)) %>% # New turnover in construction
  select(-c(BDUSMB36B, BDUSMB31B, BDUSMB37B)) %>% # Old turnover in construction
  
  select(-c(BDUSC.04O, BDUSCC04O)) %>% # Vacancies in BDQLM004O
  select(-c(BDEA4001B, BDEA4100B, BDEA4170B, BDEA4220B)) %>% # Current account in BDBGDBAQB, BDBSVBAQB, BDBI1BAQB, BDBI2BAQB
  select(-c(WGUS01NAG, WGUS02NAG, BDUSNA01G)) %>% # Unused IP
  select(-c(BDWAGES.F)) %>% # Wages

  select(-starts_with("BDIFD")) # Survey indicators




# Build dataset for factor analysis 

# Apply natural log to variables excluding interest rates and current account series
lndata  = mutate_at(data, vars(-c(BDWU0898, BDWU0899, BDWU0900, BDWU0901, BDWU0902, BDWU0903, BDWU8606, BDWU8607, BDWU8608, # Interest rates
                                 eonia, is1mo, is3mo, # Interest rates
                                 BDBGDBAQB, BDBSVBAQB, BDBI1BAQB, BDBI2BAQB)), funs(log)) # Current account

# Diff for stationarity
dlndata  = data.frame(sapply(lndata, FUN = diff))
dsurveys = data.frame(sapply(surveys, FUN = diff))

# Assign descriptive column names
varlist = read.csv("de_varlist.csv")
varnames = varlist[2]

colnames(dlndata) = varnames[[1]]


# Inquire position of the first and last row which contains no empty cells in array
ta.index = first(which( rowSums(is.na(data)) == 0 ))
te.index = last(which( rowSums(is.na(data)) == 0 ))

# Construct time range of largest complete dataset
ta = dates[ta.index]
te = dates[te.index]



# Save data to workspace
save(dates, data, lndata, dlndata, surveys, dsurveys, file = "de_data2.RData")


# Save data and varnames to csv
out.dsrequest = colnames(datastream)
out.varcodes  = colnames(data)
out.varnames  = colnames(dlndata)
out.dates     = dates[ta.index:te.index]

out.data      = data[ta.index:te.index, ] # level data
out.dlndata   = dlndata[ta.index:te.index-1, ] # diffed series are one observation shorter

out.surveys   = surveys[ta.index:te.index, ] # surveys still contain NA elements!
out.dsurveys  = dsurveys[ta.index:te.index-1, ]

# Write .csv
write.table(out.dates, file = 'de_dates.csv', row.names = FALSE, col.names = FALSE, sep = ',')
write.table(out.dsrequest, file = 'de_request.csv', row.names = FALSE, col.names = FALSE, sep = ',')
write.table(out.varnames, file = 'de_varnames.csv', row.names = FALSE, col.names = FALSE, sep = ',')

write.table(out.data, file = 'de_data.csv', row.names = FALSE, col.names = FALSE, sep = ',')
write.table(out.dlndata, file = 'de_data2.csv', row.names = FALSE, col.names = FALSE, sep = ',')

write.table(out.surveys, file = 'de_surveys.csv', row.names = FALSE, col.names = FALSE, sep = ',')
write.table(out.dsurveys, file = 'de_surveys2.csv', row.names = FALSE, col.names = FALSE, sep = ',')

#####################################################################################################
