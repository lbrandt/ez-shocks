# German data

#rm(list=ls())

require(tidyverse)
require(readxl)


# Read raw data from datastream request file
datastream = read_excel("C:/Dateien/My Dropbox/Lennart/Thesis/Data/DATASTREAM_REQUEST_MONTHLY.xls", 
                        sheet = "de", col_names = TRUE, na = "NA", skip = 1)
#View(datastream)

datastream_gs = read_excel("C:/Dateien/My Dropbox/Lennart/Thesis/Data/DATASTREAM_REQUEST_MONTHLY.xls", 
                          sheet = "de2", col_names = TRUE, na = "NA", skip = 1)



# Extract dates
dates = datastream[[1]]
ta = first(dates)
te = last(dates)


# Clean data and transform data
data = datastream %>%
  
  select(-starts_with("X__")) %>% # Remove empty columns
  select(-starts_with("Code")) %>% # Remove date columns
  select(-c(BDGDP...D,BDCNPER.D,BDCNGOV.D,BDGCMAC.D,BDGCCON.D,BDGCINT.D,BDEXNGS.D,BDIMNGS.D,BDVAPAAFE,BDVAPAECE,
            BDVAPACND,BDVAPATFD,BDVAPAICD,BDVAPAFID,BDVAPARED,BDVAPASTD,BDVAPAAHD,BDVAPAOSD,BDVAPAICB,BDVAPAFIB))
  
data2 = datastream_gs %>%
  
  select(-starts_with("X__")) %>% # Remove empty columns
  select(-starts_with("Code")) # Remove date columns


# Which series are in both datesets?
var1 = colnames(data)
var2 = colnames(data2)

test1 = matrix(NA,100,105)

for(i in 1:100){
  
  for(j in 1:105){
    
    test1[i,j] = var1[i] == var2[j]
  }
}



tindex = which(test1, arr.ind = TRUE)

ta.index = first(which( rowSums(is.na(data)) == 0 ))
te.index = last(which( rowSums(is.na(data)) == 0 ))




  # Replace R operator % in variable names
  rename(WGUN_TOTQ = "WGUN%TOTQ") %>%
  
  # Aggregate value added
  mutate(bwfinanzverm = (BDVAPAICD*BDVAPAICB + BDVAPAFID*BDVAPAFIB + BDVAPARED*BDVAPAREB + BDVAPASTD*BDVAPASTB)/
                              (BDVAPAICB + BDVAPAFIB + BDVAPAREB + BDVAPASTB)) %>%
  mutate(bwoeffprivdl = (BDVAPAAHD*BDVAPAAHB + BDVAPAOSD*BDVAPAOSB)/ (BDVAPAAHB + BDVAPAOSB)) %>%
  
  # Generate price index investment machinery and equipment
  mutate(piar = BDGCMAC.B/BDGCMAC.D*100)




# Extending time series via growth rates of related variables
source("fun.chain.R")

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




# Remove unnecessary or discontinued series from data
data = data %>%
  
  select(-c(BDVAPAICD, BDVAPAFID, BDVAPARED, BDVAPASTD, BDVAPAICB, BDVAPAFIB, BDVAPAREB, BDVAPASTB)) %>% # bwfinanzverm
  select(-c(BDVAPAAHD, BDVAPAOSD, BDVAPAAHB, BDVAPAOSB)) %>% # bwoeffprivdl
  select(-c(BDGCMAC.B, BDGCMAC.D)) %>% # piar
  
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

