# German data

#rm(list=ls())

require(tidyverse)
require(readxl)


# Read raw data from datastream request file
datastream = read_excel("C:/Dateien/My Dropbox/Lennart/Thesis/Data/DATASTREAM_REQUEST_QUARTERLY.xls", 
                        sheet = "de", col_names = TRUE, na = "NA", skip = 1)
#View(datastream)



# Extract dates
dates = datastream[[1]]
ta = first(dates)
te = last(dates)


# Clean data and transform data
data = datastream %>%
  
  select(-starts_with("X__")) %>% # Remove separator columns
  select(-starts_with("Code")) %>% # Remove date columns
  
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




# Transform data into stationary series

# Apply natural log to variables excluding interest rates and current account series
lndata  = mutate_at(data, vars(-c(BDWU0898, BDWU0899, BDWU0900, BDWU0901, BDWU0902, BDWU0903, BDWU8606, BDWU8607, BDWU8608, # Interest rates
                                 eonia, is1mo, is3mo, # Interest rates
                                 BDBGDBAQB, BDBSVBAQB, BDBI1BAQB, BDBI2BAQB)), funs(log)) # Current account

# Diff for stationarity
dlndata  = data.frame(sapply(lndata, FUN = diff))
dsurveys = data.frame(sapply(surveys, FUN = diff))




# Save data to workspace and .csv
save(dates, data, lndata, dlndata, surveys, dsurveys, file = "de_data2.RData")

out.data     = add_column(data, dates = dates, .before = 1)
out.dlndata  = add_column(dlndata, dates = dates[2:108], .before = 1)

out.surveys  = add_column(surveys, dates = dates, .before = 1)
out.dsurveys = add_column(dsurveys, dates = dates[2:108], .before = 1)


write.csv(out.data, file = 'de_data.csv', row.names = FALSE)
write.csv(out.dlndata, file = 'de_data2.csv', row.names = FALSE)

write.csv(out.surveys, file = 'de_surveys.csv', row.names = FALSE)
write.csv(out.dsurveys, file = 'de_surveys2.csv', row.names = FALSE)

