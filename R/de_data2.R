# German data

#rm(list=ls())

require(tidyverse)
require(readxl)


# Read Pirschel-Wolters data for comparison
pw_data_final <- read_excel("C:/Dateien/My Dropbox/Lennart/Thesis/Data/pw_data_final.xls", 
                            sheet = "datafinal", range = "C90:EL189", col_names = FALSE, na = "NA")


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
  
  select(-starts_with("X__")) %>% # Delete separator columns
  select(-starts_with("Code")) %>% # Delete date columns
  
  rename(WGUN_TOTQ = "WGUN%TOTQ") %>% # Replace R operator % in variable names
  
  # Aggregate value added
  mutate(bwfinanzverm = (BDVAPAICD*BDVAPAICB + BDVAPAFID*BDVAPAFIB + BDVAPARED*BDVAPAREB + BDVAPASTD*BDVAPASTB)/
                              (BDVAPAICB + BDVAPAFIB + BDVAPAREB + BDVAPASTB)) %>%
  mutate(bwoeffprivdl = (BDVAPAAHD*BDVAPAAHB + BDVAPAOSD*BDVAPAOSB)/ (BDVAPAAHB + BDVAPAOSB)) %>%
  
  # Generate price index investment machinery and equipment
  mutate(piar = BDGCMAC.B/BDGCMAC.D*100)




# Extending time series via growth rates of related variables

# Chain link UNP with UNP in Western Germany
data$unrate = fun.chain(data$BDUSCC02Q, data$WGUN_TOTQ, 0, 0)

# Chain link short term interest rates with corresponding German rates
data$eonia = fun.chain(data$BDSU0304R, data$BDSU0101, 0, 0)
data$is1mo = fun.chain(data$BDSU0310R, data$BDSU0104, 0, 0)
data$is3mo = fun.chain(data$BDSU0316R, data$BDSU0107, 0, 0)








# Extract survey indicators
surveys = select(data, starts_with("BDIFD"))


# Remove unnecessary or discontinued series from data
data = data %>%
  
  select(-c(BDVAPAICD, BDVAPAFID, BDVAPARED, BDVAPASTD, BDVAPAICB, BDVAPAFIB, BDVAPAREB, BDVAPASTB)) %>% # bwfinanzverm
  select(-c(BDVAPAAHD, BDVAPAOSD, BDVAPAAHB, BDVAPAOSB)) %>% # bwoeffprivdl
  select(-c(BDGCMAC.B, BDGCMAC.D)) %>% # piar
  
  select(-c(BDEA4001B, BDEA4100B, BDEA4170B, BDEA4220B)) %>% # Current account
  
  select(-c(BDUSCC02Q, WGUN_TOTQ)) %>% # Unemployment rate
  select(-c(BDSU0304R, BDSU0310R, BDSU0316R)) %>% # New short term interest rates
  select(-c(BDSU0101, BDSU0104, BDSU0107)) %>% # Old short term interest rates
  
  select(-c(BDUSC.04O, BDUSCC04O)) %>% # Vacancies are in series BDQLM004O

  
  select(-starts_with("BDIFD")) # Survey indicators



de.T = dim(data)[1]
de.N = dim(data)[2] 


# Check against PW



data$is1mo = fun.chain(data$BDSU0310R, data$BDSU0104, 0, 0)

eonia
is1mo

plot(0, type="n", xlab="Time", ylab="Series", xlim = c(1,de.T), ylim = c(0, 15))


lines(data$BDSU0310R, col = "blue")
lines(data$BDSU0104)
lines(data$is1mo)
lines(pw_data_final$X__45)



