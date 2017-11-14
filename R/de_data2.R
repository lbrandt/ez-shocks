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
  
  # Aggregate value added
  mutate(bwfinanzvermiet = (BDVAPAICD*BDVAPAICB + BDVAPAFID*BDVAPAFIB + BDVAPARED*BDVAPAREB + BDVAPASTD*BDVAPASTB)/
                              (BDVAPAICB + BDVAPAFIB + BDVAPAREB + BDVAPASTB)) %>%
  mutate(bwoeffprivdl = (BDVAPAAHD*BDVAPAAHB + BDVAPAOSD*BDVAPAOSB)/ (BDVAPAAHB + BDVAPAOSB)) %>%
  
  # Generate price index investment machinery and equipment
  mutate(piar = BDGCMAC.B/BDGCMAC.D*100)


# Extract survey indicators
surveys = select(data, starts_with("BDIFD"))


# Delete used series from data
data = select(data, -c(BDVAPAICD, BDVAPAFID, BDVAPARED, BDVAPASTD, BDVAPAICB, BDVAPAFIB, BDVAPAREB, BDVAPASTB)) %>%
  
  select(-c(BDVAPAAHD, BDVAPAOSD, BDVAPAAHB, BDVAPAOSB)) %>%
  select(-c(BDGCMAC.B, BDGCMAC.D)) %>%
  select(-starts_with("BDIFD"))


# Check against PW
plot(data$BDCONPRCE)
lines(pw_data_final$X__23)

