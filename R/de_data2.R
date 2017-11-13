# German data

#rm(list=ls())

require(tidyverse)

# Read raw data from datastream request file
datastream = read_excel("C:/Dateien/My Dropbox/Lennart/Thesis/Data/DATASTREAM_REQUEST_QUARTERLY.xls", 
                        sheet = "de", col_names = TRUE, na = "NA")
#View(datastream)


# Extract series codes and descriptions
names = filter(datastream, )
