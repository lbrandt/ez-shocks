# Georgiadis Jancokova Monetary Policy Shocks Database

#rm(list=ls())

require(tidyverse)
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

gj_import_q = read_dta(file.path(location.data, "Georgiadis Jancokova Shocks Database", "mpshocks_database_public.dta"))

