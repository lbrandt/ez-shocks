# Test source file location and relative paths
#rm(list=ls())

# Retrieve script location if file is called via source()
location.fromsource = dirname(sys.frame(1)$ofile)

# Retrieve script location when code is run within RStudio
location.fromscript = dirname(rstudioapi::getActiveDocumentContext()$path)


#setwd(location.fromsource)
