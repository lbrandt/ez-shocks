# -------------------------------------------------------------------------
# Calls the function svdraws to estimate a stochastic volatility model on
# prediction errors from the respective factor model estimated in MATLAB
#
# Lennart Brandt, Feb 2018
# -------------------------------------------------------------------------

# Set working directory to script file location

# Retrieve script location when code is run within RStudio
location.thisfile = dirname(rstudioapi::getActiveDocumentContext()$path)
# Retrieve script location if file is called via source()
location.thisfile = dirname(sys.frame(1)$ofile)
# Set working directory to file location
setwd(location.thisfile)

# Set location of MATLAB files relative to working directory which contains this script
location.matlab = normalizePath(file.path("..", "MATLAB"), winslash = "/")


require(stochvol)
require(h5)
options(digits = 17)
set.seed(840) # for replication

source("svdraws.R")
mcdraw = 50000
mcburn = 50000
mcthin = 10


# Eurozone stochastic volatility model ------------------------------------
in.file = h5file(file.path(location.matlab, "ez_factors_forc.mat"))
vyt = t(in.file["vyt"][]) # LASSO residuals
htvyt = t(in.file["htvyt"][]) # Hard thresholding residuals
vft = t(in.file["vft"][])
h5close(in.file)

# Stochvol draws
ez_vydraws = svdraws(vyt, mcdraw, mcburn, mcthin, outprog = TRUE)
ez_htvydraws = svdraws(htvyt, mcdraw, mcburn, mcthin, outprog = TRUE)
ez_vfdraws = svdraws(vft, mcdraw, mcburn, mcthin, outprog = TRUE)

# Save results to HDF5
out.file = h5file("ez_svresults.h5", mode = "w")
out.file["/sy"] = ez_vydraws$h
out.file["/ty"] = ez_vydraws$t
out.file["/syht"] = ez_htvydraws$h
out.file["/tyht"] = ez_htvydraws$t
out.file["/sf"] = ez_vfdraws$h
out.file["/tf"] = ez_vfdraws$t
h5close(out.file)




# Germany stochastic volatility model -------------------------------------
in.file = h5file(file.path(location.matlab, "de_factors_forc.mat"))
vyt = t(in.file["vyt"][])
vft = t(in.file["vyt"][])
h5close(in.file)

# Stochvol draws
de_vydraws = svdraws(vyt, mcdraw, mcburn, mcthin, outprog = TRUE)
de_vfdraws = svdraws(vft, mcdraw, mcburn, mcthin, outprog = TRUE)

# Save results to HDF5
out.file = h5file("de_vyresults.h5", mode = "w")
out.file["/sy"] = de_vydraws$h
out.file["/ty"] = de_vydraws$t
out.file["/sf"] = de_vfdraws$h
out.file["/tf"] = de_vfdraws$t
h5close(out.file)


