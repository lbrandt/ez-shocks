# -------------------------------------------------------------------------
# Estimate a first-order autoregressive stochastic volatility model on the
# forecast errors from an AR(4) model on each predictor (series by series)
# -------------------------------------------------------------------------

#rm(list=ls())

# Set working directory to script file location

# Retrieve script location when code is run within RStudio
location.thisfile = dirname(rstudioapi::getActiveDocumentContext()$path)
# Retrieve script location if file is called via source()
location.thisfile = dirname(sys.frame(1)$ofile)
# Set working directory to file location
setwd(location.thisfile)

# Set location of MATLAB files relative to working directory which contains this script
location.matlab = normalizePath(file.path("..", "MATLAB"), winslash = "/")




# Initialization
require(h5) # Output to HDF5
require(stochvol)
options(digits = 17)
set.seed(840) # for replication


# MATLAB file
in.file = h5file(file.path(location.matlab, "de_factors_forc.mat"))
vt = t(in.file["vft"][])
h5close(in.file)


obs.T    = dim(vt)[1]
obs.N    = dim(vt)[2]


for (i in 1:obs.N){
  
  if(min(log(vt[, i]^2)) == -Inf){
    
    vt[, i] = vt[, i] + 0.00001 #offset to avoid taking log of zero
  }
}


# Run MCMC algorithm and store draws
S    = 50000
burn = 50000 # MCMC with [50000, 50000] takes roughly 10 seconds per iteration.

h    = matrix(0, obs.T, obs.N) # Latent variable ht
t    = matrix(0, 3, obs.N) # Parameter vector theta
g    = matrix(0, 3, obs.N) # Geweke convergence statistic for each parameter

for (i in 1:obs.N){
  
  draws   = svsample(vt[, i], draws = S, burnin = burn, quiet = TRUE, thinpara = 10, thinlatent = 10)
  
  h[, i]  = colMeans(draws$latent)
  t[, i]  = colMeans(draws$para)
  g[, i]  = geweke.diag(draws$para)$z

  # Show progress in console
  cat('i =', i)
  print(Sys.time())
}


# Save results to .csv in format [estimators, variables]
# write.csv(h, file = 'gs_svflatent.csv', row.names = FALSE)
# write.csv(t, file = 'gs_svfparams.csv', row.names = FALSE)
# write.csv(g, file = 'gs_svfgeweke.csv', row.names = FALSE)


# Save results to HDF5
out.file = h5file("de_svfresults.h5", mode = "w")
out.file["/h"] = h
out.file["/t"] = t
out.file["/g"] = g
h5close(out.file)
