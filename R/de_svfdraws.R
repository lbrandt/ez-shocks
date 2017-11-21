# -------------------------------------------------------------------------
# Estimate a first-order autoregressive stochastic volatility model on the
# forecast errors from an AR(4) model on each predictor (series by series)
# -------------------------------------------------------------------------

#rm(list=ls())


# Initialization
require(stochvol)
options(digits=17)
set.seed(1000) # for replication

# KE file
vt = read.csv(header = TRUE, sep = ",", "C:/Dateien/My Dropbox/Lennart/Thesis/Code/ez-shocks/MATLAB/de_factors_vft.csv")

# Remove dates vector
vt = vt[, -1]

obs.T    = dim(vt)[1]
obs.N    = dim(vt)[2]


for (i in 1:obs.N){
  
  if(min(log(vt[, i]^2)) == -Inf){
    
    vt[, i] = vt[, i] + 0.00001 #offset to avoid taking log of zero
  }
}


# Run MCMC algorithm and store draws
S    = 50000
burn = 50000 # MCMC with [50000, 50000] takes roughly 50 seconds per iteration.

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
write.csv(h, file = 'de_svflatent.csv', row.names = FALSE)
write.csv(t, file = 'de_svfparams.csv', row.names = FALSE)
write.csv(g, file = 'de_svfgeweke.csv', row.names = FALSE)
