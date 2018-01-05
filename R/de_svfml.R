# -------------------------------------------------------------------------
# Estimate a first-order autoregressive stochastic volatility model on the
# forecast errors from an AR(4) model on each predictor (series by series)
# -------------------------------------------------------------------------

# ML Test

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



# Errors are product of first moment shock N(0,1) and volatility process logN(sigma(t), tau)
plot(1:obs.T, vt[[1]], type = "l")


# From https://support.sas.com/rnd/app/ets/examples/emmweb/index.htm

# Generate data
sim.M = 100 # Iterations
sim.T = 100 # Observations

# Params
a = -0.736
b = 0.9
s = 0.363

#
ll = sqrt( exp(a/(1-b)) )

for (i in 1:sim.T){
  
  u = rnorm(1, 0, 1)
  z = rnorm(1, 0, 1)
  
  lnssq = a + b*log(ll^2) + s*u
  
  st = sqrt(exp(lnssq))
  ll = st
  
  y[i] = st* z
}




