# -------------------------------------------------------------------------
# Estimates a first-order autoregressive stochastic volatility model on the
# forecast errors from a predictive time series model. Adapted from code
# supplied by Jurado, Ludvigson, Ng (2015). Uses the Kastner's (2016)
# implementation of the stochastic volatility model in the CRAN package
# "stochvol".
#
# Input   vt          Matrix of residuals
#         mcdraw      Number of used Monte Carlo draws
#         mcburn      Number of burn-in draws
#         mcthin      Thinning parameter (default = 1, no thinning)
#         outproc     Display progress in console (default = FALSE)
#         outfile     Name of HDF5 file to save results
#
# Output  results$h   Estimates of the latent volatility process h(t)
#         results$t   Parameter vector theta
#         results$g   Geweke convergence statistics
#         
# -------------------------------------------------------------------------
svdraws = function(vt, mcdraw, mcburn, mcthin = 1, outprog = FALSE, outfile = NULL){
  
  # Initialization
  checkInstall.stochvol = require(stochvol)
  if(checkInstall.stochvol == FALSE){
    warning("CRAN package 'stochvol' is not installed.")
  }
  
  obs.T    = dim(vt)[1]
  obs.N    = dim(vt)[2]
  
  for(i in 1:obs.N){
    if(min(log(vt[, i]^2)) == -Inf){
      vt[, i] = vt[, i] + 0.00001 #offset to avoid taking log of zero
    }
  }
  
  # Run MCMC algorithm and store draws
  draw = mcdraw
  burn = mcburn # MCMC with [50000, 50000] takes roughly 10 seconds per iteration.
  
  h    = matrix(0, obs.T, obs.N) # Latent variable ht
  t    = matrix(0, 3, obs.N) # Parameter vector theta
  g    = matrix(0, 3, obs.N) # Geweke convergence statistic for each parameter
  for(i in 1:obs.N){
    
    draws  = svsample(vt[, i], draws = mcdraw, burnin = mcburn, quiet = TRUE, thinpara = mcthin, thinlatent = mcthin)
    h[, i] = colMeans(draws$latent)
    t[, i] = colMeans(draws$para)
    g[, i] = geweke.diag(draws$para)$z
    
    if(outprog == TRUE){
      # Show progress in console
      cat('Series i =', i)
      print(Sys.time())
    }
  }
  
  # Return results
  svresults = list(h = h, t = t, g = g)
  if(is.null(outfile)){
    return(svresults)
  }else{
    return(svresults)
    
    # Save results to HDF5
    checkInstall.h5 = require(h5)
    if(checkInstall.h5 == FALSE){
      warning("CRAN package 'h5' is not installed. Needed to output results to HDF5.")
    }
    out.file = h5file(outfile, mode = "w")
    out.file["/h"] = svresults$h
    out.file["/t"] = svresults$t
    out.file["/g"] = svresults$g
    h5close(out.file)
  }
}
