# -------------------------------------------------------------------------
# Estimate a first-order autoregressive stochastic volatility model on the
# forecast errors from an AR(4) model on each predictor (series by series)
# -------------------------------------------------------------------------

# Initialization
rm(list=ls())
require(stochvol)
set.seed(1000) # for replication
options(digits=17)

vt = read.csv(header = TRUE, sep = ",", "C:/Dateien/My Dropbox/Lennart/Thesis/Code/ez-shocks/MATLAB/factors_vft.csv")

# Remove dates vector
vt = vt[, -1]

T    = dim(vt)[1]
N    = dim(vt)[2]

for (i in 1:N){
	if(min(log(vt[,i]^2))== -Inf){
		vt[,i] = vt[,i] + 0.00001 #offset to avoid taking log of zero
	}
}

# Run MCMC algorithm and store draws
S    = 500
burn = 500
m    = matrix(0,T+3,N)
g    = matrix(0,3,N)
for (i in 1:N){
	draws  = svsample(vt[,i],draws=S,burnin=burn,quiet=TRUE,thinpara=10,thinlatent=10)
	all    = cbind(draws$para,draws$latent)
	m[,i]  = colMeans(all)
	g[,i]  = geweke.diag(draws$para)$z
	name   = sprintf('svfdraws%d.txt',i)
#	write(t(all),file=name,ncolumn=dim(all)[2])
	
	# Show progress in console
	cat('i =', i)
	print(Sys.time())	
}

out = rbind(m,g) #include Geweke statistics
#write(t(out),file='svfmeans2014.txt',ncolumn=dim(out)[2])

write.csv(t(out), file = 'svfmeans2.csv', row.names = FALSE)
