function [ehat,Fhat,lamhat,ve2] = jln_factors(X,kmax,jj,DEMEAN)
% -------------------------------------------------------------------------
% Estimate latent factors from observable data X using PCA
%   Input:  X      = observable data
%           kmax   = maximum number of factors to consider
%           jj     = information criterion from Bai and Ng (2002)
%           DEMEAN = 0 - no, 1 - yes, 2 - standardize
%   Output: ehat   = idiosyncratic errors
%           Fhat   = latent factor estimates
%           lamhat = factor loadings
%           ve2    = eigenvalues of data covariance matrix
% -------------------------------------------------------------------------

[ic1,chat,fhat,eigval]  = jln_nbplog(X,kmax,jj,DEMEAN);

icstar    = ic1;
R2_static = sum(eigval(1:icstar))/sum(eigval);

if DEMEAN == 2
    [ehat,Fhat,lamhat,ve2]  = jln_pc(jln_standard(X),icstar);
end
if DEMEAN == 1
    [ehat,Fhat,lamhat,ve2]  = jln_pc(X-repmat(mean(X),size(X,1),1),icstar); 
end
if DEMEAN == 0
    [ehat,Fhat,lamhat,ve2]  = jln_pc(X,icstar); 
end

end
