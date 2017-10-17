function [ehat,fhat,lambda,ev]=jln_pc(y,nfac)
% -------------------------------------------------------------------------
% wyd 
%
%   Input
%       y       data [T x N]
%       nfac    number of extracted factors [1 x 1]
%
%   Output
%       ehat    idiosyncratic errors
%       fhat    estimator of latent factors
%       lambda  factor loadings
%       ev      eigenvalues of data matrix 
% 
%   Dependencies {source}
%       svd {vanilla}
% -------------------------------------------------------------------------

[T,N] = size(y);
yy = y'*y;
[Fhat0,eigval,Fhat1] = svd(yy);
lambda = Fhat0(:,1:nfac)*sqrt(N);
fhat = y*lambda/N;
ehat = y - fhat*lambda';

ve2 = sum(ehat' .* ehat')'/N;
ev = diag(eigval);

end