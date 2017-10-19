function [Fhat,Lhat,ehat,eigval] = factors(x, k, demean)
% -------------------------------------------------------------------------
% Estimates static factor model via Principal Component Analysis of the
% Variance-Covariance-Matrix of some data matrix.
%
% Transforms data before Eigendecomposition:
% If demean == 0, do not transform, take raw data.
% If demean == 1, only demean data.
% If demean == 2, standardise data to (0, 1).
%
%
%   Input
%       x           Matrix of observed variables [T x N]
%       nfac        Number of extracted factors [1 x 1]
%
%   Output
%       Fhat        Estimates of latent factors [T x k]
%       Lhat        Corresponding factor loadings [N x k]
%       ehat        Idiosyncratic component [T x N]
%       eigval      Unscaled Eigenvalues of data covariance matrix [N x 1]
% 
%   Dependencies {source}
%       svd {vanilla}
% -------------------------------------------------------------------------

[T, N] = size(x);

%%%%
% Transform data according to demean
if demean == 2
    xtr = (x - repmat(mean(x), T, 1))./repmat(std(x), T, 1);
    
elseif demean == 1
    xtr = x - repmat(mean(x), T, 1);
    
elseif demean == 0
    xtr = x;
end


%%%%
% Eigenvalue Decomposition of Covariance Matrix to estimate unrestricted factor space
xx = xtr'*xtr;

[EV, S, ~] = svd(xx);
eigval = diag(S);

% Extract nfac factors
Lhat = sqrt(N)*EV(:, 1:k);
Fhat = xtr*Lhat /N;
ehat = xtr - Fhat*Lhat';

end
