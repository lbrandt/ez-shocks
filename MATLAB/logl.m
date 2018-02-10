function ll = logl(e)
% -------------------------------------------------------------------------
% Computes the log-likelihood of a model assuming iid normal disturbances
% as inputs to the function. Computes based on biased, MSE-efficient
% variance estimator, i.e. S2 = mean((y-yhat)^2).
%
%   Input:  e       Vector of iid normal errors
%
%   Output:
% -------------------------------------------------------------------------
nobs = length(e);
ll = (-nobs/2)* log(2*pi*(e'*e)/nobs) - nobs/2;

end
