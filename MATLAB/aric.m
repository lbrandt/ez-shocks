function [aic, bic, hqc, ll] = aric(e, p, const, trend)
% -------------------------------------------------------------------------
% Computes information criteria of an autoregressive model assuming
% normally distributed residuals. Computes the log-likelihood based on the 
% biased, MSE-efficient variance estimator, i.e. S2 = mean((y-yhat)^2).
%
%   Input:  x       Vector of model residuals
%           p       Number of lags
%           const   Indicates intercept in regression {0, 1}
%           trend   Indicates trend in regression {0, 1}
%
%   Output:
% -------------------------------------------------------------------------
T = length(e);
k = const + trend + p;

ll = (-T/2)*log(2*pi*det((e'*e)/T)) - T/2;

aic = -2/T*ll + 2*k/T;
bic = -2/T*ll + k*log(T)/T;
hqc = -2/T*ll + 2*k*log(log(T))/T;

end