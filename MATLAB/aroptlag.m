function [pminAIC, pminBIC, pminHQC, resultstable] = aroptlag(y, pmax, const, trend, out)
% -------------------------------------------------------------------------
% Determines the optimal lag order of an autoregressive model via the AIC,
% BIC and Hannan-Quinn criterion. Outputs a structure which contains values
% of the ICs at the respective lag lengths.
%
%   Input:  y       Series of interest
%           pmax    Maximum number of lags to consider
%           const   Indicates intercept in regression {0, 1}
%           trend   Indicates trend in regression {0, 1}
%
%   Output: pminAIC, pminBIC, pminHQC
% -------------------------------------------------------------------------

T = length(y);

if const == 1
    intercept = ones(T-pmax, 1);
else
    intercept = [];
end

if trend == 1
    timeindex = (1:T-pmax)';
else
    timeindex = [];
end

% Set up RHS variables
ylags = mlag(y, pmax);
x = [intercept, timeindex, ylags(pmax+1:end, :)];

llv = zeros(pmax, 1);
aicv = zeros(pmax, 1);
bicv = zeros(pmax, 1);
hqcv = zeros(pmax, 1);
for p = 1:pmax
    % Perform univariate AR regression with respective lag order
    regvars = x(:, 1:(const + trend + p));
    arparam = (regvars'*regvars)\regvars'*y(pmax+1:end);
    arresid = y(pmax+1:end) - regvars*arparam;
    
    % Compute values of information criteria and log-likelihood
    [aic, bic, hqc, ll] = aric(arresid, p, const, trend);
    llv(p) = ll;
    aicv(p) = aic;
    bicv(p) = bic;
    hqcv(p) = hqc;
end

[minAIC, pminAIC] = min(aicv);
[minBIC, pminBIC] = min(bicv);
[minHQC, pminHQC] = min(hqcv);

% Output results if indicated
if out == 1
    fprintf('\n Optimal lag order, AIC: p = %d at AIC = %f \n', pminAIC, minAIC);
    fprintf('\n Optimal lag order, BIC: p = %d at BIC = %f \n', pminBIC, minBIC);
    fprintf('\n Optimal lag order, HQC: p = %d at HQC = %f \n', pminHQC, minHQC);
    resultstable = table(llv, aicv, bicv, hqcv);
    display(resultstable);
end

end
