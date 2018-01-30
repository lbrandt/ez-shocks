function results = irf_jorda(depvarLeads, depvarLags, x, controls, l)
% -------------------------------------------------------------------------
% Estimates an impulse response function following Jorda's local projection
% approach. (AER, 2005, Vol. 95, No. 1, pp. 161-182)
%
%   Input
%       yLeads      Matrix of leads of variable of interest [T x H]
%       yLags       Matrix of lags of variable of interest [T x P]
%       x           Series of exogenous shocks [T x 1]
%       controls    Matrix of other controls, e.g. intercept, trends, dummies
%       l           Number of lags for Newey-West correction
%       
%
%   Output
%       results     Structure containing:
%         results.T         Number of time series observations
%         results.H         Maximum horizon of impulse response
%         results.P         Number of lags of variable of interest
%         results.theta     Impact multipliers
%         results.sigma     HAC standard errors of Theta estimators
%         results.irf       Impulse response function
%
%
%   Dependencies {source}
%       nwest {ls}
%
% -------------------------------------------------------------------------

% Setup
[T, H] = size(depvarLeads);
[T1, P] = size(depvarLags);
[T2, N] = size(x);
[T3, ~] = size(controls);

if isequal(T, T1, T2, T3) == 0
    error('Dimensions must align.');
end

if N ~= 1
    error('Must supply only one shock series.');
end

% Default Newey-West lag order
if nargin < 5 || isempty(l)
    warning('Defaults to N&W (1994) rule-of-thumb if lag order not specified.')
    l = fix(4*(T/100)^(2/9));
end


% ----------------
% IRF via Jorda's local projection method
theta = zeros(H, 1);
sigma = zeros(H, 1);

for i = 1:H
    depvar = depvarLeads(:, i);
    regression = nwest(depvar, [x, depvarLags, controls], l);
    % Select impact multiplier from parameter vector
    theta(i) = regression.beta(1);
    % Recover HAC standard errors used in function via $results.tstat = results.beta./nwerr;
    sigma(i) = regression.beta(1)/regression.tstat(1);
end

results.T = T;
results.H = H;
results.P = P;
results.theta = theta;
results.sigma = sigma;
results.irf = cumsum(theta);
end
