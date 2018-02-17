function results = irf_companion(y, x, P, Q, H, controls)
% -------------------------------------------------------------------------
% Estimates an impulse response function following the approach outlined in
% Sims' (2012) "Note On Computing Impulse Responses for an AR(p)".
%
%   Input
%       y           Series of variable of interest [T x 1]
%       x           Series of exogenous shocks [T x 1]
%       P           Number of lags in y
%       Q           Number of lags in x
%       H           Maximum impulse response horizon
%       controls    Matrix of other controls, e.g. intercept, trends, dummies
%
%   Output
%       results     Structure containing:
%         results.method    'VAR'
%         results.T         Maximum horizon of impulse response
%         results.theta     Impact parameters
%         results.irf       Cumulative impulse response function
%         results.sigma     %%NOT FUNCTIONAL%%
%
%   Dependencies {source}
%       mlag {ls}
%       ols {ls}
%       companion {lb}
%
% -------------------------------------------------------------------------

T = length(y);
V = length(x);

if isequal(T, V) == 0
    error('Dimensions must align.');
end

maxlag = max([P, Q]);

% Build matrix of RHS variables
yLags = mlag(y, P);
xLags = mlag(x, Q);
X = [yLags, xLags, controls];

% OLS
linreg = ols(y((maxlag + 1):end, :), X((maxlag + 1):end, :));
% Collect beta parameters and pad with zeros if lag lengths do not align
yBetas = [linreg.beta(1:P); zeros(maxlag - P, 1)];
xBetas = [linreg.beta((P + 1):(P + Q)); zeros(maxlag - Q, 1)];

% Reorder coefficients by lag order for VAR representation
varnum = 2;
AT = zeros(varnum, varnum*maxlag);
for i = 1:maxlag
    AT(1, 2*i - 1) = yBetas(i); % Only first row because only first variable is actually endogenous
    AT(1, 2*i)     = xBetas(i);
end
% Build companion matrix
A = companion(varnum, maxlag, AT);

% Compute IRF of variable [varselect] w.r.t. shock [shockselect] up until horizon [hmax]
varselect = 1;
shockselect = 2;

theta = zeros(H, 1);
for h = 1:H
    AH = A^h;
    theta(h) = AH(varselect, shockselect); % Impact multipliers
end

results.method = 'VAR';
results.T = T;
results.theta = theta;
results.irf = cumsum(theta);
%results.sigma = []; % Compute bootstrap errors and variance
end
