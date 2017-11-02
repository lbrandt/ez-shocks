function [evar] = expectvar(st, alpha, beta, tau, h)
% -------------------------------------------------------------------------
% Computes the h-step-ahead expected variance E[exp{s(t+h)}] of the 
% observed process with respect to information set I(t) using the
% AR(1) law of motion for the unobserved (latent) process s(t). Assume:
%
%       s(t) = log{ sigma^2(t) } = alpha + beta* s(t-1) + tau* u(t)
%
%       with |beta| < 1 and u ~ N(0, 1).
%
%   Input
%       st          Vector of estimators of the latent process s(t) [T x 1]
%       alpha       Conditional mean of the logvariance process
%       beta        AR(1) parameter of the logvariance process
%       tau         Volatility of the logvariance process
%       h           Forecast horizon
%
%   Output
%       evar        Expected h-step-ahead forecast variance [T x 1]
% 
%   Dependencies {source}
% -------------------------------------------------------------------------

evar = exp( alpha* (1-beta^h)/(1-beta) + beta^h * st + 0.5*tau^2 * (1-beta^(2*h))/(1-beta^2) );

end
