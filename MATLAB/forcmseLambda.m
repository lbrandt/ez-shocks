function [lambdaopt, msemin, msevec] = forcmseLambda(y, x, lambdavec, tmin, const, roll)
% -------------------------------------------------------------------------
% Finds optimal value for LASSO regularisation parameter lambda in a given
% data set via a quasi out-of-sample forecast experiment. Picks that lambda
% which leads to minimal MSE in one-step-ahead forecast errors.
%
%   Input
%       y           Vector of dependent variable [T x 1]
%       x           Matrix of predictors [T x N]
%       lambdavec   Vector of possible parameter values [1 x M]
%       tmin        Size of minimum holdout sample
%       const       Add intercept to regressor matrix {1, 0}
%       roll        Choose rolling or expanding window {1, 0}
%
%   Output
%       lambdaopt   Optimal lambda
%       msemin      Minimal MSE at lambdaopt
%       msevec      Vector of MSE values corresponding to lambdavec
%
%   Dependencies {source}
%       forcerrorsLasso {lb}
%       solveLasso {Pendse}
%
% -------------------------------------------------------------------------
M = length(lambdavec);
lambdamin = min(lambdavec);
lambdamax = max(lambdavec);

msevec = zeros(1, M);
for i = 1:M
    lambdaforc = forcerrorsLasso(y, x, lambdavec(i), tmin, const, roll);
    msevec(i) = mean(lambdaforc.fe.^2);
end

[msemin, ind] = min(msevec);
lambdaopt = lambdavec(ind);
end
