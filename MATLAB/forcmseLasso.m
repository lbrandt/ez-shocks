function mse = forcmseLasso(lambda, y, x, tmin, const, roll)
% -------------------------------------------------------------------------
% Computes MSE on set of LASSO generated forecast errors
% -------------------------------------------------------------------------
lambdaforc = forcerrorsLasso(y, x, lambda, tmin, const, roll);
mse = mean(lambdaforc.fe.^2);

end
