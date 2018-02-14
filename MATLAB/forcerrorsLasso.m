function results = forcerrorsLasso(y, x, lambda, tmin, const, roll)
% -------------------------------------------------------------------------
% Performs a quasi out-of-sample forecast experiment with parameter
% estimators penalised via LASSO. Utilises LASSO estimation function
% "solveLasso" by Gautam V. Pendse (McLean, Harvard Med School).
%
%   Input
%       y           Vector of dependent variable [T x 1]
%       x           Matrix of predictors [T x N]
%       lambda      LASSO regularisation parameter
%       tmin        Size of minimum holdout sample
%       const       Add intercept to regressor matrix {1, 0}
%       roll        Choose rolling or expanding window {1, 0}
%
%   Output
%       results     Structure containing:
%         results.method    Type of estimation = 'LASSO'
%         results.yh        Actual values of dependent variable
%         results.ff        Forecasted values
%         results.fe        Forecast errors
%
%   Dependencies {source}
%       solveLasso {Pendse}
%
% -------------------------------------------------------------------------
results.meth = 'LASSO';

[nobs, ~] = size(y);
nfor = nobs - tmin - 1;
% Add constant to regressor matrix if required
if const == 1
    x = [ones(nobs, 1), x];
end

results.yh = zeros(nfor, 1);
results.ff = zeros(nfor, 1);
results.fe = zeros(nfor, 1);
if roll == 1 % Rolling window
    for i = 1:nfor    
        wstart = i; % Window moves each iteration and has fixed length tmin
        fstart = tmin+i;
        
        % Set sample and estimate
        ywindow = y(wstart:fstart-1, :);
        xwindow = x(wstart:fstart-1, :);
        forclasso = solveLasso(ywindow, xwindow, lambda);
        
        % Calculate one-step-ahead model forecast at fstart
        forcbetas = forclasso.beta;
        modelforc = x(fstart, :)*forcbetas;
        
        % Save results to structure
        results.yh(i) = y(fstart);
        results.ff(i) = modelforc;
        results.fe(i) = y(fstart) - modelforc;
    end
elseif roll == 0 % Expanding window
    for i = 1:nfor    
        wstart = 1; % Always start window at 1
        fstart = tmin+i;
        
        % Set sample and estimate
        ywindow = y(wstart:fstart-1, :);
        xwindow = x(wstart:fstart-1, :);
        forclasso = solveLasso(ywindow, xwindow, lambda);
        
        % Calculate one-step-ahead model forecast at fstart
        forcbetas = forclasso.beta;
        modelforc = x(fstart, :)*forcbetas;
        
        % Save results to structure
        results.yh(i) = y(fstart);
        results.ff(i) = modelforc;
        results.fe(i) = y(fstart) - modelforc;
    end
end

end
