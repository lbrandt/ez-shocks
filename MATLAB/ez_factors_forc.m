% -------------------------------------------------------------------------
% Generate forecast via factor model
%
%
%
%
% -------------------------------------------------------------------------

%clear; clc;

% Load and manipulate data in import_data. Call results here:
load ez_data


%%%%
% Find optimal number of factors according to Bai & Ng (2002)
kmax   = 30; % Max number of factors to be extracted
gnum   = 2; % ICp2 chosen in JLN2015
demean = 2; % Standardise data

bnicv = zeros(kmax,1);
for k = 1:kmax
    bnicv(k) = bnic(x, k, gnum, demean); % Compute BNIC for each k
end

bnicmin = min(bnicv);
rhat = minind(bnicv); % Optimal number of factors according to lowest IC
fprintf('\nFactors via IC(%d): rhat = %d \n', gnum, rhat);


%%%%
% Extract factors via PCA
[Fhat, LF, ef, evf] = factors(x, rhat, demean);
[Ghat, LG, eg, evg] = factors(x.^2, rhat, demean);

sumeigval = cumsum(evf)/sum(evf);
R2_static = sum(evf(1:rhat))/sum(evf);


% AR order of Fhat
pmax = 24;

const = 1;
trend = 0;

picf = zeros(rhat, 3);
icf = zeros(rhat, 3);
for i = 1:rhat
    [picf(i, 1), icf(i, 1)] = aroptlag(Fhat(:, i), pmax, 'aic', const, trend, 0);
    [picf(i, 2), icf(i, 2)] = aroptlag(Fhat(:, i), pmax, 'bic', const, trend, 0);
    [picf(i, 3), icf(i, 3)] = aroptlag(Fhat(:, i), pmax, 'hqc', const, trend, 0);
end
% Maximum lag length suggested by ICs is p = 4. Err on the side of caution?
%aroptlag(Fhat(:, 1), pmax, [], 1, 0, 1);



%%%%
% Forecast

% Build predictor set
zt       = [Fhat, Fhat(:, 1).^2, Ghat(:, 1)];
[~, M]   = size(zt);

% Set dependent variables
yt       = standardise(x);
[T, N]   = size(yt);

py       = 4; % number of depvar lags
pz       = 4; % number of predictor lags
maxlag   = max(py, pz);

L        = fix(4*(T/100)^(2/9)); % Newey-West lag length rule-of-thumb (N&W 1994)


% LASSO model selection
lambdamin = 0;
lambdamax = 3;
tmin = 100;
const = 0;
roll = 0;

% Set up optimisation problem
obfunopt = optimset('TolFun', 0.01, 'TolX', 0.01); % Increase tolerance

ylambda = zeros(1, N);
ymsemin = zeros(1, N);
ymodels = zeros(1 + py + pz*M, N); % Indicator matrix of nonzero predictors after penalisation
ybetas  = zeros(1 + py + pz*M, N); % LASSO parameter vectors of single equations in columns
yfit    = zeros(T - maxlag, N); % Fitted values 
vyt     = zeros(T - maxlag, N); % Forecast errors
for j = 1:N % Estimate system equation-by-equation
    tic;
    
    X    = [ones(T, 1), mlag(yt(:, j), py), mlag(zt, pz)]; % const + py lags of depvar + pz lags of predictors
    % Write static parameters to workspace
    yopt = yt(maxlag+1:end, j);
    xopt = X(maxlag+1:end, :);
    % Find optimal lambda
    obfunpar = @(lopt, yopt, xopt, tmin, const, roll)forcmseLasso(lopt, yopt, xopt, tmin, const, roll); % Anonymous function with all params
    obfunlam = @(lopt)obfunpar(lopt, yopt, xopt, tmin, const, roll); % Anonymous function with only one parameter, taking other params from workspace
    [lambda, msemin, ~, ~] = fminbnd(obfunlam, lambdamin, lambdamax, obfunopt);
    % Estimate full model given optimal lambda
    yLasso = solveLasso(yt(maxlag+1:end, j), X(maxlag+1:end, :), lambda);
    
    ybetas(:, j) = yLasso.beta;
    yfit(:, j) = yLasso.X*yLasso.beta;
    vyt(:, j) = yLasso.y - yLasso.X*yLasso.beta;
    
    ylambda(1, j) = lambda;
    ymsemin(1, j) = msemin;
    ymodels(:, j) = yLasso.beta ~= 0;
    
    fprintf('Series %d, Elapsed Time = %0.4f \n', j, toc);
end

% Compute MSE of LASSO generated prediction errors
vymse = mean(vyt.^2);



% Hard thresholding model selection
htmodels = zeros(1 + py + pz*M, N); % Indicator matrix of included predictors
htbetas  = zeros(1 + py + pz*M, N); % OLS parameter vectors of single equations in columns
htfit    = zeros(T - maxlag, N); % Fitted values 
htvyt    = zeros(T - maxlag, N); % Forecast errors
for j = 1:N % Estimate system equation-by-equation
    
    X    = [ones(T, 1), mlag(yt(:, j), py), mlag(zt, pz)]; % const + py lags of depvar + pz lags of predictors
    reg  = nwest(yt(maxlag+1:end, j), X(maxlag+1:end, :), L);
    pass = abs(reg.tstat(py+2:end)) > 2.575; % hard threshold 
    keep = [ones(1, py+1) == 1, pass']; % always keep const, depvar lags, F1t
    Xnew = X(:, keep);
    reg  = nwest(yt(maxlag+1:end, j), Xnew(maxlag+1:end, :), L);
    
    htbetas(keep, j) = reg.beta;
    htfit(:, j)      = reg.yhat;
    htvyt(:, j)      = reg.resid;
    
    htmodels(:, j)   = keep;
end

% Compute MSE of LASSO generated prediction errors
htvymse = mean(htvyt.^2);

% Relative MSE efficiency of LASSO models
remse = htvymse./vymse;

% Number of nonzero parameters retained under both schemes
yparnum = sum(ymodels);
hparnum = sum(htmodels);
rparnum = yparnum./hparnum;

yparbar = mean(yparnum);
hparbar = mean(hparnum);

yparbar/hparbar

summarize(vyt);


% Generate AR(4) errors for Predictor set zt


[T, R]   = size(zt);
pf       = 4;
L        = fix(4*(T/100)^(2/9));

fbetas   = zeros(1 + pf, R); % Parameter vectors of single equations in columns
ffit     = zeros(T - pf, R); % Fitted values
vft      = zeros(T - pf, R); % Prediction errors

for j = 1:R % Equation-by-equation
    X    = [ones(T, 1), mlag(zt(:, j), pf)];
    reg  = nwest(zt(pf+1:end, j), X(pf+1:end, :), L);
    
    fbetas(:, j) = reg.beta;
    ffit(:, j)   = reg.yhat;
    vft(:, j)    = reg.resid;
end


%%%%
% Save data
maxlag = max([py, pz, pf]); % Maximum lag length out of all regressions run in file
dates = dates(1+maxlag:end);

save ez_factors_forc -v7.3 dates yfit ffit ybetas fbetas vyt vft names py pz pf zt x ymodels htmodels htbetas htvyt
