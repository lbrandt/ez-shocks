% -------------------------------------------------------------------------
% Generate forecast via factor model
%
%
%
%
% -------------------------------------------------------------------------

%clear; clc;

% Load and manipulate data in import_data. Call script here:
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

ybetas   = zeros(1 + py + pz*M, N); % Parameter vectors of single equations in columns
yfit     = zeros(T - maxlag, N); % Fitted values 
vyt      = zeros(T - maxlag, N); % Forecast errors

% LASSO model selection
lambdavec = linspace(0, 4, 100);
tmin = 100;
const = 0;
roll = 0;

ylambda  = zeros(1, N);
ymodels  = zeros(1 + py + pz*M, N); % Indicator matrix of nonzero predictors after penalisation
for j = 1:N % Estimate system equation-by-equation
    tic;
    
    X    = [ones(T, 1), mlag(yt(:, j), py), mlag(zt, pz)]; % const + py lags of depvar + pz lags of predictors
    lambda = forcmseLambda(yt(maxlag+1:end, j), X(maxlag+1:end, :), lambdavec, tmin, const, roll); toc
    yLasso = solveLasso(yt(maxlag+1:end, j), X(maxlag+1:end, :), lambda);
    
    ybetas(:, j) = yLasso.beta;
    yfit(:, j) = yLasso.X*yLasso.beta;
    vyt(:, j) = yLasso.y - yLasso.X*yLasso.beta;
    
    ylambda(1, j) = lambda;
    ymodels(:, j) = yLasso.beta ~= 0;
    
    fprintf('Series %d, Elapsed Time = %0.4f \n', i, toc);
end

summarize(vyt);
summarize(vyols);





% Hard thresholding model selection
ymodels  = zeros(1 + py + pz*M, N); % Indicator matrix of included predictors
for j = 1:N % Estimate system equation-by-equation
    
    X    = [ones(T, 1), mlag(yt(:, j), py), mlag(zt, pz)]; % const + py lags of depvar + pz lags of predictors
    reg  = nwest(yt(maxlag+1:end, j), X(maxlag+1:end, :), L);
    pass = abs(reg.tstat(py+2:end)) > 2.575; % hard threshold 
    keep = [ones(1, py+1) == 1, pass']; % always keep const, depvar lags, F1t
    Xnew = X(:, keep);
    reg  = nwest(yt(maxlag+1:end, j), Xnew(maxlag+1:end, :), L);
    
    ybetas(keep, j) = reg.beta;
    yfit(:, j)      = reg.yhat;
    vyt(:, j)       = reg.resid;
    
    ymodels(:, j)   = keep;
end


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

save ez_factors_forc -v7.3 dates yfit ffit ybetas fbetas vyt vft names py pz pf zt x ymodels
