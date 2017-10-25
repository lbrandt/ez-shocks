% -------------------------------------------------------------------------
% Generate forecast via factor model
%
%
%
%
% -------------------------------------------------------------------------

clear; clc;

% Load and manipulate data in import_data. Call results here:
import_data


% Simulate data as in testnfac.m by Serena Ng
%clear; clc;
rng(999);

r = 4; 
N = 100;
T = 50;

e = randn(T,N);
f = randn(T,r);
lambda = randn(N,r);
x = f*lambda'+ e;

fprintf('Data generating process: r = %d \n', r);
fprintf('Sample: T = %d, N = %d \n', size(x));


%%%%
% Find optimal number of factors according to Bai & Ng (2002)
kmax   = 20; % Max number of factors to be extracted
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

%[ehat1, Fhat1, lamhat1, ev1]  = jln_factors(x, kmax, gnum, demean);
%[ehat1, Ghat1, lamhat1, ev1]  = jln_factors(x.^2, kmax, gnum, demean);

[Fhat, LF, ef, evf] = factors(x, rhat, demean);
[Ghat, LG, eg, evg] = factors(x.^2, rhat, demean);

sumeigval = cumsum(evf)/sum(evf);
R2_static = sum(evf(1:rhat))/sum(evf);


%%%%
% Forecast

% Build predictor set
zt       = [Fhat, Fhat(:, 1).^2, Ghat(:, 1)];
[~, M]   = size(zt);

% Set dependent variables
yt       = standardise(x(:, 1:132)); % only macro data
[T, N]   = size(yt);

py       = 4; % number of depvar lags
pz       = 2; % number of predictor lags
maxlag   = max(py, pz);

L        = fix(4*(T/100)^(2/9)); % Newey-West lag length rule-of-thumb (N&W 1994)

ybetas   = zeros(1 + py + pz*M, N); % Parameter vectors of single equations in columns
yfit     = zeros(T - maxlag, N); % Fitted values 
vyt      = zeros(T - maxlag, N); % Forecast errors

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

save factors_forc dates yfit ffit ybetas fbetas vyt vft names vartype py pz pf zt x ymodels

% Also write to .txt file for R code
%dlmwrite('factors_vyt.txt',vyt,'delimiter','\t','precision',17);
%dlmwrite('factors_vft.txt',vft,'delimiter','\t','precision',17);


% Write to .csv for further use
datetable = array2table(string(dates), 'VariableNames', {'Date'});


% vyt factor model prediction errors
ynames = varnames(1:132); % only macro variables

vytable = array2table(vyt, 'VariableNames', ynames);
vytable = [datetable, vytable];

writetable(vytable, 'factors_vyt.csv');


% vft factor AR(4) prediction errors
fstring = string(repmat('Factor', R, 1));
fnames = strcat(fstring, string((1:R)'));
fnames = cellstr(fnames);

vftable = array2table(vft, 'VariableNames', fnames);
vftable = [datetable, vftable];

writetable(vftable, 'factors_vft.csv');


