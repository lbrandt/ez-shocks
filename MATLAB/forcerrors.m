% Forecast Experiment

% -------------------------------------------------------------------------
% Test h = 1 and no contemporaneous regressors
%clear; clc;
load ez_data
dip = x(:, findstrings(names, {'EKIPTOTG'}));

plot(dates, dip)
% Static regression
pdip = aroptlag(dip, 36, 'aic', 1, 0, 0);
diplags = mlag(dip, pdip);

T = length(dip);

dipmodel = ols(dip, [ones(T, 1), diplags]);
%prt(dipmodel);
plt(dipmodel);

% Inputs
y = dip(pdip+1:end);
x = diplags(pdip+1:end, :);
lambda = 3;
tmin = 80; % Holdout sample
const = 1;
roll = 1; % Rolling window, expanding == 0
%

forclasso = forcerrorsLasso(y, x, lambda, tmin, const, roll);


[nobs, ~] = size(y);
[~, nreg] = size(x);
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
        
        % Set sample for estimation
        ywindow = y(wstart:fstart-1, :);
        xwindow = x(wstart:fstart-1, :);
        forcmodel = ols(ywindow, xwindow);
        
        % Calculate one-step-ahead model forecast at fstart
        forcbetas = forcmodel.beta;
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
        
        % Set sample for estimation
        ywindow = y(wstart:fstart-1, :);
        xwindow = x(wstart:fstart-1, :);
        forcmodel = ols(ywindow, xwindow);
        
        % Calculate one-step-ahead model forecast at fstart
        forcbetas = forcmodel.beta;
        modelforc = x(fstart, :)*forcbetas;
        
        % Save results to structure
        results.yh(i) = y(fstart);
        results.ff(i) = modelforc;
        results.fe(i) = y(fstart) - modelforc;
    end
end






regstart = 2;

h = 1;
window = 100;

neq = 1;
festart = regstart;
festop = taindex + window;

femax = regend - festop + 1 - h;

ff = zeros(femax, 1);
fe = zeros(femax, 1);
fu = zeros(femax, 1);
fl = zeros(femax, 1);


for i = 1:femax
    
    sample = festart:festop;
    beta = ols(Y(sample), X(sample), 0);
    
    ff(i) = X(sample + h)* beta;
    
    fe(i) = Y(sample + h) - forc(i);
end

ssfe = sum(fe.^2);
mse = ssfe/femax;
rmse = mse^0.5;

% Constructing Point and Interval Forecasts
coverage = 0.6;


for i = 1:femax
    
    fu(i) = ff(i) + norm_inv(0.5 + 0.5*coverage)*rmse;
    fl(i) = ff(i) - norm_inv(0.5 + 0.5*coverage)*rmse;
    
end


%%%%%%

hdip = 5;
pmax = max(pfac, pdip);

diplags = mlag(dip, pdip);
diplead = lead(dip, hdip);

faclags = mlag(Fhat, pfac);
dipvars = [ones(T, 1), diplags, faclags];

holdout = 10;

% Full sample
fullsample = 1+pmax:T-hdip;

h1model = ols(dip(fullsample), dipvars(fullsample, :));
h2model = ols(diplead(fullsample, 1), dipvars(fullsample, :));
h3model = ols(diplead(fullsample, 2), dipvars(fullsample, :));


forcerrors = zeros(holdout, 1);
for i = 1:holdout
    
    i = 10
    
    fsample = pmax+i-1:(T-holdout-hdip+i); % Rolling window
    h1model = ols(dip(fsample, :), dipvars(fsample, :));
    h2model = ols(diplead(fsample, 1), dipvars(fsample, :));
    
    forcmodel.yhat(end)
    
    matchsample(dates, dates(T-holdout))
    
    ff1 = dipvars(T-holdout+i, :)*forcmodel.beta;
    ff2 = mlag(dipvars(fsample+1, :), 1)

    forecast
end
dipmodel = ols(dip, dipvars);
prt(dipmodel);

plot(dip)
hold on
plot(dipmodel.yhat)
