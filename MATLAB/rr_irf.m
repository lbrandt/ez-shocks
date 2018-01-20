% -------------------------------------------------------------------------
% Romer & Romer 2004
% Impulse Response
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB')

% Load data
filepath = 'C:\Dateien\Dropbox\Lennart\Thesis\Data\RomerandRomerDataAppendix.xls';

[data, ~, ~] = xlsread(filepath, 'DATA BY MONTH');

[~, names, ~] = xlsread(filepath, 'DATA BY MONTH', 'B1:Q1');
[~, dates, ~] = xlsread(filepath, 'DATA BY MONTH', 'A:A');

dates = datetime(dates(2:end));

dum12(month(dates) == 12) = 1; % Monthly Dummies

% Set R&R sample from 1970:1 to 1996:12
ta = datetime('1970-01-01');
te = datetime('1996-12-01');
taindex = find(dates == ta);
teindex = find(dates == te);

sample = taindex:teindex;

% Plot shock series RESID
plot(dates(sample), data(sample, 1))

% % Load constructed series
% load rr_shocks
% msample = (find(dates == mdates(1))):(find(dates == mdates(end)));
% 
% % Compare to RESID
% figure
% plot(mdates, mshocks)
% hold on
% plot(mdates, data(msample, strcmp(names, {'RESID'})))
% % Is same!


% Select variables for estimation
depvar = data(sample, strcmp(names, {'PCIPNSA'}));

% Lag order as in R&R 2004
const = 1;
dlag  = 11;
ylag  = 24;
slag  = 36;
maxlag = max([dlag, ylag, slag]);

% Not like this! Maybe this is why RATS gives different results?
% % Sample select before lag, i.e. obs < 324
% dum12lags  = mlag(dum12(sample)', dlag);
% depvarlags = mlag(data(sample, strcmp(names, {'PCIPNSA'})), slag);
% shocklags  = mlag(data(sample, strcmp(names, {'RESID'})), slag);
% 
% % Matrix of independent variables
% X = [dum12lags, depvarlags, shocklags];

% Lag before sample select, i.e. obs = 324 because presample counts
dum12lags  = mlag(dum12', dlag);
depvarlags = mlag(data(:, strcmp(names, {'PCIPNSA'})), ylag);
shocklags  = mlag(data(:, strcmp(names, {'RESID'})), slag);

% Matrix of independent variables
X = [dum12lags(sample, :), depvarlags(sample, :), shocklags(sample, :)];

% OLS
coeffs = ols(depvar, X, const);

% The estimated betas of the shock are the same as in paper. Yay!
c = coeffs((const + dlag + ylag + 1):end);

% Append zeros to autoregressive coeffs
b = [coeffs((const + dlag + 1):((const + dlag + ylag))); zeros(maxlag - ylag, 1)];

% Dimension companion matrix (incomplete VAR) for IRF
varnum = 2;

AT = zeros(varnum, varnum*maxlag);

% Reorder coefficients by lag order
for i = 1:maxlag
    AT(1, 2*i - 1) = b(i); % Only first row because we only first variable is endogenous
    AT(1, 2*i)     = c(i);
end

% Build companion matrix
A = companion(varnum, maxlag, AT);

% Compute IRF of variable [varselect] w.r.t. shock [shockselect] up until horizon [hmax]
varselect = 1;
shockselect = 2;

hmax = 48;

impact = zeros(hmax, 1);
for h = 1:hmax
    AH = A^h;
    impact(h) = AH(varselect, shockselect); % Impact multipliers
end

irf = cumsum(impact);

% Plot on page 1071 of R&R 2004
plot(irf)
