% -------------------------------------------------------------------------
% Romer & Romer 2004
% Monetary Policy Shocks
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB')

% Load data
filepath = 'C:\Dateien\Dropbox\Lennart\Thesis\Data\RomerandRomerDataAppendix.xls';

[data, ~, ~] = xlsread(filepath, 'DATA BY MEETING');

[~, names, ~] = xlsread(filepath, 'DATA BY MEETING', 'B1:Z1');
[~, dates, ~] = xlsread(filepath, 'DATA BY MEETING', 'A:A');

dates = datetime(dates(2:end));


% Set R&R sample from 1969-01-14 to 1996-12-17 (FULL)
ta = '1969-01-14';
te = '1996-12-17';
taindex = find(dates == ta);
teindex = find(dates == te);

sample = taindex:teindex;

% Plot change of intended target DTARG
plot(dates(sample), data(sample, strcmp(names, {'DTARG'})))



% Compute shock series as residual from equation (1)

% Select variables for estimation
depvar = data(sample, strcmp(names, {'DTARG'}));

% Matrix of independent variables
X = data(sample, 5:22);

% Delete NaN rows
X(any(isnan(X), 2), :) = [];

% OLS
beta = ols(depvar, X, 1);









% Historically implemented FFR
ffrpath = 'C:\Dateien\Dropbox\Lennart\Thesis\Data\ffr.xls';

[ffr, ~, ~] = xlsread(ffrpath);
[~, ffrdates, ~] = xlsread(ffrpath, 'A:A');

ffrdates = datetime(ffrdates(12:end));

ffrtaindex = find(ffrdates == ta);
ffrteindex = find(ffrdates == te);

ffrsample = ffrtaindex:ffrteindex;





