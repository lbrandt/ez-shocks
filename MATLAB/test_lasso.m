% Variables
%clear; clc;
load ez_data
dip = x(:, findstrings(names, {'EKIPTOTG'}));
%pdip = aroptlag(dip, 36, 'aic', 1, 0, 0);
pdip = 12;
diplags = mlag(dip, pdip);

%summarize(dip);

% Inputs
y = dip(pdip+1:end);
x = diplags(pdip+1:end, :);
lambda = 0;
tmin = 120; % Holdout sample
const = 1;
roll = 0; % Rolling window, expanding == 0
%

lassoforc = forcerrorsLasso(y, x, lambda, tmin, const, roll);
yh1 = lassoforc.yh;
ff1 = lassoforc.ff;
fe1 = lassoforc.fe;
nfe = length(fe1);

figure
plot(dates(tmin+1:tmin+nfe), yh1)
hold on
plot(dates(tmin+1:tmin+nfe), ff1)

%plot(dates(tmin+1:tmin+nfe), fe1)

rmse(yh1, ff1)
std(yh1)

lambdavec = linspace(0, 2, 100);
[lambdaopt, ~, msevec] = forcmseLambda(y, x, lambdavec, tmin, const, roll);

plot(msevec(1:100))

lassoopt = forcerrorsLasso(y, x, lambdaopt, tmin, const, roll);
figure
plot(dates(tmin+1:tmin+nfe), lassoopt.yh)
hold on
plot(dates(tmin+1:tmin+nfe), lassoopt.ff)

