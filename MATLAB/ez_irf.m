% -------------------------------------------------------------------------
% Estimate impulse response functions
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB;..\..\..\Data')

% ----------------
% Load data
%h5disp('ez_vardata.h5')
vdates = datetime(h5read('ez_vardata.h5', '/dates'));
vnames = h5read('ez_vardata.h5', '/varnames');
vdata  = h5read('ez_vardata.h5', '/data')';

% MP shocks for euro area
load ez_shocks

% Match sample of ez_vardata to external shocks
bsample = matchsample(vdates, babDates);
nsample = matchsample(vdates, neuDates);
msample = matchsample(vdates, mdates);



% ----------------
% IRF via R&R incomplete VAR regression
sample = msample;
shocks = mshocks;

% Select dependent variable IP
depvar = diff(vdata(:, strcmp(vnames, {'EKIPTOT.G'})));

% Monthly dummies
mdum = (month(vdates) == 12);

% Lag order
const = 1;
dlag  = 11; 
ylag  = 12;
slag  = 24;
maxlag = max([dlag, ylag, slag]);

% Lag before sample select, such that presample counts
mdumlags  = mlag(mdum, dlag);
depvarlags = mlag(depvar, ylag);
shocklags  = mlag(shocks, slag);

% Matrix of independent variables
X = [ones(length(sample), 1), mdumlags(sample, :), depvarlags(sample, :), shocklags(:, :)];

% OLS
linreg = ols(depvar(sample), X);

% The estimated betas of the shock are the same as in paper. Yay!
c = linreg.beta((const + dlag + ylag + 1):end);

% Append zeros to autoregressive coeffs
b = [linreg.beta((const + dlag + 1):((const + dlag + ylag))); zeros(maxlag - ylag, 1)];

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

hmax = 36;

impact = zeros(hmax, 1);
for h = 1:hmax
    AH = A^h;
    impact(h) = AH(varselect, shockselect); % Impact multipliers
end

varirf = cumsum(impact);

plot(varirf)


% ----------------
% IRF via Jorda's local projection method

% Choose variable of interest from factor model dataset
selectVariables = {'EKIPMAN.G', 'EKIPTOT.G', 'EKCPHARMF', 'EMCPCOR5F', 'EKCPCOREF'};
depvar = diff(vdata(:, strcmp(vnames, {'EKIPTOT.G'})));

%aroptlag(depvar, 25, 'aic', 1, 0, 1);

% Construct leads and lags of dependent variable
hmax = 36;
ylag = 3;

yLeads = lead(depvar, hmax);
yLags  = mlag(depvar, ylag);

% Build set of controls
dlag  = 11;
ulag  = 12;

mdumlags = mlag(mdum, dlag);
%udumlags = mlag(udum(:, 1), ulag);

controls = [ones(length(vdates), 1), mdumlags];

% Choose lag order for NW correction
%nlag  = fix(4*(T/100)^(2/9)); % Rule-of-thumb (N&W 1994)
nlag = 4;

% Estimate and plot IRF
mjorda = irf_jorda(yLeads(msample(1:end-hmax), :), yLags(msample(1:end-hmax), :), mshocks(1:end-hmax), controls(msample(1:end-hmax), :), nlag);
%wjorda = irf_jorda(yLeads(wsample, :), yLags(wsample, :), wxres, controls(wsample, :), nlag);
bjorda = irf_jorda(yLeads(bsample(1:end-hmax), :), yLags(bsample(1:end-hmax), :), babShocks(1:end-hmax), controls(bsample(1:end-hmax), :), nlag);
njorda = irf_jorda(yLeads(nsample(1:end-hmax), :), yLags(nsample(1:end-hmax), :), neuShocks(1:end-hmax), controls(nsample(1:end-hmax), :), nlag);

plot(mjorda.irf)
hold on
plot(bjorda.irf)
hold on
plot(njorda.irf)
hold on
plot(wjorda.irf)
hold on
refline(0, 0)



% Desired interval coverage via HAC standard errors
coverage = 0.95;
quantile = norm_inv(coverage);

thetaup = mjorda.theta + quantile*mjorda.sigma;
thetalo = mjorda.theta - quantile*mjorda.sigma;

irfup = cumsum(thetaup);
irflo = cumsum(thetalo);

plot(mjorda.irf)
hold on
plot(irfup)
hold on
plot(irflo)
hold on
refline(0, 0)



% Interaction with uncertainty
% ----------------
% Aggregate macroeconomic uncertainty
load ez_uncertainty

% Match uncertainty to shocks
musample = matchsample(dates, mdates);


% Convert uncertainty to High-U-dummies via quantiles
Umean = mean(Ufac);
Ustd = std(Ufac);
udum = (Ufac >= Umean + 1*Ustd);

summarize(Ufac(:, 1));
summarize(Ufac(musample, 1));

tUmean = (mean(Ufac) - mean(Ufac(musample, :)))./Ustd;
% ----------------
% Center uncertainty for interaction such that average U in period
% corresponds to no conditional effect of U on Y. Take full sample Umean as
% difference between means is minimal and far from significance
centeredUfac = Ufac - Umean;
plot(centeredUfac(:, 1))

interaction = mshocks.*(centeredUfac(musample, 1));
plot(mshocks)
hold on
plot(interaction)

mintshocks = [mshocks, centeredUfac(musample, 1), interaction];

% Estimate IRF with interaction
mintjorda = irf_jorda(yLeads(msample(1:end-hmax), :), yLags(msample(1:end-hmax), :), mintshocks(1:end-hmax, :), controls(msample(1:end-hmax), :), nlag);

% Plot interaction IRF
figure
for i = 1:3
    subplot(3, 1, i)
    plot(mintjorda.irf(:, i))
end


% ----------------
% Standardise uncertainty for interaction 
standardUfac = standardize(Ufac);

plot(standardUfac(:, 1))

interaction = mshocks.*(standardUfac(musample, 1));
plot(mshocks)
hold on
plot(interaction)

mintshocks = [mshocks, standardUfac(musample, 1), interaction];

% Estimate IRF with interaction
mintjorda = irf_jorda(yLeads(msample(1:end-hmax), :), yLags(msample(1:end-hmax), :), mintshocks(1:end-hmax, :), controls(msample(1:end-hmax), :), nlag);

% Plot interaction IRF
figure
for i = 1:3
    subplot(3, 1, i)
    plot(mintjorda.irf(:, i))
end

table(mintjorda.theta)







