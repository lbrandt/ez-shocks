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

% Match sample of ez_vardata to shocks from GJ database
bsample = matchsample(vdates, babDates);
nsample = matchsample(vdates, neuDates);



% Aggregate macroeconomic uncertainty
load ez_uncertainty

% Convert uncertainty to High-U-dummies via quantiles
Umean = mean(Ufac);
Ustd = std(Ufac);
udum = (Ufac >= Umean + 1*Ustd);



% ----------------
% IRF via R&R incomplete VAR regression
sample = bsample;
shocks = babShocks;

% Select dependent variable IP
depvar = diff(vdata(:, strcmp(vnames, {'EKIPMAN.G'})));

% Monthly dummies
mdum = (month(vdates) == 12);

% Lag order
const = 1;
dlag  = 11; 
ylag  = 12;
slag  = 36;
maxlag = max([dlag, ylag, slag]);

% Lag before sample select, such that presample counts
mdumlags  = mlag(mdum, dlag);
depvarlags = mlag(dip, ylag);
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

hmax = 48;

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
depvar = diff(vdata(:, strcmp(vnames, {'EKIPMAN.G'})));

% Construct leads and lags of dependent variable
hmax = 48;
ylag = 12;

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
%wjorda = irf_jorda(yLeads(wsample, :), yLags(wsample, :), wxres, controls(wsample, :), nlag);
bjorda = irf_jorda(yLeads(bsample, :), yLags(bsample, :), babShocks, controls(bsample, :), nlag);
njorda = irf_jorda(yLeads(nsample, :), yLags(nsample, :), neuShocks, controls(nsample, :), nlag);

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

thetaup = bjorda.theta + quantile*bjorda.sigma;
thetalo = bjorda.theta - quantile*bjorda.sigma;

irfup = cumsum(thetaup);
irflo = cumsum(thetalo);

plot(bjorda.irf)
hold on
plot(irfup)
hold on
plot(irflo)
hold on
refline(0, 0)


interaction = neuShocks.*(Ufac(sample, 1)/Ustd(1));
plot(neuShocks)
hold on
plot(interaction)



%%%% LIKE THIS?!?
% Reg BIP ~ shocks

% Reg BIP ~ shocks, uncertainty, shocks*uncertainty

% Reg BIP ~ shocks, uncertainty, shocks*Udum
