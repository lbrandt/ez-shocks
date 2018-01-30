% -------------------------------------------------------------------------
% Estimate impulse response functions
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB;..\..\..\Data')

% ----------------
% Load data
load de_data

% Babecka Kucharcukova et al. MP shocks for euro area
load ez_shocks

% Aggregate macroeconomic uncertainty
load de_uncertainty

% Convert uncertainty to High-U-dummies via quantiles
Umean = mean(Ufac);
Ustd = std(Ufac);
udum = (Ufac >= Umean + 1*Ustd);

% Monthly dummies
mdum = (month(dates) == 12);

% Identify overlap in sample
% ta = babDates(1);
% te = babDates(end);
ta = neuDates(1);
te = neuDates(end);

dates = dateshift(dates, 'start', 'month');

taindex = find(dates == ta);
teindex = find(dates == te);

sample = taindex:teindex;
T = length(sample);


% ----------------
% IRF via R&R incomplete VAR regression

% Select dependent variable
deip = x(5:end, strcmp(varnames, {'BDIPTOTG'})); % Cut off leading observations because of lags in factor model
depvar = deip(sample);

% Lag order
const = 1;
dlag  = 11; 
ylag  = 12;
slag  = 36;
maxlag = max([dlag, ylag, slag]);

% Lag before sample select, such that presample counts
dum12lags  = mlag(mdum, dlag);
depvarlags = mlag(deip, ylag);
shocklags  = mlag(neuShocks, slag);

% Matrix of independent variables
X = [dum12lags(sample, :), depvarlags(sample, :), shocklags(:, :)];

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

rrirf = cumsum(impact);




% ----------------
% IRF via Jorda's local projection method

% Choose variable of interest from factor model dataset
chooseVariable = {'BDIPTOTG', 'BDUN_TOTQ', 'BDCONPRCE'};
depvar = x(5:end, strcmp(varnames, {'BDIPTOTG'})); % Cut off leading observations because of lags in factor model in order to match dimensions of U

% Choose series of shocks
shocks = neuShocks;

% Construct leads and lags of dependent variable
hmax = 48;
ylag = 12;

yLeads = lead(depvar, hmax);
yLags  = mlag(depvar, ylag);

% Build set of controls
dlag  = 11;
ulag  = 12;

mdumlags = mlag(mdum, dlag);
udumlags = mlag(udum(:, 1), ulag);

controls = [ones(T, 1), mdumlags(sample, :)];

% Choose lag order for NW correction
nlag  = fix(4*(T/100)^(2/9)); % Rule-of-thumb (N&W 1994)

% Estimate and plot IRF
jorda = irf_jorda(yLeads(sample, :), yLags(sample, :), shocks, controls, nlag);

plot(jorda.irf)
hold on
plot(rrirf)
hold on
refline(0, 0)


% Desired interval coverage via HAC standard errors
coverage = 0.90;
quantile = norm_inv(coverage);

thetaup = jorda.theta + quantile*jorda.sigma;
thetalo = jorda.theta - quantile*jorda.sigma;

irfup = cumsum(thetaup);
irflo = cumsum(thetalo);

plot(jorda.irf)
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
