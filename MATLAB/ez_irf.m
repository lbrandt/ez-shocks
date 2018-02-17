% -------------------------------------------------------------------------
% Estimate impulse response functions
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB;..\..\..\Data')

% ----------------
% Load data
load ez_data
load ez_shocks
load eu_uncertainty

% Match sample of full data set to shorter series
msample = matchsample(dates, mdates);
usample = matchsample(dates, udates);

bsample = matchsample(dates, babDates);
nsample1 = matchsample(dates, neuDates);
nsample2 = matchsample(neuDates, dates);



% Choose variables of interest from factor model dataset
selectVariables = {'EKIPMANG', 'EKIPTOTG', 'EKCPHARMF', 'EMESHARM'};
depvars = x(:, findstrings(names, selectVariables));

N = length(selectVariables);

% AR order of depvars
pmax = 36;
const = 1;
trend = 0;

picx = zeros(N, 3);
icx = zeros(N, 3);
for i = 1:N
    [picx(i, 1), icx(i, 1)] = aroptlag(depvars(:, i), pmax, 'aic', const, trend, 0);
    [picx(i, 2), icx(i, 2)] = aroptlag(depvars(:, i), pmax, 'bic', const, trend, 0);
    [picx(i, 3), icx(i, 3)] = aroptlag(depvars(:, i), pmax, 'hqc', const, trend, 0);
end

P = zeros(1, N);
% Lag length for economic activity unanimously suggested is p = 3.
P(1:2) = 3;
% Lag length for prices is not as clear. Set p = 12.
P(3:4) = 12;




% Monthly dummies
mdum = (month(dates) == 12);
mdumlags = mlag(mdum, 11);

% ----------------
% IRF via incomplete VAR companion form Sims (2012)
%P = 12;
Q = 24;
H = 36;

compirfs = zeros(H, N);
for i = 1:N
    mcompanion = irf_companion(depvars(msample, i), mshocks, P(i), Q, H, [ones(length(msample), 1), mdumlags(msample, :)]);
    compirfs(:, i) = mcompanion.irf;
end

figure
for i = 1:N
    subplot(2, 2, i)
    plot(compirfs(:, i))
    hold on
    refline(0, 0)
    title(selectVariables(i))
end


% Test shocks EKIPTOTG only
mcompanion = irf_companion(depvars(msample, 2), mshocks, P(2), Q, H, [ones(length(msample), 1), mdumlags(msample, :)]);
bcompanion = irf_companion(depvars(bsample, 2), babShocks, P(2), Q, H, [ones(length(bsample), 1), mdumlags(bsample, :)]);
ncompanion = irf_companion(depvars(nsample1, 2), neuShocks(nsample2), P(2), Q, H, [ones(length(nsample1), 1), mdumlags(nsample1, :)]);

plot(mcompanion.irf)
hold on
plot(bcompanion.irf)
hold on
plot(ncompanion.irf)




% ----------------
% IRF via Jorda's local projection method
controls = [ones(length(dates), 1), mdumlags];

jordairfs = zeros(H, N);
for i = 1:N
    % Construct leads and lags of dependent variable
    yLeads = lead(depvars(:, i), H);
    yLags  = mlag(depvars(:, i), P(i));
    
    mjorda = irf_jorda(yLeads(msample(1:end-H), :), yLags(msample(1:end-H), :), mshocks(1:end-H), controls(msample(1:end-H), :));
    jordairfs(:, i) = mjorda.irf;
end

figure
for i = 1:N
    subplot(2, 2, i)
    plot(jordairfs(:, i))
    hold on
    refline(0, 0)
    title(selectVariables(i))
end


% Test Jorda EKIPTOTG only
ylag = 12;

% Construct leads and lags of dependent variable
yLeads = lead(depvars(:, 2), H);
yLags  = mlag(depvars(:, 2), ylag);


% Estimate and plot IRF
mjorda = irf_jorda(yLeads(msample(1:end-H), :), yLags(msample(1:end-H), :), mshocks(1:end-H), controls(msample(1:end-H), :));
bjorda = irf_jorda(yLeads(bsample(1:end-H), :), yLags(bsample(1:end-H), :), babShocks(1:end-H), controls(bsample(1:end-H), :));
njorda = irf_jorda(yLeads(nsample1(1:end-H), :), yLags(nsample1(1:end-H), :), neuShocks(nsample2(1:end-H)), controls(nsample1(1:end-H), :));

plot(mjorda.irf)
hold on
plot(bjorda.irf)
hold on
plot(njorda.irf)


% Interval via HAC standard errors
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



% ----------------
% Compare IRF construction methods

figure
for i = 1:N
    subplot(2, 2, i)
    plot(jordairfs(:, i))
    hold on
    plot(compirfs(:, i))
    hold on
    refline(0, 0)
    title(selectVariables(i))
end





% Interaction with uncertainty
% ----------------
% Aggregate macroeconomic uncertainty
load ez_uncertainty

% Match uncertainty to shocks
musample = matchsample(udates, mdates);
[tu, hu] = size(ufac);


umean = mean(ufac);
ustd = std(ufac);

summarize(ufac(:, 1));
summarize(ufac(musample, 1));

tUmean = (mean(ufac) - mean(ufac(musample, :)))./ustd;
% ----------------
% Center uncertainty for interaction such that average U in period
% corresponds to no conditional effect of U on Y. Take full sample Umean as
% difference between means is minimal and far from significance
centufac = ufac - umean;
plot(centufac(:, 1))

interaction = mshocks.*(centufac(musample, 1));
plot(mshocks)
hold on
plot(interaction)

mintshocks = [mshocks, centufac(musample, 1), interaction];

mintjorda = irf_jorda(yLeads(msample(1:end-H), :), yLags(msample(1:end-H), :), mintshocks(1:end-H, :), controls(msample(1:end-H), :), nlag);

% Plot interaction IRF
figure
for i = 1:3
    subplot(3, 1, i)
    plot(mintjorda.irf(:, i))
end


% ----------------
% Standardise uncertainty for interaction 
standardufac = standardize(ufac);
%plot(udates, standardufac(:, 1))

mintparm = zeros(H, N, hu);
mintparu = zeros(H, N, hu);
for i = 1:N
    
    for j = 1:hu
        
        yLeads = lead(depvars(:, i), H);
        yLags  = mlag(depvars(:, i), P(i));

        interaction = mshocks.*standardufac(musample, j);
        mintshocks = [mshocks, standardufac(musample, j), interaction];

        mintjorda = irf_jorda(yLeads(msample(1:end-H), :), yLags(msample(1:end-H), :), mintshocks(1:end-H, :), controls(msample(1:end-H), :), nlag);
        
        mintparm(:, i, j) = mintjorda.theta(:, 1); % Partial effect parameters
        mintparu(:, i, j) = mintjorda.theta(:, 3); % Interaction parameters
    end
end

mint1sigmaup = mintparm + mintparu;
mint1sigmalo = mintparm - mintparu;



% Test IRF with interaction on EKIPTOTG and ufac at h = 1
ylag = 12;

yLeads = lead(depvars(:, 2), H);
yLags  = mlag(depvars(:, 2), ylag);

interaction = mshocks.*(standardufac(musample, 1));
mintshocks = [mshocks, standardufac(musample, 1), interaction];

mintjorda = irf_jorda(yLeads(msample(1:end-H), :), yLags(msample(1:end-H), :), mintshocks(1:end-H, :), controls(msample(1:end-H), :), nlag);

% Plot interaction IRF
figure
for i = 1:3
    subplot(3, 1, i)
    plot(mintjorda.irf(:, i))
end

table(mintjorda.theta)

mintirfup = mintjorda.irf(:, 1) + mintjorda.irf(:, 3);
mintirflo = mintjorda.irf(:, 1) - mintjorda.irf(:, 3);

plot(mintjorda.irf(:, 1))
hold on
plot(mintirfup)
hold on
plot(mintirflo)



