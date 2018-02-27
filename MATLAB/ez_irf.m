% -------------------------------------------------------------------------
% Estimate impulse response functions
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB;..\..\..\Data')

% ----------------
% Load data
load ez_data
load ez_shocks
load ez_uncertainty

% Match sample of full data set to shorter series
msample = matchsample(dates, mdates);
usample = matchsample(dates, udates);

% Monthly dummies
mdum = (month(dates) == 12);
mdumlags = mlag(mdum, 11);


% Choose variables of interest from factor model dataset
selectVariables = {'EKIPTOTG', 'EMESHARM'};
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
P(1) = 3;
% Lag length for prices is not as clear. Set p = 12 via HQC.
P(2) = 12;

% Standardise shocks such that IRF has the interpretation of the effect of
% a unit standard deviation shock on the original scale
mshocks = standardise(mshocks);

% Maximum IRF horizon. Choose fairly short because Jorda IRF requires
% direct forecasting H steps into the future.
H = 36;


% ----------------
% IRF via incomplete VAR companion form Sims (2012)
Q = 24;

compirfs = zeros(H, N);
for i = 1:N
    mcompanion = irf_companion(depvars(msample, i), mshocks, P(i), Q, H, [ones(length(msample), 1), mdumlags(msample, :)]);
    compirfs(:, i) = mcompanion.irf;
end

figure
for i = 1:N
    subplot(2, 1, i)
    plot(compirfs(:, i))
    hold on
    refline(0, 0)
    title(selectVariables(i))
end




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
    subplot(2, 1, i)
    plot(jordairfs(:, i))
    hold on
    refline(0, 0)
    title(selectVariables(i))
end




% ----------------
% Compare IRF construction methods

figure
for i = 1:N
    subplot(2, 1, i)
    plot(jordairfs(:, i))
    hold on
    plot(compirfs(:, i))
    hold on
    refline(0, 0)
    title(selectVariables(i))
end


plot(mdates, mshocks)
plot(depvars(:, 2))

summarize(depvars);

[maxshock, maxind] = max(mshocks);
mdates(maxind)

[minirf, minind] = min(jordairfs)


% Interaction with uncertainty
% ----------------
% Aggregate macroeconomic uncertainty


% Match uncertainty to shocks
musample = matchsample(udates, mdates);
[tu, hu] = size(ufac);


umean = mean(ufac);
ustd = std(ufac);

summarize(ufac(:, 1));
summarize(ufac(musample, 1));

tUmean = (mean(ufac) - mean(ufac(musample, :)))./ustd;



% ----------------
% Standardise uncertainty for interaction such that average U in period
% corresponds to no conditional effect of U on Y. Also, U_t = 1 means that
% U was one SD above mean on impact.
standardufac = standardize(ufac);


nlag = 4;
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

% Plot interaction IRF at U_t = 1, i.e. one standard deviation
mintsum = mintjorda.irf(:, 1) + mintjorda.irf(:, 3);

plot(mintsum)

% Plot decomposed IRF
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



