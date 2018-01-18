% -------------------------------------------------------------------------
% Generate monetary policy shocks for Germany
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB')

% Load data
load de_data % Agent's observable information set
%load ecb_decisions % ECB monetary policy decisions
%load ez_shadow % Euro shadow rates as actual monetary policy stance at ZLB


%%%%%%%%%%%%%%%%
% PCA

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


%%%% Override BNIC
rhat = 12
%%%%


% Extract estimated number of factors

[Fhat, LF, ef, evf] = factors(x, rhat, demean);
[Ghat, LG, eg, evg] = factors(x.^2, rhat, demean);

sumeigval = cumsum(evf)/sum(evf);
R2_static = sum(evf(1:rhat))/sum(evf);


%%%%%%%%%%%%%%%%
% Generate shocks

% Build predictor set with nonlinear terms
zt       = [Fhat, Fhat(:, 1).^2, Ghat(:, 1)];
[~, M]   = size(zt);

% Reg i* on Fhat

% Reg i* on Z

% Calculate one-step-ahead forecast errors

%u = i - ihat

%%%% ARE THESE THE CORRECT SHOCKS ALREADY??

% Save data
%save de_shocks -v7.3
