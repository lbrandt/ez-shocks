% -------------------------------------------------------------------------
% Generate forecast via factor model
%
%
%
%
% -------------------------------------------------------------------------

% Load and manipulate data in import_data. Call results here:
import_data


% Simulate data as in testnfac.m by Serena Ng
clear; clc;
rng(999);

r = 4; 
N = 100;
T = 50;

e = randn(T,N);
f = randn(T,r);
lambda = randn(N,r);
x = f*lambda'+ e;

fprintf('Data generating process: r = %d \n', r);
fprintf('Sample: T = %d, N = %d \n', size(x));

%%%%
% Find optimal number of factors according to Bai & Ng (2002)
kmax = 20; % Max number of factors to be extracted
gnum = 2; % ICp2 chosen in JLN2015
demean = 2; % Standardise data

bnicv = zeros(kmax,1);
for k = 1:kmax
    bnicv(k) = bnic(x, k, gnum, demean); % Compute BNIC for each k
end

bnicmin = min(bnicv);
rhat = minind(bnicv); % Optimal number of factors according to lowest IC
fprintf('\nFactors via IC(%d): rhat = %d \n', gnum, rhat);


%%%%
% Extract factors via PCA

%[ehat1, Fhat1, lamhat1, ev1]  = jln_factors(x, kmax, gnum, demean);
%[ehat1, Ghat1, lamhat1, ev1]  = jln_factors(x.^2, kmax, gnum, demean);

[Fhat, LF, ef, evf] = factors(x, rhat, demean);
[Ghat, LG, eg, evg] = factors(x.^2, rhat, demean);

sumeigval = cumsum(evf)/sum(evf);
R2_static = sum(evf(1:rhat))/sum(evf);





