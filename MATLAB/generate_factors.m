% -------------------------------------------------------------------------
% Generate factors
%
%
%
%
% -------------------------------------------------------------------------

%%%%
% Import data

% JLN 2015
clear; clc;
load jlndata; 
ind         = 132+(6:15); % "duplicate" series to remove    
data(:,ind) = []; 
names(ind)  = [];
x           = data;

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

% Find optimal number of factors according to Bai & Ng (2002)

kmax = 10; % Max number of factors to be extracted
gnum = 2; % ICp2 chosen in JLN2015
demean = 2; % Standardise data

bnicv = zeros(kmax,1);
for k = 1:kmax
    bnicv(k) = bnic(x, k, gnum, demean); % Compute BNIC for each k
end

bnicmin = min(bnicv);
rhat = minind(bnicv); % Optimal number of factors according to lowest IC
fprintf('\nFactors via IC: rhat = %d \n', rhat);

% Extract factors via PCA
[ehat1,Fhat1,lamhat1,ev1]  = jln_factors(x,kmax,gnum,demean);
[Fhat, Lhat, ehat, ev] = factors(x, rhat, demean);

sumeigval = cumsum(ev)/sum(ev);
R2_static = sum(ev(1:rhat))/sum(ev);

% Save factors for further analysis


