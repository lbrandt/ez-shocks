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
numfac = minind(bnicv); % Optimal number of factors according to lowest IC

% Extract factors via PCA

% Save factors for further analysis

