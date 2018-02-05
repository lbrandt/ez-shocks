% -------------------------------------------------------------------------
% Replicate Babecka Kucharcukova et al. (2016)
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB;..\..\..\Data')

% ----------------
% Load data

%h5disp('bab_data.h5')
bdates = datetime(h5read('bab_data.h5', '/dates'));
bnames = h5read('bab_data.h5', '/varnames');
bdata  = h5read('bab_data.h5', '/data')';

bta = datetime('2000-04-01');
bte = datetime('2017-12-01');
btaindex = find(bdates == bta);
bteindex = find(bdates == bte);

bsample = btaindex:bteindex;

% Build dataset
bx = bdata(:, [2:5, 1, 12:19, 21]); % Order data like in paper

bx(:, 10) = []; % Delete series with lots of zeros

% PCA
rhat = 3;
demean = 2;

[Fhat, Lhat, Ehat, eigval] = factors(bx(bsample, :), rhat, demean);

sumeigval = cumsum(eigval)/sum(eigval);
R2_static = sum(eigval(1:rhat))/sum(eigval);

summarize(Fhat);
plot(bdates(bsample), Fhat)

sign = [-1, -1, -1];
Fhat1 = sign.*standardise(Fhat);
Lhat1 = sign.*Lhat;

Chat1 = Fhat1*Lhat1';

figure
for i = 1:14
    subplot(5, 3, i)
    plot(bdates(bsample), standardise(Chat1(:, i)))
    hold on
    plot(bdates(bsample), standardise(bx(bsample, i)))
end

babLoadings = [0.97,  0.09,  0.16; ...
     0.98,  0.08,  0.12; ...
     0.97,  0.03,  0.19; ...
     0.72,  0.07,  0.40; ...
     0.96,  0.17,  0.09; ...
     0.07, -0.16,  0.58; ...
    -0.76,  0.01,  0.57; ...
    -0.83, -0.02,  0.46; ...
    -0.26,  0.93,  0.06; ...
     0.15,  0.00, -0.14; ...
     0.20,  0.82,  0.08; ...
     0.03,  0.11,  0.62; ...
    -0.09,  0.92, -0.11; ...
     0.20,  0.38, -0.19];

weights = eigval(1:rhat)/sum(eigval);
MCI = mean(bx(:, 3)) + std(bx(:, 3))*standardise(Fhat1*weights/sum(weights));

summarize(MCI);

plot(bdates(bsample), standardise(MCI))
hold on
plot(bdates(bsample), Fhat1(:, 1))
hold on
plot(bdates(bsample), standardize(bx(bsample, 3)))

% Not really what I'm looking for. Just take 3-months EURIBOR as indicator
% of monetary policy in VAR.