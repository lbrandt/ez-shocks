% -------------------------------------------------------------------------
% Compute matrix of macro uncertainty estimates for horizons 1 through 12
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB')

% Load data
load ez_factors_forc;

% New
%h5disp('ez_svresults.h5')
sf = h5read('ez_svresults.h5', '/sf')'; % Estimates of latent process s(t)
tf = h5read('ez_svresults.h5', '/tf')'; % AR(1) parameter estimators

% From LASSO model selection
sy = h5read('ez_svresults.h5', '/sy')';
ty = h5read('ez_svresults.h5', '/ty')';

% From hard thresholding model selection
syht = h5read('ez_svresults.h5', '/syht')'; 
tyht = h5read('ez_svresults.h5', '/tyht')';



%%%%
% Compute uncertainty in macro variables
hmax = 12;

uLasso = aggregateUncertainty(sy, ty, sf, tf, ybetas, py, fbetas, pf, hmax);
uHT = aggregateUncertainty(syht, tyht, sf, tf, htbetas, py, fbetas, pf, hmax);

% Test stationarity
test = full(uLasso.phic{1});
[V, D, W] = eig(test);

amat = [0.85, 0.3, 0.2; 0.125, 0.65, 0.05; -0.25, -0.3, 0.4];
[V, D, W] = eig(amat);

max(eig(test))

% Compare methods

% Plots
figure
for i = 1:hmax
    subplot(4, 3, i)
    plot(udates, uLasso.uavg(:, i))
    hold on
    plot(udates, uHT.uavg(:, i))
    title(['h = ',num2str(i)])
end

figure
for i = 1:hmax
    subplot(4, 3, i)
    plot(udates, uLasso.ufac(:, i))
    hold on
    plot(udates, uHT.ufac(:, i))
    title(['h = ',num2str(i)])
end

% Means
ubar1 = mean(uLasso.uavg);
ubar2 = mean(uHT.uavg);

ubardiff = ubar1 - ubar2;

% AR params
tf(1, :) = tf(1, :).* (1 - tf(2, :));
ty(1, :) = ty(1, :).* (1 - ty(2, :));
tyht(1, :) = tyht(1, :).* (1 - tyht(2, :));

mean(ty, 2)
mean(tyht, 2)

%%%%
% Save results
selectMethod = uLasso;
uavg = selectMethod.uavg;
ufac = selectMethod.ufac;
uing = selectMethod.uind;
phic = selectMethod.phic;
save ez_uncertainty -v7.3 udates uavg ufac uing phic


% Figure from JLN2015
selectH = [1, 3, 12];

figure
subplot(2,1,1);
for i = selectH
    plot(udates, uavg(:, i), 'DisplayName', ['h = ',num2str(i)])
    hold on
end
legend('show')
title('uavg')

subplot(2,1,2);
for i = selectH
    plot(udates, ufac(:, i), 'DisplayName', ['h = ',num2str(i)])
    hold on
end
legend('show')
title('ufac')



