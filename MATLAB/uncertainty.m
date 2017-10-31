% -------------------------------------------------------------------------
% Compute matrix of macro uncertainty estimates for horizons 1 through 12
% -------------------------------------------------------------------------

%clear; clc;

% Load data
load factors_forc;

sf = csvread('svflatent.csv', 1); % Estimates of latent process s(t)
tf = csvread('svfparams.csv', 1); % AR(1) parameter estimators

sy = csvread('svylatent.csv', 1);
ty = csvread('svyparams.csv', 1);


% Compute uf
h        = 12;
bf       = sparse(fbetas);
tf(1, :) = tf(1, :).* (1 - tf(2, :)); % Reparameterise mean



%[evarf, phif] = compute_uf(sf, tf, fb, h);


% Initialize parameters
R   = size(bf, 2);

% Create parameter matrix phif
phif_top = [];
for j = 2:pf+1
    phif_top = [phif_top,sparse(1:R,1:R,fb(:,j))]; 
end

if pf > 1
    phif_bot = [sparse(1:R*(pf-1),1:R*(pf-1),1),sparse(R*(pf-1),R,0)];
    phif     = [phif_top;phif_bot];
else
    phif = phif_top;
end




%%%%
svy = csvread('svymeans.csv', 1);


xy = svy(:, 4:621)'; 
thy = [svy(:, 1), svy(:, 2), svy(:, 3)]'; % Parameter estimators


figure
plot(dates, [vyt(:, 1), (svy(1, 1:618))'])
legend('show')


% Compute objects from predictors
h      = 12;
fb     = sparse(fbetas);
thf    = [svf(1, :).*(1-svf(2,:));svf(2,:);svf(3,:).^2];
xf     = svf(4:end-3,:);
gf     = svf(end-3+1:end,:);
[evf,phif] = compute_uf(xf,thf,fb,h);




% Compute uncertainty
[T, N] = size(vyt);
ut     = zeros(T, N, h);

for i = 1:N
    tic;
    yb    = sparse(ybetas(i,:));
    thy   = [svy(1, i).*(1 - svy(2, i)); svy(2, i); svy(3, i).^2];
    xy    = svy(4:end - 3, i);
    ut(:, i, :) = compute_uy(xy, thy, yb, py, evf, phif);
    fprintf('Series %d, Elapsed Time = %0.4f \n', i, toc);
end





gy = svy(end - 3 + 1:end, :);
save ut dates ut
save geweke dates gy gf names


