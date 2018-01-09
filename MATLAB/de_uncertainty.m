% -------------------------------------------------------------------------
% Compute matrix of macro uncertainty estimates for horizons 1 through 12
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB')

% Load data
load gs_factors_forc;

%h5disp('gs_svfresults.h5')
sf = h5read('gs_svfresults.h5', '/h')'; % Estimates of latent process s(t)
tf = h5read('gs_svfresults.h5', '/t')'; % AR(1) parameter estimators

%h5disp('gs_svyresults.h5')
sy = h5read('gs_svyresults.h5', '/h')';
ty = h5read('gs_svyresults.h5', '/t')';


%%%%
% Compute uncertainty in macro variables


% Initialisation
h  = 12;

T  = size(sy, 1);

bf = sparse(fbetas);
by = sparse(ybetas);

[numparf, R] = size(bf);
[numpary, N] = size(by);


tf(1, :) = tf(1, :).* (1 - tf(2, :)); % Reparameterise mean
ty(1, :) = ty(1, :).* (1 - ty(2, :)); % Reparameterise mean



%%%%
% Compute expected volatility in factors


% Expected h-step-ahead variance Et[sigma2_F(t+h)]
evarf = zeros(T, R, h);

for j = 1:h
    for i = 1:R
        alpha       = tf(1,i);
        beta        = tf(2,i);
        tau         = tf(3,i);
        x           = sf(:,i);
        
        evarf(:, i, j) = expectvar(x, alpha, beta, tau, j);
    end
end


% Build VAR representation of the system of factors

fvar = bf(1, :)'; % Collect intercepts
for i = 2:numparf % Append VAR parameter matrices, they are diagonal because factors are not cross-correlated
    
    fvar = [fvar, diag(bf(i, :))]; %#ok<AGROW> % Suppress preallocation error
end

% Build parameter matrix Phi of companion form
phif = companion(R, pf, fvar);


%%%%
% Compute uncertainty in macro variables


% Initialisation
evary  = zeros(T, N, h);
ut     = zeros(T, N, h);

phiy   = sparse(py, py, 0);
phi    = sparse(R*pf + py, R*pf + py, 0);

phisav = cell(N,1);

for i = 1:N
    
    tic;
    
    % Build parameter matrix for variable's FAVAR representation
    lambda = sparse(1, 1:numpary-(py+1), by(py+2:end, i), py, R*pf); % Lambda contains parameters from reg(y ~ z) on first row, then filled up with zeros to match dims
    phiy   = companion(1, py, by(1:py+1, i)');
    
    phi    = [phif, zeros(R*py, py); lambda, phiy];
    
    
    % Compute individual uncertainty via recursion
    alpha       = ty(1, i);
    beta        = ty(2, i);
    tau         = ty(3, i);
    x           = sy(:, i);
    
    
    % Compute h-step-ahead expected variance in variable i Et[sigma2_Y(t+h)]
    for j = 1:h
        evary(:, i, j) = expectvar(x, alpha, beta, tau, j); 
    end
    
    for t = 1:T
       for j = 1:h
           
           % Diagonal matrix of h-step-ahead variance forecasts for each variable/factor separately
           % Zeros on diag where depvar vector [Zt, Yjt]' contains lagged terms         
           evh = sparse(1:R*pf + py, 1:R*pf + py, [evarf(t, :, j)'; zeros(R*pf - R, 1); evary(t, i, j)'; zeros(py - 1, 1)]);
           
           if j == 1
               ui = evh; % for h=1, Omega = Et[sigma2_y(t+h)]
           else
               ui = phi* ui* phi' + evh; % for h>1, Omega = phi* Et[sigma2_Y(t+h-1)]* phi + Et[sigma2_Y(t+h)]
           end
           
           phisav{i} = phi; % save phi matrices
           
           ut(t, i, j) = ui(R*pf + 1, R*pf + 1); % Select element corresponding to yt, this is the squared uncertainty in this variable at time t looking h-steps ahead
       end
    end
    
    fprintf('Series %d, Elapsed Time = %0.4f \n', i, toc);
    
end

% Individual expected volatilities
Uind = sqrt(ut);


%%%%
% Aggregate individual uncertainty
% Simple average
Uavg  = squeeze(mean(Uind,2));


% Principal component analysis on logged variances
logut  = log(ut(:, :, :));

Upca    = zeros(T, h);
Ufac    = zeros(T, h);

for i = 1:h
    
    Upca(:, i) = factors(logut(:, :, i), 1, 2);
    
    % Flip Ufac if necessary
    rho = corr(Upca(:, i), Uavg(:, i));
    if rho < 0
        Upca = -Upca;
    end
    
    % Scale to Uavg
    Ufac(:, i) = exp( standardise(Upca(:, i))* std(log(Uavg(:, i))) + mean(log(Uavg(:, i))) );
end


%%%%
% Save results

save gs_uncertainty dates Uavg Ufac ut phisav
