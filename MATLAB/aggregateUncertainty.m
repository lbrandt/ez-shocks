function results = aggregateUncertainty(sy, ty, sf, tf, ybeta, py, fbeta, pf, h)
% -------------------------------------------------------------------------
% Generates aggregate uncertainty embodied in a set of prediction errors
% from a factor model together with their latent stochastic volatility.
% Adapted from Jurado, Ludvigson, Ng (2015).
%
%   Input
%       sy          Stochastic volatility in macro variables y [T x N]
%       ty          AR model parameters of sy process [3 x N]
%       sf          Stochastic volatility in factors y [T x R]
%       tf          AR model parameters sf process [3 x R]
%       ybeta       Parameters of prediction models for y [I x N]
%       py          Number of autoregressive lags in ybeta
%       fbeta       Parameters of prediction models for f [J x N]
%       pf          Number of autoregressive lags in fbeta
%       h           Maximum number of forecasting steps
%
%   Output
%       results     Structure containing:
%         results.uavg    Aggregate uncertainty by simple average [T x h]
%         results.ufac    Aggregate uncertainty by PCA weighting [T x h]
%         results.uind    Estimates of individual expected volatilities [T x N x h]
%         results.phic    Cell array of parameter matrices {N}
%
%   Dependencies {source}
%       companion
%       expectvar
%
% -------------------------------------------------------------------------

%%%%
% Initialisation
T = length(sy);

bf = sparse(fbeta);
by = sparse(ybeta);

[numparf, R] = size(bf);
[numpary, N] = size(by);

% Reparameterise means
tf(1, :) = tf(1, :).* (1 - tf(2, :));
ty(1, :) = ty(1, :).* (1 - ty(2, :));


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

evary = zeros(T, N, h);
ut    = zeros(T, N, h);
results.phic = cell(N, 1);
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
           
           results.phic{i} = phi; % save phi matrices
           
           ut(t, i, j) = ui(R*pf + 1, R*pf + 1); % Select element corresponding to yt, this is the squared uncertainty in this variable at time t looking h-steps ahead
       end
    end
    
    fprintf('Series %d, elapsed time = %0.4f \n', i, toc);
    
end
% Individual expected volatilities
results.uind = sqrt(ut);



%%%%
% Aggregate individual uncertainty
% Simple average
results.uavg = squeeze(mean(results.uind,2));

% Principal component analysis on logged variances
logu = log(ut(:, :, :));

results.ufac = zeros(T, h);
upca = zeros(T, h);
for i = 1:h
    upca(:, i) = factors(logu(:, :, i), 1, 2);
    
    % Flip ufac if necessary
    rho = corrcoef(upca(:, i), results.uavg(:, i));
    if rho(2, 1) < 0
        upca = -upca;
    end
    
    % Scale to Uavg
    results.ufac(:, i) = exp( standardise(upca(:, i))* std(log(results.uavg(:, i))) + mean(log(results.uavg(:, i))) );
end

end
