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

phimat = cell(N,1);

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
               %ui = evh/udiffmean(i);
           else
               ui = phi* ui* phi' + evh; % for h>1, Omega = phi* Et[sigma2_Y(t+h-1)]* phi + Et[sigma2_Y(t+h)]
           end
           
           phimat{i} = phi;
           
           ut(t, i, j) = ui(R*pf + 1, R*pf + 1); % Select element corresponding to yt, this is the squared uncertainty in this variable at time t looking h-steps ahead (?)
       end
    end
    
    fprintf('Series %d, Elapsed Time = %0.4f \n', i, toc);
    
end

Uind = sqrt(ut);



%%%%
% Aggregate individual uncertainty
% Simple average
Uavg = squeeze(mean(Uind,2));


% Principal component analysis
% Transform individual uncertainties to logdiffs
logu  = log(Uind(:, :, :));
dlogu = logu(2:end, :, :) - logu(1:end-1, :, :);

dufac = zeros(T-1, N, h);

Upca    = zeros(T, h);
Upcascaled = zeros(T, h);

Upcalog = zeros(T, h);
Upcadlog = zeros(T-1, h);

for i = 1:h
    
    Ufac(:, i)    = factors(Uind(:, :, i), 1, 2);
    Upcalog(:, i) = factors(logu(:, :, i), 1, 2);
    Upcadlog(:, i) = factors(dlogu(:, :, i), 1, 2);
    
    % Flip Ufac if necessary
    rho(i) = corr(Ufac(:, i), Uavg(:, i));
    if rho < 0
        Ufac = -Ufac;
    end
    
    % Scale to Uavg  
    % Which variance?
    Upca(:, i) = standardise(Upca(:, i))* sqrt(var(Uavg(:, i))) + mean(Uavg(:, i));
end

dUpcalog = Upcalog(2:end, :) - Upcalog(1:end-1, :);







for j = 1:h
   logu = log(Uind(:, :, j));
   dlogu(:, :, j)  = logu(2:end, :) - logu(1:end-1, :);
   
   
   [de, du, dl, dv] = jln_factors(dlogu(:, :, j), 5, 2, 1);
   
   % Rotate estimate
   rho     = corr(cumsum([0; du(:, 1)]), Uavg(:, j));
   
   if rho < 0
       du = -du; 
       dl = -dl; 
   end
   
   ufac    = cumsum([zeros(1, size(du, 2)); du]);
   dufac(:, :, j) = du;
   dlam(:, :, j)  = dl;
   deig(:, j)  = dv;
   
   % Calibrate to cross-section mean
   sd      = std(Uavg(:, j));
   mn      = mean(Uavg(:, j));
   p0      = [1, 0.5];
   opt     = optimset('tolfun', 1e-50, 'display','off');
   [p, obj] = fminsearch(@(p)calibratef(p, ufac(:, 1), sd, mn), p0, opt);
   Upca1(:, j) = exp((p(1)* ufac(:, 1) + p(2))./ 2); 
end














%%%%
% Plot

% Aggregate uncertainty, simple average
figure
for i = [1, 3, 12]
    plot(dates, Uavg(:, i))
    hold on
end
legend('show')


% Single series uncertainty
hselect = 12;

figure
for i = 1:10
    plot(dates, Uind(:, i, hselect))
    hold on
end



%%%%
% Save results




gy = svy(end - 3 + 1:end, :);
save ut dates ut
save geweke dates gy gf names


