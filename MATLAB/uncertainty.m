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
%[evarf, phif] = compute_uf(sf, tf, fb, h);



% Initialisation
h  = 12;

T  = size(sy, 1);

bf = sparse(fbetas);
by = sparse(ybetas);

[numparf, R] = size(bf);
[numpary, N] = size(by);



% Compute uf



tf(1, :) = tf(1, :).* (1 - tf(2, :)); % Reparameterise mean



% Compute expected h-step-ahead variance in factors


% Compute expected variance E[sigma2_F(t+h)]
evarf1 = cell(h,1);

for j = 1:h
    for i = 1:R
        alpha       = tf(1,i);
        beta        = tf(2,i);
        tau         = tf(3,i);
        x           = sf(:,i);
        evarf1{j}(:,i) = expectvar(x, alpha, beta, tau, j);
    end
end


% evarf2
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

isequal(evarf(:, 1, 4), evarf1{4}(:, 1))


% Build VAR representation of the system of factors

fvar = bf(1, :)'; % Collect intercepts
for i = 2:numparf % Append VAR parameter matrices
    
    fvar = [fvar, diag(bf(i, :))]; % Diag because factors are not cross-correlated
end

% Build parameter matrix Phi of companion form
phif = companion(R, pf, fvar);



%%%%
% Compute uncertainty in macro variables




phiy1 = companion(1, py, by(1:py+1, 1)');
lambda1 = [[by(py+2:end, 1); zeros(R*(pf-pz),1)], zeros(R*py, py-1)]';
lambda2(1, 1:numpary-(py+1)) = by(py+2:end, 1);


lambda2 = sparse(1, 1:numpary-(py+1), by(py+2:end, 1), py, R*pf);


ymodel1 = ymodels(6:end, 1);
ybtest1 = ybetas(6:end, 1);
lambdatest = full(lambda);




% Initialisation
evary  = zeros(T, N, h);
ut     = zeros(T, N, h);

phiy   = sparse(py, py, 0);
phi    = sparse(R*pf + py, R*pf + py, 0);

for i = 1:N
    
    
    i = 1
    
    
    % Build parameter matrix for variable's FAVAR representation
    lambda = sparse(1, 1:numpary-(py+1), by(py+2:end, i), py, R*pf); % Lambda contains parameters from reg(y ~ z) on first row, then filled up with zeros to match dims
    phiy   = companion(1, py, by(1:py+1, i)');
    
    phi    = [phif, zeros(R*py, py); lambda, phiy];
    
    
    % Compute uncertainty using the recursion
    alpha       = ty(1, i);
    beta        = ty(2, i);
    tau         = ty(3, i);
    x           = sy(:, i);
    
    % Compute h-step-ahead expected variance
    for j = 1:h
        evary(:, i, j) = expectvar(x, alpha, beta, tau, j); 
    end
    
    
    for t = 1:T
       for j = 1:h
           
           t = 600
           j = 3
           
           % Diagonal matrix of h-step-ahead variance forecasts for each variable/factor separately
           % Zeros on diag where depvar vector [Zt, Yjt]' contains lagged terms
           evh = diag([evarf(t, :, j)'; zeros(R*pf - R, 1); evary(t, i, j)'; zeros(py - 1, 1)]);
           
           evh = sparse(1:R*pf + py, 1:R*pf + py, [evarf(t, :, j)'; zeros(R*pf - R, 1); evary(t, i, j)'; zeros(py - 1, 1)]);
           
           if j == 1
               ui = evh;
           else
               ui = phi* ui* phi' + evh;
           end
           
           ut(t, i, j) = ui(R*pf + 1, R*pf + 1); % Select element corresponding to yt
        
       end
    end
    
    
    
    if pf >1; evf0 = sparse(1,R*(pf-1),0); end;
    if pf==1; evf0 = []; end;
    if py >1; evy0 = sparse(1,py-1,0); end;
    if py==1; evy0 = []; end;
    
    
    % Compute uncertainty in series i
    for t = 1:T
        for j = 1:h
            
            ev1 = sparse(1:R*pf+py,1:R*pf+py,[evarf1{j}(t,:),evf0,evary(t, i, h),evy0]);
            
            if j == 1
                u = ev;
            end
            
            if j  > 1
                u = phi*u*phi' + ev;
            end
            
            U(t,j) = u(r*pf+1,r*pf+1); % select relevant entry
        end
    end
    
    
    
    
    
end


fullphi = full(phi);


%%%%
% Compute uy



% Compute uncertainty using the recursion
alpha = thy(1);
beta  = thy(2);
tau2  = thy(3);
x     = xy;
for j = 1:h
    evy{j} = expectedvar(alpha,beta,tau2,x,j); % Et[(v^y_t)^2]
end;
for t = 1:T
    for j = 1:h
        ev = sparse(1:r*pf+py,1:r*pf+py,[evf{j}(t,:),evf0,evy{j}(t),evy0]);
        if j == 1; u = ev; end;
        if j  > 1; u = phi*u*phi' + ev; end; 
        U(t,j) = u(r*pf+1,r*pf+1); % select relevant entry    
    end 
end



%%%%




for i = 1:N
    tic;
    yb    = sparse(ybetas(i,:));
    thy   = [ty(1, i); ty(2, i); ty(3, i).^2];
    xy    = sy;
    ut(:, i, :) = compute_uy(xy, thy, yb, py, evarf, phif);
    fprintf('Series %d, Elapsed Time = %0.4f \n', i, toc);
end






%%%%%%
figure
for i = [1, 3, 12]
    plot(dates, ut(:, 1, i))
    hold on
end
legend('show')



%%%%
svy = csvread('svymeans.csv', 1);


xy = svy(:, 4:621)'; 
thy = [svy(:, 1), svy(:, 2), svy(:, 3)]'; % Parameter estimators






%%%%





gy = svy(end - 3 + 1:end, :);
save ut dates ut
save geweke dates gy gf names


