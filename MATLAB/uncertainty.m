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

% Compute uf
h        = 12;
bf       = sparse(fbetas);
tf(1, :) = tf(1, :).* (1 - tf(2, :)); % Reparameterise mean

% Compute expected h-step-ahead variance in factors
i=1;

alpha       = tf(1,i);
beta        = tf(2,i);
tau        = tf(3,i);

evarf = expectvar(sf(:, i), alpha, beta, tau, 1);


% Compute evf
evarf = cell(h,1);

for j = 1:h
    for i = 1:R
        alpha       = tf(1,i);
        beta        = tf(2,i);
        tau         = tf(3,i);
        x           = sf(:,i);
        evarf{j}(:,i) = expectvar(x, alpha, beta, tau, j); %Et[(v^f_t)^2]
    end
end



% Build VAR representation of the system of factors
R   = size(bf, 2);

fvar = bf(1, :)'; % Collect intercepts
for i = 2:size(bf, 1) % Append VAR parameter matrices
    
    fvar = [fvar, diag(bf(i, :))]; % Diag because factors are not cross-correlated
end

% Build parameter matrix Phi of companion form
phif = companion(R, pf, fvar);


%%%%
% Compute uncertainty in macro variables
[T, N] = size(vyt);
ut     = zeros(T, N, h);

yb    = sparse(ybetas);

%%%%
% Compute uy

% Initialize parameters

h  = length(evarf);
r  = size(evarf{1},2);
pf = size(phif,1)/r;
pz = (length(yb)-1-py)/r;
T  = length(xy);

% Preallocate variables
U = zeros(T,h);
if pf >1; evf0 = sparse(1,r*(pf-1),0); end;
if pf==1; evf0 = []; end;
if py >1; evy0 = sparse(1,py-1,0); end;
if py==1; evy0 = []; end;

% Construct the main phi matrix
if pf >pz; lambda_topright = sparse(1,(pf-pz)*r,0); end;
if pf==pz; lambda_topright = []; end;
if py >1;  lambda_bottom   = sparse(py-1,r*pf,0); end;
if py==1;  lambda_bottom   = []; end;
lambda   = [yb(py+2:end),lambda_topright;lambda_bottom];

phiy_top = yb(2:py+1);
if py >1; phiy_bottom = [sparse(1:py-1,1:py-1,1),sparse(py-1,1,0)]; end;
if py==1; phiy_bottom = []; end;
phiy         = [phiy_top;phiy_bottom];
phi_topright = sparse(r*pf,py,0);
phi          = [phif,phi_topright;lambda,phiy];






%%%%




for i = 1:N
    tic;
    yb    = sparse(ybetas(i,:));
    thy   = [svy(1, i).*(1 - svy(2, i)); svy(2, i); svy(3, i).^2];
    xy    = svy(4:end - 3, i);
    ut(:, i, :) = compute_uy(xy, thy, yb, py, evf, phif);
    fprintf('Series %d, Elapsed Time = %0.4f \n', i, toc);
end






%%%%%%
figure
plot(dates, [vyt(:, 1), (svy(1, 1:618))'])
legend('show')


%%%%
svy = csvread('svymeans.csv', 1);


xy = svy(:, 4:621)'; 
thy = [svy(:, 1), svy(:, 2), svy(:, 3)]'; % Parameter estimators






%%%%





gy = svy(end - 3 + 1:end, :);
save ut dates ut
save geweke dates gy gf names


