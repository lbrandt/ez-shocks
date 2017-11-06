% -------------------------------------------------------------------------
% Compute matrix of uncertainty estimates for horizons 1 through 12
% -------------------------------------------------------------------------

% Load data
clear; clc;
load jln_ferrors;
%svf = load('svfmeans.txt');
%svy = load('svymeans.txt');

%%%%
% Read in R results in different format
sf = csvread('svflatent.csv', 1); % Estimates of latent process s(t)
tf = csvread('svfparams.csv', 1); % AR(1) parameter estimators

sy = csvread('svylatent.csv', 1);
ty = csvread('svyparams.csv', 1);
%%%%


% Compute objects from predictors
h   = 12;
fb  = sparse(fbetas);
thf = [tf(1,:).*(1-tf(2,:));tf(2,:);tf(3,:).^2];
xf  = sf;
%gf  = svf(end-3+1:end,:);
[evf,phif] = jln_compute_uf(xf,thf,fb,h);

% Compute uncertainty
[T,N] = size(vyt);
ut    = zeros(T,N,h);

phimat2 = cell(N, 1);

for i = 1:N
    tic;
    yb  = sparse(ybetas(i,:));
    thy = [ty(1,i).*(1-ty(2,i));ty(2,i);ty(3,i).^2];
    xy  = sy;
    %ut(:,i,:) = jln_compute_uy(xy,thy,yb,py,evf,phif);
    [ut(:,i,:), evy2, phimat2{i}, ~, ~] = jln_compute_uy(xy,thy,yb,py,evf,phif);
    fprintf('Series %d, Elapsed Time = %0.4f \n',i,toc);
end

jlnut2 = ut;

%gy = svy(end-3+1:end,:);
save jlnresults2 dates jlnut2 evy2 phimat2
%save geweke dates gy gf names


