%%%%
% Static factor model

%clear;clc;

% Simulate static factor model with R latent factors and M observed variables
N = 100;

% Observed variables X
M = 5;

% Common factors F
R = 2;

%%%% function random() requires Statistics and Machine Learning Toolbox!
%fdist = 'beta';
%fparn = 2;
%fpara = 1;
%fparb = 4;

%%%% LeSage toolbox contains dist_rnd functions, use beta_rnd
%beta_rnd


% Vanilla MATLAB supplies Normal random numbers via randn
% Parametrisation via discrete Uniform
fmu = randi([-5, 5], 1, R);
fsd = randi(6, 1, R);

% Generate factors
F = zeros(N, R);

for j = 1:R
    for i = 1:N
        %F(i, j) = random(fdist, a, b);
        F(i, j) = fmu(j) + fsd(j)* randn;
        
    end
end

% Factor loadings
L = zeros(M, R);

for j = 1:R
    for i = 1:M
        % L(i, j) = random();
        L(i, j) = randn;
    end
end

% Idiosyncratic errors
E = randn(N, M); % structural errors, no cross correlation

%S = gallery('lehmer', M);
%U = E*S


% Generate observable variables X
X = F*L' + E;

% Descriptive statistics
fbar = mean(F);
fcov = cov(F);
frho = corrcoef(F);

fprintf('Descriptive Statistics \n');
fprintf('---------------------- \n');
fprintf('Theoretical mean fmu = \n');
disp(fmu);
fprintf('Empirical mean fbar = \n');
disp(fbar);
fprintf('Factors are uncorrelated: \n');
fprintf('frho = \n');
disp(frho);

%scatter(F(:, 1), F(:, 2))

xbar = mean(X);
xcov = cov(X);
xrho = corrcoef(X);

fprintf('Observables are not: \n');
fprintf('min(xrho) = %d, max(xrho) = %d \n', min(min(xrho)), max(max(xrho - eye(M)))); % Call function twice to obtain minmax across all rows and columns, subtract eye to get rid of ones on diagonal

%scatter(X(:, 1), X(:, 2))

%%%% corrplot requires Econometrics Toolbox!
%corrplot(X)



%%%%
% PCA
demean = 2; % Standardise variables
Rhat = 2; % Number of principal components to be extracted

[Fhat, Lhat, Ehat, eigval] = factors(X, Rhat, demean);

% Analyse results
% Matrix of correlations between true factors and estimated factors
fcorrmat = zeros(R, Rhat);
frmse    = zeros(R, Rhat);

for i = 1:R
    for j = 1:Rhat
        corrcoeffs     = corrcoef(F(:, i), Fhat(:, j));
        fcorrmat(i, j) = corrcoeffs(2, 1); % Extract correlation
        
        frmse(i, j)    = rmse(standardise(F(:, i)), standardise(Fhat(:, j)));
    end
end

figure
subplotindex = 1;
for i = 1:R
    for j = 1:Rhat
        subplot(R, Rhat, subplotindex)
        
        plot(standardise(F(:, i)))
        hold on
        plot(standardise(Fhat(:, j)))
        title(['F(', num2str(i), '), Fhat(', num2str(j),'); rho = ', num2str(fcorrmat(i, j)), ', RMSE = ', num2str(frmse(i, j))])
        
        subplotindex = subplotindex + 1;
    end
end

% Did PCA recover some rotation of the underlying factor space?
[frmse1, fbias, fvar] = rmse(F(:, 1), Fhat(:, 1));


