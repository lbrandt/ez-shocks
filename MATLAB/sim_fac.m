%%%%
% Static factor model

%clear;clc;
rng(58);

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

% RMSE(F_1, Fhat_1) = 0.3277
% RMSE(F_2, Fhat_2) = 0.7149


%%%%
% Duplicate series
% Factor loadings
Lhat

% Here, first observable X_1 is highly correlated with second factor. Thus,
% saturate X matrix with X_1 and check stability of factor estimators. When
% X_1 is simply duplicated (rank(X) = M) and when Normal disturbance is
% added (rank(X2) = M + duplicate), identification eventually switches. 
% The first factor is then represented by the second extracted component 
% and the second factor is represented by the first extracted component. 
% This is likely due to its variance eventually dominating the eigenvalues 
% of the covariance matrix. This can lead to problems when too few
% principal components are extracted. Then the true structural factor space
% will not be recovered but instead only that part which is dominating the 
% sample variance. It might therefore be sensible to choose a more lenient
% information criterion for estimating R.

X2 = X;

duplicate = 20;

fcorrmat2 = zeros(R, Rhat, duplicate);
frmse2    = zeros(R, Rhat, duplicate);

for s = 1:duplicate
    
    X2 = [X2, X2(:, 1) + randn(N, 1)];
    
    [Fhat2, Lhat2, Ehat2, eigval2] = factors(X2, Rhat, demean);
    
    for i = 1:R
        for j = 1:Rhat
            corrcoeffs2 = corrcoef(F(:, i), Fhat2(:, j));
            fcorrmat2(i, j, s) = corrcoeffs2(2, 1); % Extract correlation
            
            frmse2(i, j, s)    = rmse(standardise(F(:, i)), standardise(Fhat2(:, j)));
        end
    end
    
end
    
figure
subplotindex = 1;
for i = 1:R
    for j = 1:Rhat
        subplot(R, Rhat, subplotindex)
        
        plot(squeeze(abs(fcorrmat2(i, j, :))))
        %plot(squeeze(frmse2(i, j, :)))
        title(['F(', num2str(i), '), Fhat(', num2str(j), ')'])
        
        subplotindex = subplotindex + 1;
    end
end




%%%%
% Add disturbed copies of raw factors to data
% Does this approximate including survey indicators?

% When a factor is directly included in the data set (plus some Normal
% disturbance, its representation by the first extracted component improves
% almost monotonously. The identification of the second factor fluctuates
% but does not suffer systematically. Therefore one might assume that it is
% not harmful to identification to include survey indicators in the data
% set even if they might already represent the underlying factors.
% A linear combination of the two factors does not significantly worsen
% convergence either. However, a product of the two factors apparently
% dilutes the data in such a way that identification deteriorates a lot. 
X3 = X;

addfactor = 20;

fcorrmat3 = zeros(R, Rhat, addfactor);
frmse3    = zeros(R, Rhat, addfactor);

for s = 1:addfactor
    
    X3 = [X3, F(:, 1) + F(:, 2) + randn(N, 1)];
    
    [Fhat3, Lhat3, Ehat3, eigval3] = factors(X3, Rhat, demean);
    
    for i = 1:R
        for j = 1:Rhat
            corrcoeffs3 = corrcoef(F(:, i), Fhat3(:, j));
            fcorrmat3(i, j, s) = corrcoeffs3(2, 1); % Extract correlation
            
            frmse3(i, j, s)    = rmse(standardise(F(:, i)), standardise(Fhat3(:, j)));
        end
    end
    
end
    
figure
subplotindex = 1;
for i = 1:R
    for j = 1:Rhat
        subplot(R, Rhat, subplotindex)
        
        plot(squeeze(abs(fcorrmat3(i, j, :))))
        %plot(squeeze(frmse3(i, j, :)))
        title(['F(', num2str(i), '), Fhat(', num2str(j), ')'])
        
        subplotindex = subplotindex + 1;
    end
end
