%%%%
% Static factor model

% Simulate static factor model with R latent factors and M observed variables
N = 100;

% Observed variables X
M = 8;

% Common factors F
R = 4;

%%%% function random() requires Statistics and Machine Learning Toolbox!
%fdist = 'beta';
%fparn = 2;
%fpara = 1;
%fparb = 4;

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
Rhat = 4; % Number of principal components to be extracted

[Fhat, Lhat, Ehat, eigval] = factors(X, Rhat, demean);

% Analyse results
% Matrix of correlations between true factors and estimated factors
fcorrmat = zeros(R, Rhat);

for i = 1:R
    for j = 1:Rhat
        corrcoeffs     = corrcoef(F(:, i), Fhat(:, j));
        fcorrmat(i, j) = corrcoeffs(2, 1);
    end
end

figure
subplotindex = 1;
for i = 1:R
    for j = 1:Rhat
        subplot(R, Rhat, subplotindex)
        
        plot(standardise(F(:, i)))
        hold on
        plot(Fhat(:, j))
        title(['F(', num2str(i), '), Fhat(', num2str(j),'), rho = ', num2str(fcorrmat(i, j))])
        
        subplotindex = subplotindex + 1;
    end
end
        
plot(standardise(F(:, 1)))
%hold on
%plot(Fhat(:, 1))
hold on
plot(Fhat(:, 2))


% Structural errors ~N(0, s^2)
S = [2, 0; 0, 1]';
E = normrnd(0, 1, [N, 2]) *S;

% Exogeneous relationship Beta
B = [0.4, 0; 0, -0.6];

% Endogeneous relationship Gamma
G = [1, 0; 0, 1];

% Model GY' = BX' + E
Y = (B*X' + E')';

scatter(X(:, 1), Y(:, 1))

bhat = inv(X'*X)* X'*Y;




