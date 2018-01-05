%%%%
% Linear regression model

% Simulate bivariate regression model Y = X*B
N = 100;

% Exogeneous variables
M = 2;

exodist = 'beta';
exoparn = 2;

exopars = unidrnd(6, [exoparn, M]);

% Exogeneous variables
X = zeros(N, M);

for j = 1:M
    for i = 1:N
        X(i, j) = random(exodist, exopars(1, j), exopars(2, j));
        
    end
end

X = X*100;


xbar = mean(X);
xcov = cov(X(:, 1), X(:, 2));

scatter(X(:, 1), X(:, 2))

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




