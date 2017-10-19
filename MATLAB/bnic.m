function [BNIC, Fhat, Lhat, eigval] = bnic(x, k, gnum, demean)
% -------------------------------------------------------------------------
% Compute value of some chosen information criterion as in Bai & Ng (2002).
% Adapted from nbplog.m by Serena Ng.
%
% Choose Information Criterion by penalty term g(N, T) via gnum. Right now,
% only ICs suggested by BN2002 are properly implemented (JLN2015 use ICp2).
% 1. ICp1
% 2. ICp2
% 3. ICp3
%
% Transforms data before Eigendecomposition:
% If demean == 0, do not transform, take raw data.
% If demean == 1, only demean data.
% If demean == 2, standardise data to (0, 1).
%
%
%   Input
%       x           Matrix of observed variables [T x N]
%       k        	Number of factors
%       gnum        Indicate penalty term to be used {1, .., 10}
%       demean      Indicate transformation of data {0, 1, 2}
%
%   Output
%       BNIC        Value of chosen IC at value k         
%       Fhat        Factors corresponding to first k Eigenvalues
%       Lhat        Factor loadings corresponding to k factors
%       eigval      Vector of unscaled Eigenvalues [min(T, N)]
%
%   Dependencies {source}
%       svd {vanilla}
%
% -------------------------------------------------------------------------

[T, N] = size(x);

% Transform data according to demean
if demean == 2
    xtr = (x - repmat(mean(x),T,1))./repmat(std(x),T,1);
    
elseif demean == 1
    xtr = x - repmat(mean(x),T,1);
    
elseif demean == 0
    xtr = x;
end

%%%%
% Eigenvalue Decomposition of Covariance Matrix to estimate unrestricted factor space
if T < N
    % Determine F and Lambda via Ftilde = sqrt(T) * EV(k)
    [EV, S, ~] = svd(xtr*xtr');
    eigval = diag(S);
    Fhat0 = sqrt(T)*EV;
    Lhat0 = xtr'*Fhat0 /T;
    
else
    % Determine F and Lambda via Lbar = sqrt(N) * EV(k)
    [EV, S, ~] = svd(xtr'*xtr);
    eigval = diag(S);
    Lhat0 = sqrt(N)*EV;
    Fhat0 = xtr*Lhat0 /N;

end


%%%%
% Information Criteria
NTprd = N*T;
NTsum = N+T;

CNT = min([sqrt(N), sqrt(T)]);


% Choice of unscaled penalty term g(N, T) at k:

if gnum == 1
    kgNT = k* NTsum/NTprd*log(NTprd/NTsum); % p1
    
elseif gnum == 2
    kgNT = k* (NTsum/NTprd)*log(CNT^2); % p2

elseif gnum == 3
    kgNT = k* log(CNT^2)/(CNT^2); % p3

elseif gnum == 4
    kgNT = k* 2/T; % AIC1
    
elseif gnum == 5
    kgNT = k* log(T)/T; % BIC1
    
elseif gnum == 6
    kgNT = k* 2/N; % AIC2
    
elseif gnum == 7
    kgNT = k* log(N)/N; % BIC2

elseif gnum == 8
    kgNT = k* 2*(NTsum - k)/NTprd; % AIC3

elseif gnum == 8
    kgNT = k* (NTsum - k)*log(NTprd)/NTprd; % BIC3

elseif gnum == 10
    kgNT = k* 2*(sqrt(N)+sqrt(T))^2/(NTprd); % JLN: new modified CP

else
    text = ['To compute proper Information Criteria, "gnum" can only take up integers {1, .., 10}.' ...
            'In other cases, number of factors is chosen via explained variance >= 0.5.'];
    warning(text)
end



if gnum <= 10
    
    % Compute chosen Information Criterion for k factors
    Fhat = Fhat0(:, 1:k); % Select k factors
    Lhat = Lhat0(:, 1:k); % Select corresponding factor loadings
    
    Chat = Fhat*Lhat'; % Compute common component
    ehat = xtr - Chat; % Compute idiosyncratic component
    
    VF = mean(sum(ehat.*ehat /T)); % Sigmahat = V(k, Fk)
    BNIC = log(VF) + kgNT; % Value of chosen IC
    
else
    BNIC = NaN;
end


end
