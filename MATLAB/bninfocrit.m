function [numfac, IC, Fhat, Lhat, eigval] = bninfocrit(x, kmax, gnum, demean)
% -------------------------------------------------------------------------
% Estimates number of factors according to some information criterion as in
% Bai & Ng (2002) under asymptotic N: number of series, T: length of series.
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
%       kmax        Maximum number of factors to consider
%       gnum        Indicate penalty term to be used {1, .., 10}
%       demean      Indicate transformation of data {0, 1, 2}
%
%   Output
%       numfac      Number of factors according to chosen IC {1, .., kmax}         
%       Fhat        Factors corresponding to first kmax Eigenvalues
%       Lhat        Factor loadings corresponding to kmax factors
%       eigval      Vector of unscaled Eigenvalues [min(T, N)]
%
%   Dependencies {source}
%       svd {vanilla}
%       minind {lb}
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
    sumeigval = cumsum(diag(S))/sum(diag(S)); % cumsum of scaled eigenvalues 
    Fhat0 = sqrt(T)*EV;
    Lhat0 = xtr'*Fhat0 /T;
    
else
    % Determine F and Lambda via Lbar = sqrt(N) * EV(k)
    [EV, S, ~] = svd(xtr'*xtr);
    sumeigval = cumsum(diag(S))/sum(diag(S)); % cumsum of scaled eigenvalues 
    Lhat0 = sqrt(N)*EV;
    Fhat0 = xtr*Lhat0 /N;

end


%%%%
% Information Criteria
NTprd = N*T;
NTsum = N+T;

CNT = min([sqrt(N), sqrt(T)]);


% Choice of unscaled penalty term g(N, T) at each value k:
k = 1:1:kmax;

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





VF = zeros(1, kmax+1);
IC = zeros(1, kmax+1);

if gnum <= 10
    
    for k = kmax:-1:1 % count down from kmax
        % Compute chosen Information Criterion for each k
        Fhat = Fhat0(:, 1:k); % Select first k factors
        Lhat = Lhat0(:, 1:k); % Select corresponding factor loadings
        
        Chat = Fhat*Lhat'; % Compute common component
        ehat = xtr - Chat; % Compute idiosyncratic component
        
        VF(k) = mean(sum(ehat.*ehat /T)); % Sigmahat = V(k, Fk)
        IC(:, k) = log(VF(k)) + kgNT(k);
    end
    
    % Non-penalty case at kmax+1
    VF(kmax+1) = mean(sum(xtr.*xtr /T));
    IC(kmax+1) = log(VF(kmax+1));
    
    % Select k with optimal IC as consistent estimator of r
    numfac = minind(IC,2);
    numfac = numfac .*(numfac <= kmax);
    
else  % Number of factors by explained variance
    for j = 1:rows(sumeigval)
        if sumeigval(j) >= 0.5
            numfac=j;
            break;
        end
    end
end

% Output kmax-restricted factor space and loadings
IC = min(IC);
Fhat = Fhat0(:, 1:kmax);
Lhat = Lhat0(:, 1:kmax);
eigval = diag(S); % unscaled Eigenvalues

   

end
