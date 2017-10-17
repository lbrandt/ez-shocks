function [ic1, chat,Fhat,eigval] = bninfocrit(x, kmax, gnum, demean)
% -------------------------------------------------------------------------
% Estimates number of factors according to some information criterion as in
% Bai & Ng (2002) under asymptotic N,T.
%
% Information criteria
%
%
% If demean == 0, do not transform, take raw data.
% If demean == 1, demean data to (0, sigma).
% If demean == 2, standardise data to (0, 1).
%
% (Bai & Ng, 2002)
%
%   Input
%       x           Matrix of observed variables [T x N]
%       kmax        Maximum number of factors to consider
%       gnum        Indicate penalty term to be used {1, .., 8}
%       demean      Indicate transformation of data {0, 1, 2}
%
%   Output
%       ic1         
%       chat
%       Fhat
%       eigval 
%
%   Dependencies {source}
%       standardize {oxford}
%       svd {vanilla}
%       minindc {jln}
%
% -------------------------------------------------------------------------

[T, N] = size(x);

NTprd = N*T;
NTsum = N+T;

CNT = min([sqrt(N), sqrt(T)]);



% Transform data according to demean
if demean == 2
    xt = standardize(x);
end

if demean == 1
    xt = x-repmat(mean(x),T,1);
end

if demean == 0
    xt = x;
end




kgNT = zeros(1,kmax);

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
end





   

end