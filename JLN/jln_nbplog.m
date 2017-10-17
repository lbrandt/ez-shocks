function [ic1, chat,Fhat,eigval]=jln_nbplog(x,kmax,jj,DEMEAN)
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
%       x       Matrix of observed variables [T x N]
%       kmax    Maximum number of factors to consider
%       jj      Indicate which penalty term to use {1, .., 8}
%       demean  Indicate transformation of data {0, 1, 2}
%
%   Output
%       icl
%       chat
%       Fhat
%       eigval 
%
%   Dependencies {source}
%       standard {jln}
%       svd {vanilla}
%       minindc {jln}
%
% -------------------------------------------------------------------------

%function [ic1, lambda,Fhat,IC1]=nbplog(x,kmax,jj,DEMEAN);
T = size(x,1);
N = size(x,2);
NT = N*T;
NT1 = N+T;
CT = zeros(1,kmax);
ii = 1:1:kmax;

if jj ==1
    CT(1,:) = log(NT/NT1)*ii*NT1/NT;
end

if jj==2
    CT(1,:)=(NT1/NT)*log(min([N;T]))*ii;
end

GCT=min([N;T]);

if jj==3
    CT(1,:)=ii*log(GCT)/GCT;
end

if jj==4
    CT(1,:)=2*ii/T;
end

if jj==5
    CT(1,:)=log(T)*ii/T;
end
if jj==6
    CT(1,:)=2*ii*NT1/NT;
end
if jj==7
    CT(1,:)=log(NT)*ii*NT1/NT;
end
if jj==8
    CT(1,:)= 2*ii*(sqrt(N)+sqrt(T))^2/(NT);
end % new modified CP


if DEMEAN ==2
    X=jln_standard(x);
end

if DEMEAN ==1
    X=x-repmat(mean(x),T,1);
end

if DEMEAN==0
    X=x;
end

IC1 = zeros(size(CT,1),kmax+1);
Sigma = zeros(1,kmax+1);

if T < N
    [ev,eigval,ev1] = svd(X*X');
    sumeigval = cumsum(diag(eigval))/sum(diag(eigval));
    Fhat0 = sqrt(T)*ev;
    Lambda0 = X'*Fhat0/T;
else
    [ev,eigval,ev1] = svd(X'*X);
    sumeigval=cumsum(diag(eigval))/sum(diag(eigval));
    Lambda0=sqrt(N)*ev;
    Fhat0=X*Lambda0/N;
end



if jj <= 8
    for i=kmax:-1:1 % count down from kmax
        Fhat=Fhat0(:,1:i);
        %lambda=Fhat'*X;
        lambda=Lambda0(:,1:i);
        chat=Fhat*lambda';
        %disp([i sumeigval(i) sum(sum(chat.*chat))/sum(sum(X.*X))]);
        ehat=X-chat;
        Sigma(i)=mean(sum(ehat.*ehat/T));
        IC1(:,i)=log(Sigma(i)) + CT(:,i);
    end
    Sigma(kmax+1)=mean(sum(X.*X/T));
    IC1(:,kmax+1)=log(Sigma(kmax+1));
    ic1 = jln_minindc(IC1')';
    ic1 = ic1 .*(ic1 <= kmax);
end

if jj==9
    for j=1:rows(sumeigval)
        if sumeigval(j) >= .5
            ic1=j;
            break;
        end
    end
end

Fhat = [];
Fhat = Fhat0(:,1:kmax);
Lambda = Lambda0(:,1:kmax);
chat = Fhat*Lambda';
eigval = diag(eigval);

end