% Test factors

% Load data from generate_ferrors
clear; clc;
load jlndata; 
ind         = 132+(6:15); % "duplicate" series to remove    
data(:,ind) = []; 
names(ind)  = [];
xt          = data;

% Estimate factors
[e,fhat,lf,vf] = jln_factors(xt,20,2,2);
[e,ghat,lg,vg] = jln_factors(xt.^2,20,2,2);


jln_fhat = fhat;


[ic1,chat,fhat,eigval]  = jln_nbplog(xt,20,2,2);

icstar    = ic1;

[ehat,Fhat,lamhat,ve2]  = jln_pc(jln_standard(xt),icstar);




%%%%%%%
% from testnfac.m
clear; clc;

r=4;
N=100;
T=50;
randn('state',999);

e=randn(T,N);
f=randn(T,r);
lambda=randn(N,r);
x=f*lambda'+e;

rmax=10;
DEMEAN=2;
fprintf('Demean %d \n',DEMEAN);
disp('Determining number of factors');
fprintf('T= %d N= %d r = %d \n', size(x), r)

for i=1:8
  disp(jln_nbplog(x,rmax,i,DEMEAN));
end;  


[T, N] = size(x);


x = randn(4,3)
xx = x'*x;


[U, S, V] = svd(xx);
[E, D] = eigsort(xx);

disp([U, E])
disp([S, D])


