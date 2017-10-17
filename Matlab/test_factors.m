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
R2_static = sum(eigval(1:icstar))/sum(eigval);

[ehat,Fhat,lamhat,ve2]  = jln_pc(jln_standard(xt),icstar);




%%%%%%%
[T, N] = size(xt);




