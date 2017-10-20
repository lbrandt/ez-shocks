% -------------------------------------------------------------------------
% Generate forecast errors
% Adapted from JLN 2015 for use with vanilla MATLAB
% -------------------------------------------------------------------------

% Load data
clear; clc;
load jlndata; % includes all 132 macro series & 157 financial series 1960:01 .. 2011:12
ind         = 132+(6:15); % column index of series to be removed, 6:15 in fin block    
data(:,ind) = []; % delete series in data matrix
names(ind)  = []; % delete series in name vector
xt          = data;

% Estimate factors
[e,fhat,lf,vf] = jln_factors(xt,20,2,2);
[e,ghat,lg,vg] = jln_factors(xt.^2,20,2,2);

% Calculate marginal R2 (unused)
outf     = 'mR2_fhat.out';
outg     = 'mR2_ghat.out';

[R2,mR2] = mrsq(fhat,lf,vf,names,vartype,outf);
[R2,mR2] = mrsq(ghat,lg,vg,names,vartype,outg);

% Build predictor set
ft       = [fhat,fhat(:,1).^2,ghat(:,1)];