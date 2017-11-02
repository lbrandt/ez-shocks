function [cmatphi] = companion(numvars, numlags, params)
% -------------------------------------------------------------------------
% Builds the parameter matrix of the companion form for a general VAR
% model with N variables and P lags. The univariate AR(P) is the special
% case with N=1.
%
% XYXY Constants!!
%
%   Input
%       numvars     Number of variables in the system N
%       numlags     Number of lags P
%       params      Parameter vector [3 x 1]
%
%   Output
%       cmatphi  	Parameter matrix of VAR companion form [NP x NP]
% 
%   Dependencies {source}
%
% -------------------------------------------------------------------------

% Read dimensions
cmatdim = numvars*numlags;

% Arrange VAR parameters
cmattop = params(:, end-cmatdim+1:end);

% Build lower section
cmatbot = [eye(cmatdim-numvars), zeros(cmatdim-numvars, numvars)];

% Combine
cmatphi = [cmattop; cmatbot];

end
