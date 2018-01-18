% -------------------------------------------------------------------------
% Estimate impulse response functions
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB')

% Load data
%load de_shocks
%load de_uncertainty

% Convert uncertainty to High-U-dummies via quantiles
%if U > F_U(0.6)
%    Udum = 1
%else
%    Udum = 0

% IRF via Jorda's local projection method

%%%% LIKE THIS?!?
% Reg BIP ~ shocks

% Reg BIP ~ shocks, uncertainty, shocks*uncertainty

% Reg BIP ~ shocks, uncertainty, shocks*Udum
