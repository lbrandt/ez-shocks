% -------------------------------------------------------------------------
% Import data
%
% 
%
%
% -------------------------------------------------------------------------

% JLN 2015
%clear; clc;
load jlndata; 
ind         = 132+(6:15); % "duplicate" series to remove    
data(:,ind) = []; 
names(ind)  = [];
x           = data;

fprintf('Sample: T = %d, N = %d \n', size(x));
