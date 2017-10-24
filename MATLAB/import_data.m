% -------------------------------------------------------------------------
% Import JLN data
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

[T, N]   = size(x);

ta = datetime(1960, 03, 15);
te = datetime(2011, 12, 15);
dates = (ta:calmonths(1):te)';

% Format names such that they are permissible MATLAB varnames
varnames = strrep(names, ' ', '_');
varnames = strcat('var_', varnames);
varnames = strrep(varnames, '&', '_');
varnames = strrep(varnames, '/', '_');
varnames = strrep(varnames, '\', '_');
varnames = strrep(varnames, ':', '');
varnames = strrep(varnames, '.', '');
varnames = strrep(varnames, ',', '');
varnames = strrep(varnames, '+', '');
varnames = strrep(varnames, '-', '');
varnames = strrep(varnames, '<', '');
varnames = strrep(varnames, '>', '');
varnames = strrep(varnames, '(', '');
varnames = strrep(varnames, ')', '');
varnames = strrep(varnames, '’', '');


fprintf('Monthly series from %s to %s \n', datestr(ta), datestr(te));
fprintf('Sample: T = %d, N = %d \n', T, N);

