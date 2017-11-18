% -------------------------------------------------------------------------
% Import German data
%
% 
%
%
% -------------------------------------------------------------------------

% German data
%clear; clc;

x = csvread('de_data2.csv', 1, 1); % logdiffs

names = []



[T, N]   = size(x);

ta = datetime(1991, 06, 30);
te = datetime(2017, 03, 31);
dates = dateshift((ta:calquarters(1):te)', 'end', 'month');


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




fprintf('Quarterly series from %s to %s \n', datestr(ta), datestr(te));
fprintf('Sample: T = %d, N = %d \n', T, N);
