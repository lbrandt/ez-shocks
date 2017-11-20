% -------------------------------------------------------------------------
% Import German data
%
% 
%
%
% -------------------------------------------------------------------------

% German data
%clear; clc;

[~, names, ~] = xlsread('de_varnames.csv');
[~, dates, ~] = xlsread('de_dates.csv');

x = csvread('de_data2.csv'); % logdiffs

[T, N]   = size(x);
ta       = dates{2}; % diffed series lack one observation compared to raw dataset
te       = dates{end};


% Dates in MATLAB datetime format
%ta = datetime(1991, 06, 30);
%te = datetime(2017, 03, 31);
%dates = dateshift((ta:calquarters(1):te)', 'end', 'month');


% Format names such that they are permissible MATLAB varnames
varnames = strrep(names, ' ', '_');
%varnames = strcat('var_', varnames);
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


%fprintf('Quarterly series from %s to %s \n', datestr(ta), datestr(te));
fprintf('Quarterly series from %s to %s \n', ta, te);
fprintf('Sample: T = %d, N = %d \n', T, N);
