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
dates    = dates(2:end); % diffed series lack one observation compared to raw dataset

% Dates in MATLAB datetime format
ta = datetime(dates{1}, 'InputFormat', 'QQ yyyy', 'Format', 'QQQ yyyy');
te = datetime(dates{end}, 'InputFormat', 'QQ yyyy', 'Format', 'QQQ yyyy');
dates = (ta:calquarters(1):te)';


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

save de_data dates varnames x
