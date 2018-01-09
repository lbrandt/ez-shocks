% -------------------------------------------------------------------------
% Import German data
%
% 
%
%
% -------------------------------------------------------------------------

% German data
%clear; clc;
%addpath('..\R;..\MATLAB')


% x = h5read('gs_data.h5', '/dlndata');
% varnames = h5read('gs_data.h5'); %ATTRIBUTE?!?
% dates =  h5read('gs_data.h5'); %ATTRIBUTE?!?


gsdata   = table2array(readtable('gs_data.csv', 'ReadVariableNames', false));
varnames = table2array(readtable('gs_varnames.csv', 'ReadVariableNames', false));
dates    = datetime(table2array(readtable('gs_dates.csv', 'ReadVariableNames', false)));

[T, N]   = size(gsdata);

fprintf('Monthly series from %s to %s \n', datestr(dates(1)), datestr(dates(end)));
fprintf('Sample: T = %d, N = %d \n', T, N);

% Format names such that they are permissible MATLAB varnames
varnames = strrep(varnames, ' ', '_');
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

save gs_data varnames dates gsdata




%%%% OLD QUARTERLY
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
