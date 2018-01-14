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


%h5disp('gs_data.h5')
varnames = h5read('de_data.h5', '/varnames');
dates    = datetime(h5read('de_data.h5', '/dates'));
x        = h5read('de_data.h5', '/dlndata')';

[T, N]   = size(x);

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

% Shorten dates of diffed series
dates = dates(2:end);

save de_data varnames dates x




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
