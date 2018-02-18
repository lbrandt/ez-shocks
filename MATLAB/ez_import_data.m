% -------------------------------------------------------------------------
% Imports Eurozone data that was previously cleaned and transformed in R
% into MATLAB and saves balanced panel of macroeconomic series to .mat
% workspace (version 7.3).
%
% Lennart Brandt, Feb 2018
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB;')


%h5disp('ez_data.h5')
names = h5read('ez_data.h5', '/names');
dates = datetime(h5read('ez_data.h5', '/dates'), 'InputFormat', 'yyyy-MM-dd', 'Format', 'yyyy-MM-dd');
x     = h5read('ez_data.h5', '/tdata')';

[T, N]   = size(x);

fprintf('Monthly series from %s to %s \n', datestr(dates(1)), datestr(dates(end)));
fprintf('Sample: T = %d, N = %d \n', T, N);

% Format names such that they are permissible MATLAB varnames
names = strrep(names, ' ', '_');
names = strrep(names, '&', '_');
names = strrep(names, '/', '_');
names = strrep(names, '\', '_');
names = strrep(names, ':', '');
names = strrep(names, '.', '');
names = strrep(names, ',', '');
names = strrep(names, '+', '');
names = strrep(names, '-', '');
names = strrep(names, '<', '');
names = strrep(names, '>', '');
names = strrep(names, '(', '');
names = strrep(names, ')', '');
names = strrep(names, '’', '');

% Shorten diffed series
dates = dates(2:end);
x = x(2:end, :);

save ez_data -v7.3 names dates x
