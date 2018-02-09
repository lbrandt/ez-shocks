% -------------------------------------------------------------------------
% Import Eurozone data
%
%
%
%
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB')

% EZ data monthly

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

save ez_data names dates x




% EZ data QUARTERLY
%h5disp('ez_data90_q.h5')
names = h5read('ez_data90_q.h5', '/varnames');
dates    = datetime(h5read('ez_data90_q.h5', '/dates'), 'InputFormat', 'yyyyQQQ', 'Format', 'yyyyQQQ');
x        = h5read('ez_data90_q.h5', '/dlndata')';

[T, N]   = size(x);

fprintf('Quarterly series from %s to %s \n', datestr(dates(1)), datestr(dates(end)));
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

% Shorten dates of diffed series
dates = dates(2:end);

save ez_data90_q names dates x
