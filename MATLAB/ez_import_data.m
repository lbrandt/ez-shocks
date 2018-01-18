% -------------------------------------------------------------------------
% Import Eurozone data
%
%
%
%
% -------------------------------------------------------------------------

% EZ data QUARTERLY
%clear; clc;
%addpath('..\R;..\MATLAB')


%h5disp('ez_data90_q.h5')
varnames = h5read('ez_data90_q.h5', '/varnames');
dates    = datetime(h5read('ez_data90_q.h5', '/dates'), 'InputFormat', 'yyyyQQQ', 'Format', 'yyyyQQQ');
x        = h5read('ez_data90_q.h5', '/dlndata')';

[T, N]   = size(x);

fprintf('Quarterly series from %s to %s \n', datestr(dates(1)), datestr(dates(end)));
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

save ez_data90_q varnames dates x
