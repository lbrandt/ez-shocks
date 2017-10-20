function x = standardise(y)
% -------------------------------------------------------------------------
% Standardises a matrix to have zero means and unit variances in columns.
%
%   Input
%       y       Matrix with different means and variances [M x N]
%
%   Output
%       x       Matrix of (0,1) columns
% 
%   Dependencies {source}
%
% -------------------------------------------------------------------------


[M, ~] = size(y);

x = (y - repmat(mean(y), M, 1))./repmat(std(y), M, 1);
    
end
