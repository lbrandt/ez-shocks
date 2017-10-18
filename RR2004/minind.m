function pos = minind(x, dim)
% -------------------------------------------------------------------------
% Indicates the position of the smallest value in columns or rows
% of an arbitrary matrix x.
% If dim == 1, search in columns, return row index.
% If dim == 2, search in rows, return column index.
% If dim is not supplied, default to dim = 1.
%
%   Input
%       x       Matrix [M x N]
%       dim     Search dimension to be passed to min() {1, 2}
%
%   Output
%       pos     Column vector with position indices
% 
%   Dependencies {source}
%       {}
%
% -------------------------------------------------------------------------

if nargin < 2 || isempty(dim) % default to dim = 1 if not supplied.
    dim = 1;
end

if dim == 1
    [pos, ~] = find(x == min(x,[],dim)); % search for min in columns and return vector of row indices.
elseif dim == 2
    [~, pos] = find(x == min(x,[],dim)); % search for min in rows and return vector of column indices.
else
    warning('"dim" can only take up values 1 or 2')
end

end
