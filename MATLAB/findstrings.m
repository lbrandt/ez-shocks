function [index] = findstrings(where, what)
% -------------------------------------------------------------------------
% Searches for for multiple cells in cell of strings and returns their
% position. Works like strcmp(where, what) but allowing [what] to contain
% multiple cell of strings.
%
%   Input:  where       x
%           what        x
%   Output: indices     x
% -------------------------------------------------------------------------

indexFind = zeros(length(where), length(what));

for i = 1:length(what)
    indexFind(:, i) = strcmp(where, what(i));
    if sum(indexFind(:, i)) == 0
        warning(['Could not find ', what{i}])
    end
end

[index, ~] = find(indexFind);
end
