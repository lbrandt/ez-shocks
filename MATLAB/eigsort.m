function [E, D] = eigsort(A)
% -------------------------------------------------------------------------
% Calculates Eigendecomposition of a square matrix via MATLAB's eig
% but returns eigenvalues and left eigenvectors in descending order.
% (adapted from eigsort.m by D. R. Bohnenstiel (NCSU))
%
%   Input
%       A           Square Matrix [M x M]
%
%   Output
%       E           Left Eigenvectors of A
%       D           Eigenvalues of A sorted in descending order
% 
%   Dependencies {source}
%       {}
%
% -------------------------------------------------------------------------

if nargin > 1 % eig(A,B) not implemented
    warning('Input only one matrix, eig(A,B) is not implemented.')
end

[E, D] = eig(A);
[D, index] = sort(diag(D), 'descend');

D = diag(D);
E = E(:, index);

end
