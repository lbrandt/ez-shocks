function z = lead(x , k)
% -------------------------------------------------------------------------
% Create a matrix of k leads of a vector or matrix
%   Input:  x = matrix or vector, (nobs x nvar)
%           k = number of leads, (default = 1)
%   Output: z = matrix (or vector) of leads (nobs x nvar*k)
% -------------------------------------------------------------------------

if nargin == 1
    k = 1;
end

[nobs, nvar] = size(x);

z = zeros(nobs, nvar*k);

for j = 1:k
    z(1:nobs-j, nvar*(j-1)+1:j*nvar) = x(1+j:nobs, :);
end

end
