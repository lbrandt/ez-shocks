function [RMSE, Bias, Variance] = rmse(x, xhat)
% -------------------------------------------------------------------------
% Computes the Root Mean Squared Error between a matrix of true, observed
% data and a matrix of estimates. Mean is computed over rows, ergo data
% should be input [N x M], i.e. [observations x variables]. Additionally
% computes bias and error variance (n-1).
%
%
%   Input
%       x           Matrix of observed variables [N x M]
%       xhat        Matrix of estimates of x [N x M]
%
%   Output
%       RMSE        Root Mean Squared Error
%       Bias        Bias(x, xhat)
%       Variance    Residual Variance
% 
%   Dependencies {source}
%
% -------------------------------------------------------------------------

if size(x) ~= size(xhat)
    warning('Dimensions of x and xhat must align.');
elseif size(x) == size(xhat)
    
    errors   = x - xhat;
    
    sqerrors = errors.^2;

    RMSE     = sqrt(mean(sqerrors));
    
    Bias     = mean(errors);
    Variance = var(errors);
end

end
