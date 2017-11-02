function [evf,phif] = compute_uf(xf,tf,fb,h)
% -------------------------------------------------------------------------
% Compute expected volatility of predictors up to horizon h,and construct
% coefficient matrix phiF
%
%   Input
%       xf          Latent?
%       tf          Parameter vector of logVar process [3 x 1]
%       fb          ?
%       h           Max forecast horizon of expected volatility
%
%   Output
%       U        	?
% 
%   Dependencies {source}
%
% -------------------------------------------------------------------------

% Initialize parameters
R   = size(xf, 2);
pf  = size(fb, 2)-1;

% Create phif matrix
phif_top = [];
for j = 2:pf+1
    phif_top = [phif_top,sparse(1:R,1:R,fb(:,j))]; 
end

if pf > 1
    phif_bot = [sparse(1:R*(pf-1),1:R*(pf-1),1),sparse(R*(pf-1),R,0)];
    phif     = [phif_top;phif_bot];
else
    phif = phif_top;
end

% Compute evf
evf = cell(h,1);

for j = 1:h
    for i = 1:R
        alpha       = tf(1,i);
        beta        = tf(2,i);
        tau2        = tf(3,i);
        x           = xf(:,i);
        evf{j}(:,i) = expectvar(alpha,beta,tau2,x,j); %Et[(v^f_t)^2]
    end
end

end
