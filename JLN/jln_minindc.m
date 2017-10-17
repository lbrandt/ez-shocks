function pos=jln_minindc(x)
% -------------------------------------------------------------------------
% Indicates the position of the smallest value in columns of a matrix.
% GAUSS compatibility function.
%
%   %%%% USE MININD(x,1) INSTEAD
%
%   Input
%       x       arbitrary matrix [M x N]
%
%   Output
%       pos     vector with position indices [N x 1]
% 
%   Dependencies {source}
%       seqa {ls}
%
% -------------------------------------------------------------------------

ncols=size(x,2);
nrows=size(x,1);
pos=zeros(ncols,1);
seq=seqa(1,1,nrows);

for i=1:ncols
    dum = min(x(:,i));
    dum1 = seq .* ( (x(:,i)-dum) ==0);
    pos(i) = sum(dum1);
end

end
