function x=jln_standard(y)
% -------------------------------------------------------------------------
% wyd 
%
%   Input
%       
%
%   Output
%
% 
%   Dependencies {source}
%
%
% -------------------------------------------------------------------------
T=size(y,1);
N=size(y,2);
my=repmat(mean(y),T,1);
sy=repmat(std(y),T,1);
x=(y-my)./sy;

%x=(y-kron(mean(y),ones(rows(y),1)))./kron(std(y),ones(rows(y),1));
end