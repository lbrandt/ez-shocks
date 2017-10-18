
tic

x = randn(10,10);

toc



tic

%pos = jln_minindc(x);
pos = minind(x,1);

toc

[row, col] = find(x == min(x,[],1))

% minindc

tic

ncols = size(x,2);
nrows = size(x,1);
pos = zeros(ncols,1);
seq = seqa(1,1,nrows);


for i=1:ncols
    dum = min(x(:,i));
    dum1 = seq .* ( (x(:,i)-dum) ==0);
    pos(i) = sum(dum1);
end

toc