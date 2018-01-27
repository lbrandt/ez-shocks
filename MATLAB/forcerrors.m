% Forecast Evaluation
% 

regstart = 2;

h = 1;
window = 100;

neq = 1;
festart = regstart;
festop = taindex + window;

femax = regend - festop + 1 - h;

ff = zeros(femax, 1);
fe = zeros(femax, 1);
fu = zeros(femax, 1);
fl = zeros(femax, 1);


for i = 1:femax
    
    sample = festart:festop;
    beta = ols(Y(sample), X(sample), 0);
    
    ff(i) = X(sample + h)* beta;
    
    fe(i) = Y(sample + h) - forc(i);
end

ssfe = sum(fe.^2);
mse = ssfe/femax;
rmse = mse^0.5;

% Constructing Point and Interval Forecasts
coverage = 0.6;


for i = 1:femax
    
    fu(i) = ff(i) + norm_inv(0.5 + 0.5*coverage)*rmse;
    fl(i) = ff(i) - norm_inv(0.5 + 0.5*coverage)*rmse;
    
end
