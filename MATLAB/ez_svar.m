% EZ SVAR


%clear; clc;
%addpath('..\R;..\MATLAB;..\..\..\Data')

% ----------------
% Load data


%h5disp('ez_vardata.h5')
vdates = datetime(h5read('ez_vardata.h5', '/dates'));
vnames = h5read('ez_vardata.h5', '/varnames');
vdata  = h5read('ez_vardata.h5', '/data')';

% Choose variables for small monetary VAR
prices = {'EKCPHARMF', 'EMCPCOR5F', 'EKCPCOREF'};
activity = {'EKIPMAN.G', 'EKIPTOT.G'};
interest = {'EMPRATE.', 'EMIBOR3.', 'EMIBOR1Y'};
money = {'EMM1....B', 'EMM2....B', 'EMM3....B'};
exchange = {'USEURSP'};

figure
subplot(2, 2, 1)
plot(vdates, vdata(:, findstrings(vnames, prices)))
subplot(2, 2, 2)
plot(vdates, vdata(:, findstrings(vnames, activity)))
subplot(2, 2, 3)
plot(vdates, vdata(:, findstrings(vnames, interest)))
subplot(2, 2, 4)
plot(vdates, vdata(:, findstrings(vnames, money)))


selectVariables = {'EKIPMAN.G', 'EKCPHARMF', 'EMIBOR3.', 'EMM1....B'}; % Ordered like CEE1999
y = vdata(:, findstrings(vnames, selectVariables));
dy = diff(y);
[T, N] = size(dy);

nlag = 12;

var1 = vare(dy, nlag);
%plt_var(var1,  char(selectVariables))

% Reduced form residuals in matrix
ehat = zeros(T-nlag, N);
for i = 1:N
    ehat(:, i) = var1(i).resid;
end

% Covariance matrix of reduced form sigma_eta
sigma = ehat'*ehat/(T-nlag-(N*nlag+1));

% Identify contemporaneous parameter matrix via Cholesky
chols = chol(sigma, 'lower');

% Normalise diagonal of Hjj to unity
%H = chols*diag(diag(chols).^2)^(-1/2);
H = chols./diag(chols)';

% Identify structural shocks from reduced form residuals
epsilon = ehat*inv(H)';
summarize(epsilon);
plot(epsilon)


eps2 = ehat*inv(chols)';
summarize(eps2);
plot(vdates((2+nlag):end), eps2)



% Euro OIS
%h5disp('ois_data.h5')
dates = datetime(h5read('ois_data.h5', '/dates'));
names = h5read('ois_data.h5', '/varnames');
x  = h5read('ois_data.h5', '/data')';
