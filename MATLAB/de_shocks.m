% -------------------------------------------------------------------------
% Generate monetary policy shocks for Germany
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB;..\..\..\Data')

% Load de factor model data
load de_data % Agent's observable information set

% Load ECB monetary policy decisions
[ddata, ~, ~] = xlsread('ez_announce.xlsx');
[~, dname, ~] = xlsread('ez_announce.xlsx', 'B1:D1');
[~, ddate, ~] = xlsread('ez_announce.xlsx', 'A:A');

announceDates = datetime(ddate(2:end));

% Flip ECB announcement series
ydate = flipud(announceDates);

ta = datetime('1999-01-04');
te = datetime('2017-12-14');
taindex = find(ydate == ta);
teindex = find(ydate == te);

announceFull = taindex:teindex;

dydate = ydate(2) - ydate(1);

y = flipud(ddata);
% Fill NaN for continuous series
for i = 1:length(ydate)
    for j = 1:3
        if isnan(y(i, j))
            y(i, j) = y(i-1, j);
        end
    end
end

% Euro Wu & Xia shadow rates as actual monetary policy stance at ZLB
load ez_wxrate
wxdstr = str2num(strcat(num2str(shadowrate(:, 1)), num2str(15))); %#ok<ST2NM> % Suppress str2double conversion advice
wxdate = datetime(wxdstr(1:end-1), 'ConvertFrom', 'yyyymmdd', 'Format', 'yyyy-MM-dd');
wxrate = shadowrate(1:end-1, 2);
wxdiff = diff(wxrate, 1);

wta = wxdate(2);
wte = wxdate(end);
wtaindex = find(dates == wta);
wteindex = find(dates == wte);

wsample = wtaindex:wteindex;
wxDates = wxdate(2:end);


% ----------------
% Choose variable of interest from factor model dataset
chooseVariable = {'eonia', 'BDUN_TOTQ', 'BDCONPRCE'};
eonia = x(1:end, strcmp(varnames, {'eonia'}));
is3mo = x(1:end, strcmp(varnames, {'is3mo'}));

plot(dates, is3mo)
hold on
plot(dates, [zeros(wtaindex-1, 1); wxdiff])

% ----------------
% Short sample including shadow rate
xwx = [wxdiff, x(wsample, :)];

[T, N] = size(xwx);


% PCA

% Find optimal number of factors according to Bai & Ng (2002)
kmax   = 20; % Max number of factors to be extracted
gnum   = 1; % ICp2 chosen in JLN2015
demean = 2; % Standardise data

bnicv = zeros(kmax,1);
for k = 1:kmax
    bnicv(k) = bnic(xwx, k, gnum, demean); % Compute BNIC for each k
end

bnicmin = min(bnicv);
rhat = minind(bnicv); % Optimal number of factors according to lowest IC
fprintf('\nFactors via IC(%d): rhat = %d \n', gnum, rhat);


%%%% Override BNIC
%rhat = 12
%%%%


% Extract estimated number of factors

[Fhat, LF, ef, evf] = factors(xwx, rhat, demean);
[Ghat, LG, eg, evg] = factors(xwx.^2, rhat, demean);

sumeigval = cumsum(evf)/sum(evf);
R2_static = sum(evf(1:rhat))/sum(evf);

% ----------------
% In order to "name" one factor the MP factor, check correlation with shadow rate:
F = xwx(:, 1);

% Matrix of correlations between true factors and estimated factors
[~, r] = size(F);

fcorrmat = zeros(r, rhat);
frmse    = zeros(r, rhat);
for i = 1:r
    for j = 1:rhat
        corrcoeffs     = corrcoef(F(:, i), Fhat(:, j));
        fcorrmat(i, j) = corrcoeffs(2, 1); % Extract correlation
        
        if fcorrmat(i, j) <= 0
            Fhat(:, j) = -Fhat(:, j); % Flip
        end
        
        frmse(i, j)    = rmse(standardise(F(:, i)), standardise(Fhat(:, j)));
    end
end

figure
subplotindex = 1;
for i = 1:r
    for j = 1:rhat
        subplot(rhat, r, subplotindex)
        
        plot(standardise(F(:, i)))
        hold on
        plot(standardise(Fhat(:, j)))
        title(['F(', num2str(i), '), Fhat(', num2str(j),'); rho = ', num2str(fcorrmat(i, j)), ', RMSE = ', num2str(frmse(i, j))])
        
        subplotindex = subplotindex + 1;
    end
end

% ----------------
% Since shadow rate is highly correlated with first factor, restrict Lambda
% according to Stock & Watson (2016) p. 474

% Check if cov(Fhat) is full rank, i.e. r = q
qhat = rank(cov(Fhat));

%Lrestrict = zeros(N, rhat);
Lrestrict = LF; %XYXY Not true. Need to minimise LS

Lrestrict(1, :) = [1, 0, 0]; % Normalise first row

Lrestrict(2:end, 1) = Lrestrict(2:end, 1)./LF(1, 1); %XYXY "Take loadings on F1 from variables"
%Frestrict = Fhat*LF(1, 1); %XYXY "Give loadings back to factors"
%test = Frestrict*Lrestrict';
test = Fhat*Lrestrict';

plot(xwx(:, 1))
hold on
plot(test(:, 1))


%%%%%%%%%%%%%%%%
% Generate shocks
nlag  = fix(4*(T/100)^(2/9));

% AR model
ar = [ones(T, 1), mlag(xwx(:, 1), 12)];

wareg = nwest(xwx(:, 1), ar, nlag);

wares = wareg.resid;

% Fit
plot(wareg.y)
hold on
plot(wareg.yhat)

% Errors
plot(wareg.resid)



% Linear predictors
xt = [ones(T, 1), mlag(xwx(:, 1), 1), mlag(Fhat, 1)];

% Predict shadow rate
%wxreg = ols(xwx(:, 1), xt);
wxreg = nwest(xwx(:, 1), xt, nlag);

wxres = wxreg.resid;

% Fit
plot(wxreg.y)
hold on
plot(wxreg.yhat)

% Errors
plot(wxreg.resid)


% Nonlinear predictors
z      = [Fhat, Fhat(:, 1).^2, Ghat(:, 1)];
zt     = [ones(T, 1), mlag(z, 1)];
[~, M] = size(zt);

% Predict shadow rate
%wxregz = ols(xwx(:, 1), zt);
wzreg = nwest(xwx(:, 1), zt, nlag);

wzres = wzreg.resid;

% Fit
plot(wzreg.y)
hold on
plot(wzreg.yhat)

% Errors
plot(wxres)
hold on
plot(wares)

summarize(wxres);
summarize(wares);

histogram(wares)

save de_shocks -v7.3 wxDates wxres wzres
