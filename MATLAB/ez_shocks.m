% -------------------------------------------------------------------------
% Monetary Policy Shocks in Eurozone
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB;..\..\..\Data')

% Wu & Xia (2017)
load ez_wxrate
wxdstr = str2num(strcat(num2str(shadowrate(:, 1)), num2str(15))); %#ok<ST2NM> % Suppress str2double conversion advice
wxdate = datetime(wxdstr, 'ConvertFrom', 'yyyymmdd', 'Format', 'MMM yyyy');
wxrate = shadowrate(:, 2);

% ----
% Cloyne & Huertgen

% 1. Forecast macro data on assumed information set

% Load factor model data
load ez_data

% Find optimal number of factors according to Bai & Ng (2002)
kmax   = 20; % Maximum number of factors considered
gnum   = 2; % Choose ICp2
demean = 2; % Standardise data before PCA

bnicv = zeros(kmax,1);
for k = 1:kmax
    bnicv(k) = bnic(x, k, gnum, demean); % Compute BNIC for each k
end
rhat = minind(bnicv);

% Extract rhat factors via PCA
[Fhat, Lhat, ef, evf] = factors(x, rhat, demean);

sumeigval = cumsum(evf)/sum(evf);
R2_static = sum(evf(1:rhat))/sum(evf);



% ----------------
% IP Forecast Model
dip = x(:, findstrings(names, {'EKIPTOTG'}));

% Set depvar p = 3 via AIC
pdip = aroptlag(dip, 36, 'aic', 1, 0, 0);

dipvars = [dip, Fhat];
%corrplot(dipvars)

% Static regression
dipmodel = vare(dipvars, pdip);
dipnames = char('dIP', 'F1', 'F2', 'F3', 'F4');
%prt_var(dipmodel, dipnames, fopen('varout_dip.txt', 'w'));
%plt_var(dipmodel, dipnames);

[nobs, ~] = size(dipvars);

hmax = 6; % Maximum forecast horizon
tmin = 80; % Holdout sample
nfor = nobs - tmin - hmax;

yh = zeros(nfor, hmax);
ff = zeros(nfor, hmax);
fe = zeros(nfor, hmax);
for i = 1:nfor
    
    fstart = tmin+i;
    modelforc = varf(dipvars, pdip, hmax, fstart);
    
    % Actual values of variable of interest (1)
    actual = dipvars(fstart:fstart+hmax-1, 1);
    for j = 1:hmax
        yh(i, j) = actual(j);
        ff(i, j) = modelforc(j, 1);
        fe(i, j) = modelforc(j, 1) - actual(j);
    end
end

plot(dates(tmin+1:tmin+nfor), fe(:, 1))


% ----------------
% Inflation Forecast Model
inf = x(:, findstrings(names, {'EMESHARM'}));
%plot(dates, inf)

% Set depvar p = 12 via HQC, AIC not stable.
pinf = aroptlag(inf, 24, 'hqc', 1, 0, 0);

infvars = [inf, Fhat];
%corrplot(infvars)

% Static regression
infmodel = vare(infvars, pinf);
infnames = char('Infl', 'F1', 'F2', 'F3', 'F4');
%prt_var(infmodel, infnames, fopen('varout_inf.txt', 'w'));
%plt_var(infmodel, infnames);

[nobs, ~] = size(infvars);

hmax = 6; % Maximum forecast horizon

% Single forecast experiment
% begf = 180; % beginning forecast period
% fcasts = varf(infvars, pinf, hmax, begf);
% actual = infvars(begf:begf+hmax-1,:);
% 
% plotstart = 160;
% plot([infvars(plotstart:begf-1, 1); fcasts(:, 1)])
% hold on
% plot(infvars(plotstart:begf+hmax-1, 1))

tmin = 80; % Holdout sample
nfor = nobs - tmin - hmax;

yh = zeros(nfor, hmax);
ff = zeros(nfor, hmax);
fe = zeros(nfor, hmax);
for i = 1:nfor
    
    fstart = tmin+i;
    modelforc = varf(infvars, pinf, hmax, fstart);

    % Actual values of variable of interest (1)
    actual = infvars(fstart:fstart+hmax-1, 1);
    for j = 1:hmax
        yh(i, j) = actual(j);
        ff(i, j) = modelforc(j, 1);
        fe(i, j) = modelforc(j, 1) - actual(j);
    end
end

plot(dates(tmin+1:tmin+nfor), fe(:, 1))


% 2. Match forecasts to announcement dates

% ECB MP announcements
[ddata, ~, ~] = xlsread('ez_announce.xlsx');
[~, dname, ~] = xlsread('ez_announce.xlsx', 'B1:D1');
[~, ddate, ~] = xlsread('ez_announce.xlsx', 'A:A');


% Flip ECB announcement series
announceDates = flipud(datetime(ddate(2:end), 'Format', 'yyyy-MM-dd'));
announceRates = flipud(ddata);

anum = length(announceDates);

% Fill NaN for continuous series, i.e. no change
for i = 1:anum
    for j = 1:3
        if isnan(announceRates(i, j))
            announceRates(i, j) = announceRates(i-1, j);
        end
    end
end

% Plot ECB monetary policy corridor by meeting
figure
for i = 1:3
    plot(announceDates(:), announceRates(:, i))
    hold on
end
plot(announceDates(:), zeros(1, anum), 'black') % Add zero line
hold off


sampleStart = matchsample(announceDates, {'2006-07-06'});
sampleDates = announceDates(sampleStart:end-2); % Start later to allow for estimation of initial forecast and match end of factor data
sampleRates = diff(announceRates(sampleStart-1:end-2, 2)); % ECB decisions MRO
snum = length(sampleDates);

% Alter observation 2006-08-31
sampleDates(matchsample(sampleDates, '2006-08-31')) = '2006-09-01';

% Forecast series based on assumed info set at meeting date
hmax = 6;

dipforcs = zeros(snum, hmax);
infforcs = zeros(snum, hmax);
for i = 1:snum
    
    announceMonth = dateshift(sampleDates(i), 'start', 'month');
    announceIndex = matchsample(dates, announceMonth); % Finds index of announcement month in data set
    dipforcvar = varf(dipvars, pdip, hmax, announceIndex);
    infforcvar = varf(infvars, pinf, hmax, announceIndex);
    
    dipforcs(i, :) = dipforcvar(:, 1)';
    infforcs(i, :) = infforcvar(:, 1)';
end

% Forecast revision, set first revision to zero
dipforcrev = [zeros(1, hmax); diff(dipforcs)];
infforcrev = [zeros(1, hmax); diff(infforcs)];


% 3. Build predictor set at every meeting
i3m = x(:, findstrings(names, {'EMIBOR3'})); % Short term interest rate
unp = x(:, findstrings(names, {'EKESTUNPO'})); % UNP logdiffs

% Match monthly data to meetings via announceIndex, i.e. x(1:index-1).
% at meeting month m --> info(month(dates) < m)

i3minfo = zeros(snum, 1);
dipinfo = zeros(snum, 1);
infinfo = zeros(snum, 1);
unpinfo = zeros(snum, 3);
for i = 1:snum
    
    announceMonth = dateshift(sampleDates(i), 'start', 'month');
    announceIndex = matchsample(dates, announceMonth); % Finds index of announcement month in data set
    
    i3minfo(i) = i3m(announceIndex-1);
    dipinfo(i) = dip(announceIndex-1);
    infinfo(i) = inf(announceIndex-1);
    unpinfo(i, :) = [unp(announceIndex-1), unp(announceIndex-2), unp(announceIndex-3)];
end

infoset = [dipforcs, infforcs, dipforcrev, infforcrev,... % Forecasts and revisions
    i3minfo, dipinfo, infinfo, unpinfo]; % Information about general state of economy
% How large is dataset? Feasible with sample size? If not, shorten forecast
% horizon to 4, i.e. current month and one quarter ahead. And/or shorten
% presample.

plot(sampleDates, sampleRates)

announceReg = ols(sampleRates, [ones(snum, 1), infoset]);
prt(announceReg);
plt(announceReg);

announceShocks = announceReg.resid;


% 4. Transform into monthly series

% Reuse code from R&R replication
sampleMonths = dateshift(sampleDates, 'start', 'month');
mta = sampleMonths(1);
mte = sampleMonths(end);
mdates = (mta:calmonths(1):mte)';

mT = length(mdates);

shocks = announceShocks;
mshocks = zeros(mT, 1);
istep = 1;
jstep = 1;
while istep <= mT-1
    while jstep <= snum-1
        if month(mdates(istep)) == month(sampleDates(jstep))  % If months align,
            mshocks(istep) = shocks(jstep);                 % assign value,
            
            if month(mdates(istep)) == month(sampleDates(jstep + 1))          % If also next month aligns,
                mshocks(istep) = mshocks(istep) + shocks(jstep + 1);        % add value to previous,
                jstep = jstep + 1;                                          % move up one meeting,
            end
            
            jstep = jstep + 1;                              % and move up one meeting.
        else                                            % If months do not align,
            mshocks(istep) = 0;                             % assign zero,
        end
        istep = istep + 1;                                  % move up one month.
    end
end
% Hard code final month. How to solve in loop?
mshocks(end) = shocks(end);

% Replace NaN by zeros
mshocks(isnan(mshocks)) = 0;



% Save results as well as alternative shocks for exploration
% ----------------
% From GJ MP Shocks Database

% Babecka Kucharcukova et al. (2016)
[babShocks, ~, ~] = xlsread('gj_shocks_m.xlsx', 'bab');
[~, bdates, ~] = xlsread('gj_shocks_m.xlsx', 'bab', 'A:A');
babDates = datetime(bdates(3:end));

% Neuenkirch (2013)
[neuShocks, ~, ~] = xlsread('gj_shocks_m.xlsx', 'neu');
[~, ndates, ~] = xlsread('gj_shocks_m.xlsx', 'neu', 'A:A');
neuDates = datetime(ndates(3:end));

% ----------------
save ez_shocks -v7.3 mdates mshocks babDates babShocks neuDates neuShocks

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Matrix of independent variables
X = data(sample, 5:22);

% Delete NaN rows
nanindex = any(isnan(X), 2);

X(nanindex, :) = [];
depvar(nanindex, :) = [];

% OLS
beta = ols(depvar, X, 1);
yhat = [ones(length(X), 1), X]* beta;

u = depvar - yhat;







% Set R&R sample from 1969-01-14 to 1996-12-17 (FULL)
ta = datetime('1969-01-14');
te = datetime('1996-12-17');
taindex = find(dates == ta);
teindex = find(dates == te);

sample = taindex:teindex;

[T, N] = size(data(sample, :));

% Plot change of intended target DTARG
plot(dates(sample), data(sample, strcmp(names, {'DTARG'})))



% Compute shock series as residual from equation (1)

% Select variables for estimation
depvar = data(sample, strcmp(names, {'DTARG'}));

% Matrix of independent variables
X = data(sample, 5:22);

% Delete NaN rows
nanindex = any(isnan(X), 2);

X(nanindex, :) = [];
depvar(nanindex, :) = [];

% OLS
beta = ols(depvar, X, 1);
yhat = [ones(length(X), 1), X]* beta;

u = depvar - yhat;

% Reinsert NaN
shocks = zeros(T, 1);
step = 1;

for i = 1:T
    if nanindex(i) == 1
        shocks(i) = NaN;
    else
        shocks(i) = u(step);
        step = step + 1;
    end
end

% Compare to R&R shocks by meeting
figure
plot(dates(sample), data(sample, strcmp(names, {'RESID'})))
hold on
plot(dates(sample), shocks)


% Transform into monthly series and sum if several meetings in same month
mta = dateshift(ta, 'start', 'month');
mte = dateshift(te, 'start', 'month');
mdates = (mta:calmonths(1):mte)';

mT = length(mdates);

mshocks = zeros(mT, 1);
istep = 1;
jstep = 1;

while istep <= mT-1
    
    while jstep <= T-1
        if month(mdates(istep)) == month(dates(jstep))  % If months align,
            mshocks(istep) = shocks(jstep);                 % assign value,
            
            if month(mdates(istep)) == month(dates(jstep + 1))          % If also next month aligns,
                mshocks(istep) = mshocks(istep) + shocks(jstep + 1);        % add value to previous,
                jstep = jstep + 1;                                          % move up one meeting,
            end
            
            jstep = jstep + 1;                              % and move up one meeting.
        else                                            % If months do not align,
            mshocks(istep) = 0;                             % assign zero,
        end
            
        istep = istep + 1;                          % move up one month.
    end   
end
% Hard code final month. How to solve in loop?
mshocks(end) = shocks(end);

% Replace NaN by zeros
mshocks(isnan(mshocks)) = 0;

save rr_shocks -v7.3 mdates mshocks

%%%%
% Historically implemented FFR
ffrpath = 'C:\Dateien\Dropbox\Lennart\Thesis\Data\ffr.xls';

[ffr, ~, ~] = xlsread(ffrpath);
[~, ffrdates, ~] = xlsread(ffrpath, 'A:A');

ffrdates = datetime(ffrdates(12:end));

ffrtaindex = find(ffrdates == ta);
ffrteindex = find(ffrdates == te);

ffrsample = ffrtaindex:ffrteindex;





