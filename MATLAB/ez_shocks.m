% -------------------------------------------------------------------------
% Euro OIS
% Monetary Policy Shocks
% -------------------------------------------------------------------------

%clear; clc;
%addpath('..\R;..\MATLAB;..\..\..\Data')

% Load data
[ddata, ~, ~] = xlsread('ez_announce.xlsx');

[~, dname, ~] = xlsread('ez_announce.xlsx', 'B1:D1');
[~, ddate, ~] = xlsread('ez_announce.xlsx', 'A:A');

ddate = datetime(ddate(2:end));

%h5disp('ois_data.h5')
dates = datetime(h5read('ois_data.h5', '/dates'));
names = h5read('ois_data.h5', '/varnames');
data  = h5read('ois_data.h5', '/data')';


% % Wu & Xia (2017)
% load ez_wxrate
% wxdstr = str2num(strcat(num2str(shadowrate(:, 1)), num2str(15))); %#ok<ST2NM> % Suppress str2double conversion advice
% wxdate = datetime(wxdstr, 'ConvertFrom', 'yyyymmdd', 'Format', 'MMM yyyy');
% wxrate = shadowrate(:, 2);



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





