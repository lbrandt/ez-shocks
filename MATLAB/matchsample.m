function [sample1, sample2, ta, te] = matchsample(dates1, dates2)
% -------------------------------------------------------------------------
% Finds the largest set of shared observations between two time series of
% the same frequency but differing range.
%
%   Input:  dates1  = Vector of dates [T x 1]
%           dates2  = Vector of dates [U x 1]
%   Output: sample1 = Vector of common observations expressed in indices of
%                     first time series
%           sample2 = Vector of common observations expressed in indices of
%                     second time series
%           ta      = Datetime of first common observation
%           te      = Datetime of last common observation
% -------------------------------------------------------------------------

% Coerce inputs into dates
d1 = datetime(dates1);
d2 = datetime(dates2);

% Inquire range of time series
ta1 = d1(1);
ta2 = d2(1);
te1 = d1(end);
te2 = d2(end);

% Set range of shared sample
if ta1 <= ta2
    ta = ta2;
elseif ta1 > ta2
    ta = ta1;
end

if te1 <= te2
    te = te1;
elseif te1 > te2
    te = te2;
end

% Construct overlap
taindex1 = find(dates1 == ta);
teindex1 = find(dates1 == te);

sample1  = taindex1:teindex1;

taindex2 = find(dates2 == ta);
teindex2 = find(dates2 == te);

sample2  = taindex2:teindex2;
end
