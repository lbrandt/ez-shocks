%%%%%%%%
% Results

%%%%
% Compare to JLN aggu results
%clear; clc;

load jlnresults % aggu
load uncertainty

load ez_uncertainty
load de_uncertainty
load gs_uncertainty

%%%%
% Descriptive statistics
ulognorm = log(uind.^2); % Normal

ulognormavg = squeeze(mean(ulognorm,2)); % Normal
summarize(ulognormavg);
histogram(ulognormavg(:, 1))

expuavg = exp(ulognormavg); % LogNormal
summarize(expuavg);

expuavg2 = expuavg.^2; % ??
summarize(expuavg2);


upca2 = factors(uind(:, :, 1), 1, 2);

summarize(upca2);


[T, N, h] = size(uind(:, :, :));
% Correct aggu?
logu = log(uind.^2); % Normal
loguavg = squeeze(mean(logu,2)); % Normal
uavg2 = exp(loguavg);
uavg3 = sqrt(uavg2);


histogram(uavg2(:, 1))
[jbdecision, ~,jbstat, jbcrit] = jbtest(uavg2(:, 1))



ufac1 = zeros(T, h); % Normal
ufac2 = zeros(T, h); % LogNormal
ufac3 = zeros(T, h); % ?
upca = zeros(T, h); % Normal
for i = 1:h
    upca(:, i) = factors(logu(:, :, i), 1, 2); % Normal?
    
    % Flip ufac if necessary
    rho = corrcoef(upca(:, i), loguavg(:, i));
    if rho(2, 1) < 0
        upca = -upca;
    end
    
    % Scale to Uavg
    ufac1(:, i) = standardise(upca(:, i))* std(loguavg(:, i)) + mean(loguavg(:, i));
    ufac2(:, i) = exp(ufac1(:, i));
    ufac3(:, i) = sqrt(ufac2(:, i));
end

summarize(uavg);

uavg4 = exp(0.5*loguavg);

plot(standardise(uavg(:, 1)))
hold on
plot(loguavg(:, 1)) % Normal
hold on
plot(uavg2(:, 1))
hold on
plot(uavg3(:, 1))
hold on
plot(uavg4(:, 1))

plot(upca(:, 1));
hold on
plot(ufac1(:, 1))
hold on
plot(ufac2(:, 1))
hold on
plot(ufac3(:, 1))

plot(loguavg(:, 1)) % Normal
hold on
plot(ufac1(:, 1))

plot(uavg2(:, 1)) % LogNormal
hold on
plot(ufac2(:, 1))

plot(uavg3(:, 1)) % ??
hold on
plot(ufac3(:, 1))




%%%%
% Diagram from JLN 2015
figure
subplot(2,1,1);
for i = [1, 3, 12]
    plot(dates, Uavg(:, i), 'DisplayName', ['h = ',num2str(i)])
    hold on
end
legend('show')
title('Uavg')

subplot(2,1,2);
for i = [1, 3, 12]
    plot(dates, Ufac(:, i), 'DisplayName', ['h = ',num2str(i)])
    hold on
end
legend('show')
title('Ufac')

% Print to pdf
%print('gs_uncertainty_plot', '-dpdf', '-bestfit')
% Print to .eps
%print('gs_ufac_plot', '-depsc')



% Single series uncertainty
hselect = 12;
vselect = 2;

figure
for i = 1:hselect
    plot(dates, uind(:, vselect, i))
    hold on
end



%%%%
% Compare alternative aggregate uncertainties

figure
for i = 1:hselect
    
    subplot(4, 3, i)
    plot(dates, Ufac(:, i))
    hold on
    plot(dates, Uavg(:, i))
    
    title(['h = ',num2str(i)])
end
legend('show')


% Plot standardised measures
figure
for i = 1:hselect
    
    subplot(4, 3, i)
    plot(dates(1:end), standardise(Ufac(:, i)))
    hold on
    plot(dates(1:end), standardise(Uavg(:, i)))
    title(['h = ',num2str(i)])
end


% Does it hold that U is larger for larger h?
tselect = 150;
vselect = 1:9;

figure
for i = vselect
    plot(1:hselect, squeeze(Uind(tselect, i, :)), 'DisplayName', strcat('Var', num2str(i)))
    hold on
end
legend('show')

legend(names(vselect))
