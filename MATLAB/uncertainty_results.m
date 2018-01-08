%%%%%%%%
% Results

%%%%
% Compare to JLN aggu results
%clear; clc;

load jlnresults % aggu
load uncertainty

load de_uncertainty
load gs_uncertainty




utcsa1 = squeeze(mean(jlnut,2));
%isequal(utcsa, utcsa1) % Why different? utcsa should be "correct" though.

%%%%
% Descriptive statistics
%jlndesc = summarize(utcsa);
%usdesc1 = summarize(Uavg);
%usdesc2 = summarize(Ufac);

%array2table([jlndesc, usdesc1, usdesc2], 'VariableNames', {'N', 'mean', 'sd', 'min', 'max'});


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
%print('gs_uncertainty_plot', '-depsc')



% Single series uncertainty
hselect = 12;
vselect = 2;

figure
for i = 1:hselect
    plot(dates, Uind(:, vselect, i))
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
tselect = 250;
vselect = 1:9;

figure
for i = vselect
    plot(1:hselect, squeeze(Uind(tselect, i, :)), 'DisplayName', strcat('Var', num2str(i)))
    hold on
end
legend('show')

legend(names(vselect))
