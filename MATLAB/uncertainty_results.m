%%%%%%%%
% Results

%%%%
% Compare to JLN aggu results
%clear; clc;

load jlnresults % aggu
load uncertainty




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
    plot(dates, Uavg(:, i))
    hold on
end
legend('show')

subplot(2,1,2);
for i = [1, 3, 12]
    plot(dates, utcsa(:, i))
    hold on
end
legend('show')




% Single series uncertainty
hselect = 12;
vselect = 1;

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
tselect = 600;
vselect = 1:9;

figure
for i = vselect
    plot(1:hselect, squeeze(Uind(tselect, i, :)), 'DisplayName', strcat('Var', num2str(i)))
    hold on
end
legend('show')

legend(names(vselect))
