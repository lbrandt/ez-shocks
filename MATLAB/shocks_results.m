% -------------------------------------------------------------------------
% ez-shocks
% Results
% -------------------------------------------------------------------------

% ----------------
% Uncertainty plot
figure('name', 'Uncertainty in German data');
subplot(2,1,1);
for i = [1, 3, 12]
    plot(dates, Uavg(:, i), 'DisplayName', ['h = ', num2str(i)], ...
    'LineWidth', 1.5)
    hold on
end
legend('show')
title('Uavg')

subplot(2,1,2);
for i = [1, 3, 12]
    plot(dates, Ufac(:, i), 'DisplayName', ['h = ', num2str(i)], ...
    'LineWidth', 1.5)
    hold on
end
legend('show')
title('Ufac')

% Print to pdf
print('de_uncertainty_plot', '-dpdf', '-bestfit')


% ----------------
% IRF plot

varirf = [bvarirf, nvarirf, wvarirf, wzarirf];
jorirf = [bjorda.irf, njorda.irf, wjorda.irf];

birf = [bvarirf, bjorda.irf];
nirf = [nvarirf, njorda.irf];
wirf = [wvarirf, wjorda.irf];



figure('name', 'Euro Area MP Shocks');
subplot(3, 1, 1);
plot(wirf, 'LineWidth', 1.5)
title('Shadow rate -->  German IP')
subplot(3, 1, 2)
plot(nirf, 'LineWidth', 1.5)
title('Communication --> German IP')
subplot(3, 1, 3)
plot(birf, 'LineWidth', 1.5)
title('MCI --> German IP')

% Print to pdf
print('de_irf_plot', '-dpdf', '-bestfit')

