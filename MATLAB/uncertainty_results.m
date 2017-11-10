%%%%%%%%
% Results

%%%%
% Compare to JLN aggu results
%clear; clc;

%save jlnresults jlnut utcsa utpca
load jlnresults % aggu
load jlnresults2 % from run code
load ut


Uind = sqrt(ut);
Uavg = squeeze(mean(Uind,2));


mean(Uind);
summarize(Uind);

scatter(1:N, mean(Uind(:,:,1)))




utcsa1 = squeeze(mean(jlnut,2));
utcsa2 = squeeze(mean(jlnut2,2));

lbsum = summarize(Uavg);
jlnsum = summarize(utcsa);

array2table([lbsum, zeros(12,1), jlnsum]);

%%%%
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



h = 12;
vselect = 4;

figure
for i = 1:h
    
    subplot(4, 3, i)
    plot(dates, ut(:, vselect, i))
    hold on
    plot(dates, jlnut2(:, vselect, i))
    title(['h = ',num2str(i)])
end





% Single series uncertainty
hselect = 12;

figure
for i = 1:10
    plot(dates, Uind(:, i, hselect))
    hold on
end



% Aggregate uncertainty, simple average
figure
for i = [1, 3, 12]
    plot(dates, Uavg(:, i))
    hold on
end
legend('show')




%%%%
% Compare alternative aggregate uncertainties
h = 12;


Upcadlogsum = zeros(T, h);

for i = 2:T
    
    Upcadlogsum(i, :) = Upcadlogsum(i-1, :) + Upcadlog(i-1, :);
end


figure
for i = 1:h
    
    subplot(4, 3, i)
    plot(dates, Upcascaled(:, i))
    hold on
    plot(dates, Uavg(:, i))
    %hold on
    %plot(dates(1:end), Upcadlogsum(:, i))
    title(['h = ',num2str(i)])
end

% Plot diffed measures
figure
for i = 1:h
    
    subplot(4, 3, i)
    plot(dates(1:end), standardise(utcsa(:, i)))
    hold on
    plot(dates(1:end), standardise(Uavg(:, i)))
    title(['h = ',num2str(i)])
end


% Does it hold that U is larger for larger h?
tselect = 600;
vselect = 1:9;

figure
for i = vselect
    plot(1:h, squeeze(Uind(tselect, i, :)), 'DisplayName', strcat('Var', num2str(i)))
    hold on
end
legend('show')

legend(names(vselect))









