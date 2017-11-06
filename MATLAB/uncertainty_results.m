

%%%%%%%%
% Results

test1 = sqrt(ut(:, 1, 3));
test2 = test1.^2;

isequal(test1, Uind(:, 1, 3))
isequal(test2, ut(:, 1, 3))




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


%%%%
% Compare to JLN aggu results
%clear; clc;

%save jlnresults jlnut utcsa utpca
load jlnresults % aggu
load jlnresults2 % from run code
load ut


Uind = sqrt(ut);
Uavg = squeeze(mean(Uind,2));

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


