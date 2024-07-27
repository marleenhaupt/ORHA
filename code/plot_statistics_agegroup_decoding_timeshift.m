function [] = plot_statistics_agegroup_decoding_timeshift(result)

data=result';
n_perm=100; %originally 10000  
significance_thr = 0.05;
output_path = '../output/';

%% Statistics
[SignificantVariables_all,~,adjusted_pvalues_all] = fdr_corrected_perm_test_1sample(data, n_perm, significance_thr, 'right');

ya=-15;
significant_time_points_all = find(SignificantVariables_all>0);
y_significants_all = repmat(ya, size(significant_time_points_all,2),1)';

%% Plotting
data = data-50;
x = 1:numel(mean(data,1));
x2 = [x, fliplr(x)];
c1 = [0/255 0/255 0/255];
SEM_all = std(data)/sqrt(size(data,1));

figure
plot(mean(data),'Color',c1, 'LineWidth', 1.6)
hold on
upper = mean(data) + SEM_all;
lower = mean(data) - SEM_all;
inBetween = [upper, fliplr(lower)];
fill(x2, inBetween, c1, 'FaceAlpha', 0.16, 'LineStyle', 'none');
hold on;
plot(significant_time_points_all, y_significants_all,'.','Color',c1)
hold on
title('Age group decoding accounting for timeshift')
xlabel('time (ms)')
ylabel('classification accuracy - chance level (%)')
xticks([20 40 60])
set(gca, 'XTickLabel', [100 200 300])
yline(0,'color', 'black' ,'LineStyle','--', 'LineWidth', 2);
xline(33, 'color', 'red', 'LineStyle','-', 'LineWidth', 2);
xlim([20,60])
ylim([-20,50])
set(gca,'box','off')

%% save
saveas(gca,sprintf('%sAgegroup_decoding_timeshift.png',output_path));
end

