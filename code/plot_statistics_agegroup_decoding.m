function [] = plot_statistics_agegroup_decoding(data)

data=data(:,:)'; 
n_perm=100; %originally 10000  
significance_thr = 0.05;
output_path = '../output/';
    
%% Statistics
[SignificantVariables_all,~,adjusted_pvalues_all] = fdr_corrected_perm_test_1sample(data, n_perm, significance_thr,'right');

ya=-20;
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
title('Age group decoding based on EEG RDMs')
xlabel('time (ms)')
ylabel('classification accuracy - chance level (%)')
xticks([0 40 80 120 160 200 240])
set(gca, 'XTickLabel', [-200 0 200 400 600 800 1000])
yline(0,'color', '#808080' ,'LineStyle','--', 'LineWidth', 2);
ylim([-25,50])

%% save
saveas(gca,sprintf('%sAgegroup_decoding.png',output_path));

end

