function [] = plot_statistics_decoding(matrix_avg_allsubs,decodinglevel)

% run cluster based permutation test
n_perm = 100; %originally 10000
cluster_thr = 0.05;
significance_thr = 0.05;
output_path = '../output/';

[SignificantVariables_young,significantVarMax_young,pValWei_young,pValMax_young,clusters_young] = permutation_cluster_1sample(matrix_avg_allsubs(1:21,:),n_perm,cluster_thr,significance_thr,'right');
[SignificantVariables_old,significantVarMax_old,pValWei_old,pValMax_old,clusters_old] = permutation_cluster_1sample(matrix_avg_allsubs(22:43,:), n_perm, cluster_thr, significance_thr,'right');
[SignificantVariables_agediff,significantVarMax_agediff,pValWei_agediff,pValMax_agediff,clusters_agediff] = permutation_cluster_2sample(matrix_avg_allsubs(1:21,:),matrix_avg_allsubs(22:43,:), n_perm, cluster_thr, significance_thr,'right');

% find significant time points for plot
yb=-6;
yw=-7;
ys=-8;
    
significant_time_points_young = find(SignificantVariables_young>0);
y_significants_young = repmat(yb, size(significant_time_points_young,2),1)';

significant_time_points_old = find(SignificantVariables_old>0);
y_significants_old = repmat(yw, size(significant_time_points_old,2),1)';

significant_time_points_agediff = find(SignificantVariables_agediff>0);
y_significants_agediff = repmat(ys, size(significant_time_points_agediff,2),1)';

% boostrap the peaks
bootstrap_samples = 100; %originally 10000
    
%for younger
peak_latency_samples_y = NaN(bootstrap_samples,2);

rng('shuffle')
for bs = 1:bootstrap_samples
	bootstrapped_datasets = datasample(matrix_avg_allsubs(1:21,:),length(matrix_avg_allsubs(1:21,1)),1);
	avg_datasets = squeeze(mean(bootstrapped_datasets,1));
	[peaks_x, peaks_y] = find(avg_datasets==max(avg_datasets,[],'all'));
	if size(peaks_x,2) > 1
        peaks_x = round(mean(peaks_x));
        peaks_y = round(mean(peaks_y));
    end
    peak_latency_samples_y(bs,1)=peaks_x;
    peak_latency_samples_y(bs,2)=peaks_y;
end

peak_latency_y = round(mean(peak_latency_samples_y,1));
CI_low_y = prctile(peak_latency_samples_y,2.5);
CI_high_y = prctile(peak_latency_samples_y,97.5);
    
%for older
peak_latency_samples_o = NaN(bootstrap_samples,2);

rng('shuffle')
for bs = 1:bootstrap_samples
	bootstrapped_datasets = datasample(matrix_avg_allsubs(22:43,:),length(matrix_avg_allsubs(22:43,1)),1);
	avg_datasets = squeeze(mean(bootstrapped_datasets,1));
	[peaks_x, peaks_y] = find(avg_datasets==max(avg_datasets,[],'all'));
	if size(peaks_x,2) > 1
        peaks_x = round(mean(peaks_x));
        peaks_y = round(mean(peaks_y));
    end
	peak_latency_samples_o(bs,1)=peaks_x;
	peak_latency_samples_o(bs,2)=peaks_y;
end
peak_latency_o = round(mean(peak_latency_samples_o,1));
CI_low_o = prctile(peak_latency_samples_o,2.5);
CI_high_o = prctile(peak_latency_samples_o,97.5);
    
%age group difference
peak_latency_samples_d=peak_latency_samples_o-peak_latency_samples_y;
peak_latency_d = round(mean(peak_latency_samples_d,1));
CI_low_d = prctile(peak_latency_samples_d,2.5);
CI_high_d = prctile(peak_latency_samples_d,97.5);
    
diff_peak(1,1) = peak_latency_d(1,2);
diff_peak(2,2) = CI_low_d(1,2);
diff_peak(3,3) = CI_high_d(1,2);

pvals = sum((peak_latency_samples_d(:,2)<=0)) / (size(peak_latency_samples_d,1));

x = 1:numel(mean(matrix_avg_allsubs,1));
x2 = [x, fliplr(x)];
%c1 = [0/255 0/255 0/255];
c2= [62/255 73/255 137/255];
c3 = [181/255 222/255 43/255];
c4 = [128/255 128/255 128/255];

%substract chance level
matrix_avg_allsubs = matrix_avg_allsubs-50;

SEM_young = std(matrix_avg_allsubs(1:21,:))/sqrt(size(matrix_avg_allsubs(1:21,:),1));
SEM_old = std(matrix_avg_allsubs(22:43,:))/sqrt(size(matrix_avg_allsubs(22:43,:),1));

figure
a=plot(mean(matrix_avg_allsubs(1:21,:)),'Color',c2, 'LineWidth', 2)
hold on
b=plot(mean(matrix_avg_allsubs(22:43,:)),'Color',c3,'LineWidth', 2)
hold on
c=plot(mean(matrix_avg_allsubs(1:21,:))-mean(matrix_avg_allsubs(22:43,:)),'Color',c4, 'LineWidth', 2)
hold on
upper = mean(matrix_avg_allsubs(1:21,:)) + SEM_young;
lower = mean(matrix_avg_allsubs(1:21,:)) - SEM_young;
inBetween = [upper, fliplr(lower)];
fill(x2, inBetween, c2, 'FaceAlpha', 0.155, 'LineStyle', 'none');
hold on;
upper = mean(matrix_avg_allsubs(22:43,:)) + SEM_old;
lower = mean(matrix_avg_allsubs(22:43,:)) - SEM_old;
inBetween = [upper, fliplr(lower)];
fill(x2, inBetween, c3, 'FaceAlpha', 0.15, 'LineStyle', 'none');
hold on;
plot(significant_time_points_agediff, y_significants_agediff,'.','Color',c4)
hold on
plot(significant_time_points_young, y_significants_young,'.', 'Color',c2)
hold on
plot(significant_time_points_old, y_significants_old,'.', 'Color',c3)
hold on
plot([CI_low_y(1,2);CI_high_y(1,2)],[49;49],'Color',c2, 'LineWidth', 2) 
hold on
plot(peak_latency_y(1,2),49,'.','MarkerSize',10,'Color',c2)
hold on
plot([CI_low_o(1,2);CI_high_o(1,2)],[48;48],'Color',c3, 'LineWidth', 2) 
hold on
plot(peak_latency_o(1,2),48,'.','MarkerSize',10,'Color',c3)
hold on
title(sprintf("%s decoding in time", decodinglevel))
xlabel('time (ms)')
ylabel('classification accuracy - chance level (%)')
xticks([0 40 80 120 160 200 240])
set(gca, 'XTickLabel', [-200 0 200 400 600 800 1000])
yline(0,'color', '#808080' ,'LineStyle','--', 'LineWidth', 2);
xline(40, 'color', '#808080', 'LineStyle','--', 'LineWidth', 2);
xlim([0,240])
ylim([-10,50])
legend([a b c], 'younger', 'older','diff')
set(gca,'box','off')
legend('boxoff')
saveas(gca,sprintf('%s%s_decoding_time.png',output_path, decodinglevel));

end


