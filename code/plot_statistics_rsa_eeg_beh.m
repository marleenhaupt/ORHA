function [] = plot_statistics_rsa_eeg_beh(matrix_allsub)

n_perm=100; %originally 10000
cluster_thr = 0.05;
significance_thr = 0.05;
output_path = '../output/';


%% Statistics
%statistics code subtracts chance level; as these are corr values, we
%have to add 50 here
rdm_rsa = matrix_allsub+50;

[SignificantVariables_young,significantVarMax_young,pValWei_young,pValMax_young,clusters_young] = permutation_cluster_1sample(rdm_rsa(1:21,:), n_perm, cluster_thr, significance_thr,'right');
[SignificantVariables_old,significantVarMax_old,pValWei_old,pValMax_old,clusters_old] = permutation_cluster_1sample(rdm_rsa(22:43,:), n_perm, cluster_thr, significance_thr,'right');
[SignificantVariables_agediff,significantVarMax_agediff,pValWei_agediff,pValMax_agediff,clusters_agediff] = permutation_cluster_2sample(rdm_rsa(1:21,:),rdm_rsa(22:43,:),n_perm, cluster_thr, significance_thr,'both');

%plotting specifications
y_a= -0.06;
y_y= -0.05;
y_o= -0.052;
y_d= -0.054;
    
significant_time_points_y = find(SignificantVariables_young>0);
y_significants_y = repmat(y_y, size(significant_time_points_y,2),1)';
    
significant_time_points_o = find(SignificantVariables_old>0);
y_significants_o = repmat(y_o, size(significant_time_points_o,2),1)';

significant_time_points_d = find(SignificantVariables_agediff>0);
y_significants_d = repmat(y_d, size(significant_time_points_d,2),1)';

SEM_y   = std(matrix_allsub(1:21,:))./sqrt(size(matrix_allsub(1:21,:),1));
SEM_o   = std(matrix_allsub(22:43,:))./sqrt(size(matrix_allsub(22:43,:),1));
    
%% Bootstrapping peak

bootstrap_samples = 100; %originally 10000


%young
peak_latency_samples_y = NaN(bootstrap_samples,2);

rng('shuffle')
for bs = 1:bootstrap_samples
	bootstrapped_datasets = datasample(matrix_allsub(1:21,:),length(matrix_allsub(1:21,1)),1);
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
    
%old
peak_latency_samples_o = NaN(bootstrap_samples,2);

rng('shuffle')
for bs = 1:bootstrap_samples
	bootstrapped_datasets = datasample(matrix_allsub(22:42,:),length(matrix_allsub(22:42,1)),1);
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
    
%diff
peak_latency_samples_d=peak_latency_samples_o-peak_latency_samples_y;
peak_latency_d = round(mean(peak_latency_samples_d,1));
CI_low_d = prctile(peak_latency_samples_d,2.5);
CI_high_d = prctile(peak_latency_samples_d,97.5);
    
diff_peak(1,1) = peak_latency_d(1,2);
diff_peak(2,2) = CI_low_d(1,2);
diff_peak(3,3) = CI_high_d(1,2);
    
% prerequisites for plotting
%c1 = [0/255 0/255 0/255];
c2= [62/255 73/255 137/255];
c3 = [181/255 222/255 43/255];
c4 = [128/255 128/255 128/255];
x_y = 1:numel(mean(matrix_allsub(1:21,:),1));
x2_y = [x_y, fliplr(x_y)];
x_o = 1:numel(mean(matrix_allsub(22:42,:),1));
x2_o = [x_o, fliplr(x_o)];

% plot
figure
young = plot(mean(matrix_allsub(1:21,:)),'Color',c2, 'LineWidth', 1.6)
hold on
old = plot(mean(matrix_allsub(22:42,:)),'Color',c3, 'LineWidth', 1.6)
hold on
diff = plot(mean(matrix_allsub(1:21,:))-mean(matrix_allsub(22:42,:)),'Color',c4, 'LineWidth', 1.6)
hold on
plot(significant_time_points_y, y_significants_y,'.','Color',c2)
hold on
plot(significant_time_points_o, y_significants_o,'.','Color',c3)
hold on
plot(significant_time_points_d, y_significants_d,'.','Color',c4)
hold on
%plot([CI_low_y(1,2);CI_high_y(1,2)],[0.08;0.08],'Color',c2, 'LineWidth', 1.6) 
%hold on
%plot(peak_latency_y(1,2),0.08,'.','MarkerSize',10,'Color',c2)
%hold on
%plot([CI_low_o(1,2);CI_high_o(1,2)],[0.07;0.07],'Color',c3, 'LineWidth', 1.6) 
%hold on
%plot(peak_latency_o(1,2),0.07,'.','MarkerSize',10,'Color',c3)
%hold on
upper = mean(matrix_allsub(1:21,:)) + SEM_y;
lower = mean(matrix_allsub(1:21,:)) - SEM_y;
inBetween = [upper, fliplr(lower)];
fill(x2_y, inBetween, c2, 'FaceAlpha', 0.155, 'LineStyle', 'none');
hold on;
upper = mean(matrix_allsub(22:42,:)) + SEM_o;
lower = mean(matrix_allsub(22:42,:)) - SEM_o;
inBetween = [upper, fliplr(lower)];
fill(x2_o, inBetween, c3, 'FaceAlpha', 0.155, 'LineStyle', 'none');
hold on;
xlabel('time (ms)')
ylabel('Spearman correlation coefficient')
xticks([20 40 60 80 100 120 140 160 180 200 220 240])
set(gca, 'XTickLabel', [-100 0 100 200 300 400 500 600 700 800 900 1000])
yline(0,'LineStyle','--', 'LineWidth', 1.5);
xline(40, 'LineStyle','--', 'LineWidth', 1.5);
ylim([-0.06,0.06])
set(gca,'box','off')
legend([young old diff],'younger', 'older','diff')
legend('boxoff')

%% save
saveas(gca,sprintf('%sRSA_EEG_beh.png',output_path));
