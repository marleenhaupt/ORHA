function [] = plot_statistics_fusion(matrix_allsub)

roi_labels = {'V1','V2','V3','FFA','PPA','LOC'};
n_perm=100; %originally 10000
cluster_thr = 0.05;
significance_thr = 0.05;
output_path = '../output/';

figure

% roi loop      
for iroi = 1:size(matrix_allsub,3) %roi dimension  
    roi = string(roi_labels(1,iroi));

    %% Statistics
    
    %statistics code subtracts chance level; as these are corr values, we
    %have to add 50 here
    rsa=squeeze(matrix_allsub(:,:,iroi))+50;
    
    [SignificantVariables_young,significantVarMax_young,pValWei_young,pValMax_young,clusters_young] = permutation_cluster_1sample(rsa(1:21,:), n_perm, cluster_thr, significance_thr,'right');
    [SignificantVariables_old,significantVarMax_old,pValWei_old,pValMax_old,clusters_old] = permutation_cluster_1sample(rsa(22:42,:), n_perm, cluster_thr, significance_thr,'right');
    [SignificantVariables_agediff,significantVarMax_agediff,pValWei_agediff,pValMax_agediff,clusters_agediff] = permutation_cluster_2sample(rsa(1:21,:),rsa(22:42,:), n_perm, cluster_thr, significance_thr,'right');

    %plotting specifications
    y_y= -0.05;
    y_o= -0.06;
    y_d= -0.07;

    significant_time_points_y = find(SignificantVariables_young>0);
    y_significants_y = repmat(y_y, size(significant_time_points_y,2),1)';

    significant_time_points_o = find(SignificantVariables_old>0);
    y_significants_o = repmat(y_o, size(significant_time_points_o,2),1)';

    significant_time_points_d = find(SignificantVariables_agediff>0);
    y_significants_d = repmat(y_d, size(significant_time_points_d,2),1)';
  
    SEM_y   = std(matrix_allsub(1:21,:,iroi))./sqrt(size(matrix_allsub(1:21,:,iroi),1));
    SEM_o   = std(matrix_allsub(22:42,:,iroi))./sqrt(size(matrix_allsub(22:42,:,iroi),1));
    
    %% Bootstrapping peak
    bootstrap_samples = 100; %originally 10000
    
    peak_latency_samples_y = NaN(bootstrap_samples,2);
    peak_latency_samples_o = NaN(bootstrap_samples,2);

    rng('shuffle');
    for bs = 1:bootstrap_samples
        bootstrapped_datasets_y = datasample(rsa(1:21,:),size(rsa(1:21,:),1),1);
        bootstrapped_datasets_o = datasample(rsa(22:42,:),size(rsa(22:42,:),1),1);
        avg_datasets_y = squeeze(mean(bootstrapped_datasets_y,1));
        avg_datasets_o = squeeze(mean(bootstrapped_datasets_o,1));
        [peaks_x_y, peaks_y_y] = find(avg_datasets_y==max(avg_datasets_y,[],'all'));
        [peaks_x_o, peaks_y_o] = find(avg_datasets_o==max(avg_datasets_o,[],'all'));
        if size(peaks_x_y,2) > 1
        peaks_x_y = round(mean(peaks_x_y));
        peaks_y_y = round(mean(peaks_y_y));
        end
        if size(peaks_x_y,2) > 1
        peaks_x_o = round(mean(peaks_x_o));
        peaks_y_o = round(mean(peaks_y_o));
        end
        peak_latency_samples_y(bs,1)=peaks_x_y;
        peak_latency_samples_y(bs,2)=peaks_y_y;
        peak_latency_samples_o(bs,1)=peaks_x_o;
        peak_latency_samples_o(bs,2)=peaks_y_o;
        peak_latency_samples_d(bs,1)=peaks_x_y-peaks_x_o;
        peak_latency_samples_d(bs,2)=peaks_y_y-peaks_y_o;
        peak_latency_perroi_d(bs,iroi)=peaks_y_o-peaks_y_y; 
        peak_latency_perroi_y(bs,iroi)=peaks_y_y;
        peak_latency_perroi_o(bs,iroi)=peaks_y_o;
    end
    peak_latency_y = round(mean(peak_latency_samples_y,1));
    peak_latency_o = round(mean(peak_latency_samples_o,1));
    peak_latency_d = round(mean(peak_latency_samples_d,1));

    CI_low_y = prctile(peak_latency_samples_y,2.5);
    CI_high_y = prctile(peak_latency_samples_y,97.5);
    CI_low_o = prctile(peak_latency_samples_o,2.5);
	CI_high_o = prctile(peak_latency_samples_o,97.5);
    CI_low_d = prctile(peak_latency_samples_d,2.5);
    CI_high_d = prctile(peak_latency_samples_d,97.5);
    
    pvals_new(iroi) = sum((peak_latency_samples_d(:,2)>=0)) / (size(peak_latency_samples_d,1));
    delay(iroi) = peak_latency_d(2);

    
    % prerequisites for plotting
    %c1 = [0/255 0/255 0/255];
    c2= [62/255 73/255 137/255];
    c3 = [181/255 222/255 43/255];
    %c4 = [128/255 128/255 128/255];

    x_y = 1:numel(mean(rsa(1:21,:),1));
    x2_y = [x_y, fliplr(x_y)];
    x_o = 1:numel(mean(rsa(22:42,:),1));
    x2_o = [x_o, fliplr(x_o)];


    subplot(2,3,iroi)
    young = plot(mean(matrix_allsub(1:21,:,iroi)),'Color',c2, 'LineWidth', 1.6)
    hold on
    old = plot(mean(matrix_allsub(22:42,:,iroi)),'Color',c3, 'LineWidth', 1.6)
    hold on
    plot(significant_time_points_y, y_significants_y,'.','Color',c2)
    hold on
    plot(significant_time_points_o, y_significants_o,'.','Color',c3)
    hold on
    plot([CI_low_y(1,2);CI_high_y(1,2)],[0.29;0.29],'Color',c2, 'LineWidth', 1.6) 
    hold on
    plot(peak_latency_y(1,2),0.29,'.','MarkerSize',10,'Color',c2)
    hold on
    plot([CI_low_o(1,2);CI_high_o(1,2)],[0.28;0.28],'Color',c3, 'LineWidth', 1.6) 
    hold on
    plot(peak_latency_o(1,2),0.28,'.','MarkerSize',10,'Color',c3)
    hold on
    upper = mean(matrix_allsub(1:21,:,iroi)) + SEM_y;
    lower = mean(matrix_allsub(1:21,:,iroi)) - SEM_y;
    inBetween = [upper, fliplr(lower)];
    fill(x2_y, inBetween, c2, 'FaceAlpha', 0.155, 'LineStyle', 'none');
    hold on;
    upper = mean(matrix_allsub(22:42,:,iroi)) + SEM_o;
    lower = mean(matrix_allsub(22:42,:,iroi)) - SEM_o;
    inBetween = [upper, fliplr(lower)];
    fill(x2_o, inBetween, c3, 'FaceAlpha', 0.155, 'LineStyle', 'none');
    hold on;
    xlabel('time (ms)')
    ylabel('Spearman R')
    xticks([20 40 60 80 100 120 140 160 180 200 220 240])
    set(gca, 'XTickLabel', [-100 0 100 200 300 400 500 600 700 800 900 1000])
    yline(0,'LineStyle','--', 'LineWidth', 1.5);
    xline(40, 'LineStyle','--', 'LineWidth', 1.5);
    ylim([-0.1,0.3])
    title(sprintf("%s", roi))
    set(gca,'box','off')
    legend([young old],'younger', 'older')
    legend('boxoff')
end 

%% save
%keyboard %manually adjust the plot size
saveas(gca,sprintf('%sfusion.png',output_path));
end 
