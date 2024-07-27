function [] = plot_statistics_timegen(matrix_avg_allsubs,decodinglevel)

%% run cluster based permutation test
n_perm = 100; %originally 10000
cluster_thr = 0.05;
significance_thr = 0.05;
output_path = '../output/';
       
[SignificantVariables_young,significantVarMax_young,pValWei_young,pValMax_young,clusters_young] = permutation_cluster_1sample(matrix_avg_allsubs(1:21,:,:), n_perm, cluster_thr, significance_thr,'right');
[SignificantVariables_old,significantVarMax_old,pValWei_old,pValMax_old,clusters_old] = permutation_cluster_1sample(matrix_avg_allsubs(22:43,:,:), n_perm, cluster_thr, significance_thr,'right');
[SignificantVariables_agediff,significantVarMax_agediff,pValWei_agediff,pValMax_agediff,clusters_agediff] = permutation_cluster_2sample(matrix_avg_allsubs(1:21,:,:),matrix_avg_allsubs(22:43,:,:), n_perm, cluster_thr, significance_thr,'right');

%% specify axes for plotting because data for image decoding is further downsampled
if strcmp(decodinglevel, 'image')==1  
    xticks_vals=[0 10 20 30 40 50 60];
    yticks_vals=[0 10 20 30 40 50 60];
    zeropoint=10;
else
    xticks_vals=[0 40 80 120 160 200 240];
    yticks_vals=[0 40 80 120 160 200 240];
    zeropoint=40;
end

%% plot younger adults

data = matrix_avg_allsubs(1:21,:,:)-50;
mean_data = squeeze(mean(data,1));
mean_data = (mean_data+mean_data')/2;
SignificantVariables = (SignificantVariables_young+SignificantVariables_young')/2;


figure
imagesc(mean_data)
hold on
[B,~] = bwboundaries(SignificantVariables);
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end
set(gca,'YDir','normal');
title(sprintf("Time generalization for %s in younger adults", decodinglevel))
xlabel('training time (ms)')
ylabel('testing time (ms)')
xticks(xticks_vals)
set(gca, 'XTickLabel', [-200 0 200 400 600 800 1000])
yticks(yticks_vals)
set(gca, 'YTickLabel', [-200 0 200 400 600 800 1000])
colormap(viridis)
colorbar
%a = colorbar('Limits', [-6 9]);
a = colorbar
a.Label.String = 'classification accuracy - chance level (%)'
caxis([-5 50]) 
yline(zeropoint,'color', '#FFFFFF' ,'LineStyle','--', 'LineWidth', 1.5);
xline(zeropoint,'color', '#FFFFFF' ,'LineStyle','--', 'LineWidth', 1.5);
hline = refline([1 0]);
hline.Color = '#FFFFFF';
hline.LineStyle = '--';
hline.LineWidth = 1.5;
set(gca,'box','off')
axis square
saveas(gca,sprintf('%s%s_timegen_younger.png',output_path, decodinglevel));

%% plot older adults

data = matrix_avg_allsubs(22:43,:,:)-50;
mean_data = squeeze(mean(data,1));
mean_data = (mean_data+mean_data')/2;
SignificantVariables = (SignificantVariables_old+SignificantVariables_old')/2;

figure
imagesc(mean_data)
hold on
[B,~] = bwboundaries(SignificantVariables);
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end
set(gca,'YDir','normal');
title(sprintf("Time generalization for %s in older adults", decodinglevel))
xlabel('training time (ms)')
ylabel('testing time (ms)')
xticks(xticks_vals)
set(gca, 'XTickLabel', [-200 0 200 400 600 800 1000])
yticks(yticks_vals)
set(gca, 'YTickLabel', [-200 0 200 400 600 800 1000])
colormap(viridis)
colorbar
%a = colorbar('Limits', [-6 9]);
a = colorbar
a.Label.String = 'classification accuracy - chance level (%)'
caxis([-5 50]) 
yline(zeropoint,'color', '#FFFFFF' ,'LineStyle','--', 'LineWidth', 1.5);
xline(zeropoint,'color', '#FFFFFF' ,'LineStyle','--', 'LineWidth', 1.5);
hline = refline([1 0]);
hline.Color = '#FFFFFF';
hline.LineStyle = '--';
hline.LineWidth = 1.5;
set(gca,'box','off')
axis square
saveas(gca,sprintf('%s%s_timegen_older.png',output_path, decodinglevel));

%% plot age group difference
data = mean(matrix_avg_allsubs(1:21,:,:),1)-mean(matrix_avg_allsubs(22:43,:,:),1);
mean_data = squeeze(mean(data,1));
mean_data = (mean_data+mean_data')/2;
SignificantVariables = (SignificantVariables_agediff+SignificantVariables_agediff')/2;

figure
imagesc(mean_data)
hold on
[B,~] = bwboundaries(SignificantVariables);
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end
set(gca,'YDir','normal');
title(sprintf("Time generalization for %s: difference between age groups", decodinglevel))
xlabel('training time (ms)')
ylabel('testing time (ms)')
xticks(xticks_vals)
set(gca, 'XTickLabel', [-200 0 200 400 600 800 1000])
yticks(yticks_vals)
set(gca, 'YTickLabel', [-200 0 200 400 600 800 1000])
colormap(viridis)
colorbar
%a = colorbar('Limits', [-6 9]);
a = colorbar
a.Label.String = 'classification accuracy - chance level (%)'
caxis([-15 15]) %5-40 for non diff
yline(zeropoint,'color', '#FFFFFF' ,'LineStyle','--', 'LineWidth', 1.5);
xline(zeropoint,'color', '#FFFFFF' ,'LineStyle','--', 'LineWidth', 1.5);
hline = refline([1 0]);
hline.Color = '#FFFFFF';
hline.LineStyle = '--';
hline.LineWidth = 1.5;
set(gca,'box','off')
axis square
saveas(gca,sprintf('%s%s_timegen_agediff.png',output_path, decodinglevel));

end

