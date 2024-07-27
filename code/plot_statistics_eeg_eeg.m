function [] = plot_statistics_eeg_eeg(matrix_allsub)

n_perm=100; %originally 10000
significance_thr = 0.05;
cluster_thr = 0.05;
output_path = '../output/';


%% Statistics entire matrix
%statistics code subtracts chance level; as these are corr values, we
%have to add 50 here
data = matrix_allsub+50;

[SignificantVariables_all,significantVarMax_all,pValWei_all,pValMax_all,clusters_all] = permutation_cluster_1sample(data, n_perm, cluster_thr, significance_thr,'right');


%% Plotting entire matrix
rsa_time_avg=squeeze(nanmean(matrix_allsub,1));
xticks_vals=[0 40 80 120 160 200 240];
yticks_vals=[0 40 80 120 160 200 240];

figure
imagesc(rsa_time_avg') 
hold on
[B,~] = bwboundaries(SignificantVariables_all');
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end
set(gca,'YDir','normal');
title('RSA EEG young - EEG old')
xlabel('time younger (ms)')
ylabel('time older (ms)')
xticks(xticks_vals)
set(gca, 'XTickLabel', [-200 0 200 400 600 800 1000])
yticks(yticks_vals)
set(gca, 'YTickLabel', [-200 0 200 400 600 800 1000])
colormap(viridis)
colorbar
a = colorbar
yline(40,'color', '#FFFFFF' ,'LineStyle','--', 'LineWidth', 1.5);
xline(40,'color', '#FFFFFF' ,'LineStyle','--', 'LineWidth', 1.5);
hline = refline([1 0]);
hline.Color = '#FFFFFF';
hline.LineStyle = '--';
hline.LineWidth = 1.5;
set(gca,'box','off')
axis square

%% Save entire matrix
saveas(gca,sprintf('%sRSA_EEG_EEG.png',output_path));

%% Statistics difference

for i = 1:240
    for j = 1:240
        for isub = 1:22
        % Subtract the element below the diagonal from the element above the diagonal
        rsa_time_diff_persub(isub,i,j) = matrix_allsub(isub,i,j) - matrix_allsub(isub,j,i);
        end
    end
end

data = rsa_time_diff_persub+50;
[SignificantVariables_all,significantVarMax_all,pValWei_all,pValMax_all,clusters_all] = permutation_cluster_1sample(data, n_perm, cluster_thr, significance_thr,'right');

rsa_time_diff=squeeze(nanmean(rsa_time_diff_persub,1));
mask = triu(true(size(rsa_time_diff)), 1);  % Upper triangular matrix without the diagonal
rsa_time_diff(~mask) = NaN;% Set elements below the diagonal to NaN


%% Plotting difference
xticks_vals=[0 40 80 120 160 200];
yticks_vals=[0 40 80 120 160 200];
zeropoint =40;

figure
imagesc(rsa_time_diff') 
hold on
[B,~] = bwboundaries(SignificantVariables_all');
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end
set(gca,'YDir','normal');
xlabel('time young(ms)')
ylabel('time old(ms)')
xticks(xticks_vals)
set(gca, 'XTickLabel', [-200 0 200 400 600 800 1000])
yticks(yticks_vals)
set(gca, 'YTickLabel', [-200 0 200 400 600 800 1000])
colormap(viridis)
colorbar
a = colorbar
yline(zeropoint,'color', '#FFFFFF' ,'LineStyle','--', 'LineWidth', 1.5);
xline(zeropoint,'color', '#FFFFFF' ,'LineStyle','--', 'LineWidth', 1.5);
hline = refline([1 0]);
hline.Color = '#FFFFFF';
hline.LineStyle = '--';
hline.LineWidth = 1.5;
set(gca,'box','off')
axis square

%% Saving difference plot
saveas(gca,sprintf('%sRSA_EEG_EEG_diff.png',output_path));
end