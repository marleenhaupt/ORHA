function [] = plot_statistics_rsa_fmri_beh(matrix_allsub)

n_perm=100; %originally 10000
significance_thr = 0.05;
output_path = '../output/';


%% Statistics
%statistics code subtracts chance level; as these are corr values, we
%have to add 50 here
rsa = matrix_allsub+50;
rsa(40,:)=[];

[SignificantVariables_agediff,~,adjusted_pvalues_agediff] = fdr_corrected_perm_test_2sample(rsa(1:21,:),rsa(22:42,:),n_perm, significance_thr,'right');
[SignificantVariables_young,~,adjusted_pvalues_young] = fdr_corrected_perm_test_1sample(rsa(1:21,:),n_perm, significance_thr,'right');
[SignificantVariables_old,~,adjusted_pvalues_old] = fdr_corrected_perm_test_1sample(rsa(22:42,:),n_perm, significance_thr,'right');

% prerequisites for plotting
%c1 = [0/255 0/255 0/255];
c2= [62/255 73/255 137/255];
c3 = [181/255 222/255 43/255];
%c4 = [128/255 128/255 128/255];

%% Simple plot
roi_labels = {'V1','V2','V3','FFA','PPA','LOC'};
figure 
       
for iroi = 1:length(roi_labels)  
    roi = string(roi_labels(1,iroi));
    subplot(2,3,iroi)
    values = [mean(matrix_allsub(1:21,iroi));mean(matrix_allsub(22:42,iroi))];
    h = bar(values,'FaceColor','flat','EdgeColor','flat');
    h.CData(1,:) = c2; 
    h.CData(2,:) = c3;
    ylim([-0.05,0.1])
    xticklabels({'younger','older'})
    %xlabel('Group');
    ylabel('Spearman correlation coefficient');
    title(sprintf("%s", roi))
    set(gca,'box','off')
end 

%% save
saveas(gca,sprintf('%sRSA_fMRI_beh.png',output_path));
