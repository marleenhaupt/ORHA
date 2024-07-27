function [rdm_rsa_allroi] = rsa_time_roi_persub(i_sub)

subj = sprintf('sub%02d',i_sub);
rois={'EVC','V1','V2','V3','LOC','FFA','OFA','STS','PPA','TOS','RSC'};
selected_rois= [2,3,4,6,9,5];


for i = 1:length(selected_rois)
    iroi = selected_rois(i);
    roi = rois{iroi};
        
    %load fMRI file
    fmrifile = fullfile(sprintf('../data/fMRI/fMRI_RDM/%s/res_accuracy_matrix_minus_chance.mat',subj));
    load(fmrifile);
    fmri_rdm(:,:)=squeeze(results.accuracy_matrix_minus_chance.output{iroi,1});
    indicator = triu(true(size(fmri_rdm(:,:))),1);
    fmri_vector = fmri_rdm(indicator);
    
    %load eeg file
    eegfile = fullfile(sprintf('../data/EEG/EEG_RDM/%s/EEG_RDM.mat',subj));
    load(eegfile);
    eeg_vector = decodingAccuracy_matrixup;
        
    % preallocate result matrix
    num_tp = length(eeg_vector(1,:));
    rdm_rsa = NaN(1,num_tp);

    %perfom RSA at each EEG timepoint
    for time = 1:num_tp
        tp = time;
        rdm_rsa(time) = corr(eeg_vector(:,tp),fmri_vector,'type','Spearman');
    end

rdm_rsa_allroi(:,i)=rdm_rsa;
end

end