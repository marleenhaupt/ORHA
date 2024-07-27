function [rdm_rsa] = rsa_roi_beh_behavg(i_sub)

subj = sprintf('sub%02d',i_sub);
    
%load avg beh rdm
if i_sub < 22 % younger participants
	behfile = fullfile('../data/beh/beh_RDM/beh_rdm_young_avg.mat');
	load(behfile);  
	beh_rdm=beh_rdm_young_avg';
elseif i_sub > 21 % older participants
	behfile = fullfile('../data/beh/beh_RDM/beh_rdm_old_avg.mat');
	load(behfile);  
	beh_rdm=beh_rdm_old_avg';
end

%load fmri
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

    %perform RSA 
    rdm_rsa(i) = corr(beh_rdm,fmri_vector,'type','Spearman');
end 

end