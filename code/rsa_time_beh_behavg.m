function [rdm_rsa] = rsa_time_beh_behavg(i_sub)

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

%load eeg file
eegfile = fullfile(sprintf('../data/EEG/EEG_RDM/%s/EEG_RDM.mat',subj));
load(eegfile);
eeg_vector = decodingAccuracy_matrixup;
        
% preallocate result matrix
num_tp = length(eeg_vector(1,:));

%perfom RSA at each EEG timepoint
    for time = 1:num_tp
        rdm_rsa(time) = corr(beh_rdm,eeg_vector(:,time),'type','Spearman');
    end
end