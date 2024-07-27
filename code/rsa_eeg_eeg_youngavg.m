function [rsa_time] = rsa_eeg_eeg_youngavg(i_sub)

%load average RDM young
load('../data/EEG/EEG_RDM/decodingAcc_matrix_young.mat');
rdm_young=decodingAcc_young(:,:);

%preallocate matrix
rsa_time=NaN(size(rdm_young,2),size(rdm_young,2));

%% RSA
subj = sprintf('sub%02d',i_sub);

%load RDM for specific older subject
load(sprintf('../data/EEG/EEG_RDM/%s/EEG_RDM.mat',subj));
rdm_old=decodingAccuracy_matrixup(:,:);
for time_young = 1:size(rdm_young,2)
    for time_old = 1:size(rdm_old,2)
        rsa_time(time_young, time_old) = corr(rdm_young(:,time_young),rdm_old(:,time_old),'Type','Spearman');
	end
end  

end