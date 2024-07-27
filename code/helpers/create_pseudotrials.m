function [pseudo_trials] = create_pseudotrials(num_conditions, num_trials_per_bin,n_pseudotrials, data)
%UNTITLED4 Summary of this function goes here
    

    %preallocate for pseudotrials 
    pseudo_trials=NaN(num_conditions,n_pseudotrials,size(data,3),size(data,4)); %pre-allocate memory for matrix saving pseudotrials
    % average min_num_of_trials into pseudotrials 
    
    for pseudoX=1:n_pseudotrials-1 %average trials such that we get 6 bins 
        trial_selector=(1+(pseudoX-1)*num_trials_per_bin):(num_trials_per_bin+(pseudoX-1)*num_trials_per_bin); %select trials to be averaged
         %pseudo_trialD % permutedD averaged into pseudo trials, i.e. of dimensions M * L * T
         pseudo_trials(:,pseudoX,:,:)= nanmean(data(:,trial_selector,:,:),2); %assign pseudo trial to pseudo_trial_D
    end
       % pseudo_trials(:,pseudoX,:,:)= mean(data_both_categories_standard(:,trial_selector,:,:),2); %assign pseudo trial to pseudo_trial_D

       pseudo_trials(:,n_pseudotrials,:,:) = nanmean(data(:,(1+(n_pseudotrials-1)*num_trials_per_bin):end,:,:),2);
end

