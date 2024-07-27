function [decodingAccuracy] = animacy_decoding(i_sub)

%add fieldtrip
%ft_defaults
%add libsvm
%addpath(genpath('../toolboxes/libsvm'));

% define variables
n_permutations = 5; %originally 100
n_pseudotrials = 8;
n_conditions = 2; 

% read in preprocessed EEG data
subj = sprintf('sub%02d',i_sub);
filepath_preprocessed_data = sprintf('../data/EEG/EEG_preprocessed/%s/stepxx_timelock.mat',subj);
load(filepath_preprocessed_data)
data = data_rej_channel_interpolated_timelocked;

% specify output directory
results_dir = sprintf('../data/EEG/EEG_RDM/%s',subj);
if ~isfolder(results_dir)
    mkdir(results_dir);
end

%define required information 
time_points = size(data.time,2);
cfg = [];
data = ft_selectdata(cfg, data);
    
% get minimum number of trials 
% coded in column3 with animate=1 and inanimate=2
number_of_trial_animate = sum(data.trialinfo(:,3)==1,'all');
number_of_trial_inanimate = sum(data.trialinfo(:,3)==2,'all');

min_number_of_trials = min([number_of_trial_animate,number_of_trial_inanimate]);

% Preallocate 
decodingAccuracy=NaN(n_permutations,time_points);

%% do the actual decoding 
for perm = 1:n_permutations

    %create data matrix for MVNN
    data_MVNN = create_data_matrix_MVNN_MH(n_conditions, min_number_of_trials, data, 'animinanim');
    
    % get inverted covariance matrix
    [data_MVNN, ~] = multivariate_noise_normalization(data_MVNN); 

    %% split data into animate and inanimate trials
    data_animate = data_MVNN(1,:,:,:); 
    data_inanimate = data_MVNN(2,:,:,:);
    
    %% create pseudotrials
    %Permute the trials
    data_animate_permuted = data_animate(:, randperm(size(data_animate,2)),:,:);
    data_inanimate_permuted = data_inanimate(:,randperm(size(data_inanimate,2)),:,:);
    
    %Put both categories into one matrix
    data_both_categories = NaN([size(data_animate_permuted)]);
    data_both_categories(1,:,:,:) = data_animate_permuted;
    data_both_categories(2,:,:,:) = data_inanimate_permuted;
    
    %Split into pseudotrials
    num_trials_per_bin = round(min_number_of_trials/n_pseudotrials);
    pseudo_trials = create_pseudotrials(n_conditions, num_trials_per_bin, n_pseudotrials, data_both_categories);
    
    %% do the actual decoding 
    for time = 1:time_points
        
        % split into trainung and testing 
        training_data=[squeeze(pseudo_trials(1,1:end-1,:,time)) ; squeeze(pseudo_trials(2,1:end-1,:,time))];
        testing_data=[squeeze(pseudo_trials(1,end,:,time))' ; squeeze(pseudo_trials(2,end,:,time))'];
        
        % create labels for the SVM 
        labels_train  = [ones(1,n_pseudotrials-1) 2*ones(1,n_pseudotrials-1)];
        labels_test   = [1 2];
        
        %disp('Train the SVM');
        train_param_str=  '-s 0 -t 0 -b 0 -c 1 -q';
        model= svmtrain(labels_train',training_data,train_param_str); 
        
        %disp('Test the SVM');
        [~, accuracy, ~] = svmpredict(labels_test',testing_data,model,'-q');  
        decodingAccuracy(perm,time)=accuracy(1);     
        
    end
end

% Average over permutations
decodingAccuracy = squeeze(mean(decodingAccuracy,1)); 
     
end

