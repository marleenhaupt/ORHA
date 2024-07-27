function [decodingAccuracy_matrix] = image_decoding(i_sub)

%add fieldtrip
%ft_defaults
%add libsvm
%addpath(genpath('../toolboxes/libsvm'));

% define variables
n_permutations = 5; %originally 100
n_pseudotrials = 8;
n_conditions = 64; 

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
[min_number_of_trials_all, individual_objects] = get_min_trial_per_object(data);
min_number_of_trials = min(min_number_of_trials_all);

% preallocate result matrix
%decodingAccuracy=NaN(n_permutations, n_conditions, n_conditions, time_points);

%% do the actual decoding 
for perm = 1:n_permutations
    
    %create data matrix for MVNN
    data_matrix_MVNN = create_data_matrix_MVNN_MH(n_conditions, min_number_of_trials, data, 'image',individual_objects); 
    % get inverted covariance matrix
    [data_catA_catB_MVNN, ~] = multivariate_noise_normalization(data_matrix_MVNN);
    num_trials_per_bin = round(min_number_of_trials/n_pseudotrials);
    pseudo_trials = create_pseudotrials(n_conditions, num_trials_per_bin, n_pseudotrials, data_catA_catB_MVNN);
    
    for objA = 1:n_conditions
        
        for objB = 1:n_conditions 
            
            for time = 1:time_points
               
                training_data =[squeeze(pseudo_trials(objA,1:end-1,:,time)) ; squeeze(pseudo_trials(objB,1:end-1,:,time))];
                testing_data  =[squeeze(pseudo_trials(objA,end,:,time))' ; squeeze(pseudo_trials(objB,end,:,time))'];
                labels_train  = [ones(1,n_pseudotrials-1) 2*ones(1,n_pseudotrials-1)];
                labels_test   = [1 2];

                %train
                train_param_str=  '-s 0 -t 0 -b 0 -c 1 -q';
                model=svmtrain(labels_train',training_data,train_param_str); 

                %test
                [~, accuracy, ~] = svmpredict(labels_test',testing_data,model,'-q');  
                decodingAccuracy_image(perm,objA, objB, time)=accuracy(1); 
            
            end
            
        end
        
    end
    
end

% Average over permutations
decodingAccuracy = squeeze(mean(decodingAccuracy_image,1)); 
    
% Extract RDM vector/matrix
indicator = triu(true(size(decodingAccuracy(:,:,1))),1);
for time = 1:time_points
    helper = decodingAccuracy(:,:,time);
    decodingAccuracy_matrix(:,time) = helper(indicator);
end

end

