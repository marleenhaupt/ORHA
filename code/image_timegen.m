function [decodingAccuracy_avg] = image_timegen(i_sub)
% define variables
n_permutations = 5; %originally 100 
n_pseudotrials = 8;
n_conditions = 64; 

% read in preprocessed EEG data
subj = sprintf('sub%02d',i_sub);
filepath_preprocessed_data = sprintf('../data/EEG/EEG_preprocessed/%s/stepxx_timelock.mat',subj);
load(filepath_preprocessed_data)
data = data_rej_channel_interpolated_timelocked;

%define required information 
time_points = size(data.time,2);

%downsampling from 200Hz to 50Hz for time and memory efficiency
cfg=[];
cfg.resamplefs=50;
data = ft_resampledata(cfg,data);

%% define required information

time_points = size(data.time,2);

% get minimum number of trials 
[min_number_of_trials_all, individual_objects] = get_min_trial_per_object(data);
min_number_of_trials = min(min_number_of_trials_all);

% Preallocate
decodingAccuracy_objects_time_time=NaN(n_permutations,n_conditions,n_conditions,time_points,time_points);


for perm = 1:n_permutations
    %% MVNN 
    data_matrix_MVNN = create_data_matrix_MVNN_MH(n_conditions, min_number_of_trials, data, 'image',individual_objects); 
    % get inverted covariance matrix
    [data_catA_catB_MVNN, ~] = multivariate_noise_normalization(data_matrix_MVNN);
    num_trials_per_bin = round(min_number_of_trials/n_pseudotrials);
    pseudo_trials = create_pseudotrials(n_conditions, num_trials_per_bin, n_pseudotrials, data_catA_catB_MVNN);
    
    for objA = 1:n_conditions
        for objB = 1:n_conditions
            for time1 = 1:time_points
                training_data =[squeeze(pseudo_trials(objA,1:end-1,:,time1)) ; squeeze(pseudo_trials(objB,1:end-1,:,time1))];
                labels_train  = [ones(1,n_pseudotrials-1) 2*ones(1,n_pseudotrials-1)];
                labels_test   = [1 2];
                %disp('Train the SVM');
                train_param_str=  '-s 0 -t 0 -b 0 -c 1 -q';
                model=svmtrain(labels_train',training_data,train_param_str);
                for time2 = 1:time_points
                    %disp('Test the SVM');
                    testing_data  =[squeeze(pseudo_trials(objA,end,:,time2))' ; squeeze(pseudo_trials(objB,end,:,time2))'];
                    [~, accuracy, ~] = svmpredict(labels_test',testing_data,model, '-q');
                    decodingAccuracy_objects_time_time(perm,objA,objB,time1,time2)=accuracy(1);
                end
            end
        end
    end
end    

%average over permutations
decodingAccuracy_matrix = squeeze(mean(decodingAccuracy_objects_time_time,1)); 
indicator = triu(true(size(decodingAccuracy_matrix(:,:,1,1))),1);
for time1 = 1:time_points
    for time2 = 1:time_points
        helper = decodingAccuracy_matrix(:,:,time1,time2);
        decodingAccuracy_avg(:,time1,time2) = mean(helper(indicator),'all');
    end
end

  
end