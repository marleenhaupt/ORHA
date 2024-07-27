function [decodingAccuracy_avg] = category_timegen(i_sub)

% define variables
n_permutations = 5; %originally 100
n_pseudotrials = 8;
n_conditions = 4; 

% read in preprocessed EEG data
subj = sprintf('sub%02d',i_sub);
filepath_preprocessed_data = sprintf('../data/EEG/EEG_preprocessed/%s/stepxx_timelock.mat',subj);
load(filepath_preprocessed_data)
data = data_rej_channel_interpolated_timelocked;

%define required information 
time_points = size(data.time,2);
cfg = [];
data = ft_selectdata(cfg, data);

% get minimum number of trials 
number_of_trial_faces = sum(data.trialinfo(:,2)==1,'all');
number_of_trial_animals = sum(data.trialinfo(:,2)==2,'all');
number_of_trial_places = sum(data.trialinfo(:,2)==3,'all');
number_of_trial_objects = sum(data.trialinfo(:,2)==4,'all');

min_number_of_trials = min([number_of_trial_faces,number_of_trial_animals,number_of_trial_places,number_of_trial_objects]);

% Preallocate
decodingAccuracy_category_time_time=NaN(n_permutations,n_conditions,n_conditions,time_points,time_points);
size(decodingAccuracy_category_time_time)

for perm = 1:n_permutations
    %% MVNN 
    data_matrix_MVNN = create_data_matrix_MVNN_MH(n_conditions, min_number_of_trials, data, 'category'); 
    % get inverted covariance matrix
    [data_catA_catB_MVNN, ~] = multivariate_noise_normalization(data_matrix_MVNN);
    num_trials_per_bin = round(min_number_of_trials/n_pseudotrials);
    pseudo_trials = create_pseudotrials(n_conditions, num_trials_per_bin, n_pseudotrials, data_catA_catB_MVNN);
    
    %looping through category dimension 1
    for objA = 1:n_conditions
        
        %looping through category dimension 2
        for objB = 1:n_conditions
            
            %training on all time points
            for time1 = 1:time_points
                training_data =[squeeze(pseudo_trials(objA,1:end-1,:,time1)) ; squeeze(pseudo_trials(objB,1:end-1,:,time1))];
                labels_train  = [ones(1,n_pseudotrials-1) 2*ones(1,n_pseudotrials-1)];
                labels_test   = [1 2];
                
                %disp('Train the SVM');
                train_param_str=  '-s 0 -t 0 -b 0 -c 1 -q';
                model=svmtrain(labels_train',training_data,train_param_str);
                
                %testing on all time points
                for time2 = 1:time_points
                    %disp('Test the SVM');
                    testing_data  =[squeeze(pseudo_trials(objA,end,:,time2))' ; squeeze(pseudo_trials(objB,end,:,time2))'];
                    [~, accuracy, ~] = svmpredict(labels_test',testing_data,model,'-q');
                    decodingAccuracy_category_time_time(perm,objA, objB, time1,time2)=accuracy(1);
                end
            end
        end
    end
end

%average over permutations
decodingAccuracy_matrix = squeeze(mean(decodingAccuracy_category_time_time,1)); 
indicator = triu(true(size(decodingAccuracy_matrix(:,:,1,1))),1);
for time1 = 1:time_points
    for time2 = 1:time_points
        helper = decodingAccuracy_matrix(:,:,time1,time2);
        decodingAccuracy_avg(:,time1,time2) = mean(helper(indicator),'all');
    end
end
 
end
