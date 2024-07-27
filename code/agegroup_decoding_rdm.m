function [decodingAccuracy_avgperm] = agegroup_decoding_rdm()

n_perm=5; %originally 1000

%load rdms
for i_sub = 1:43
    subj = sprintf('sub%02d',i_sub);
    eegfile = fullfile(sprintf('../data/EEG/EEG_RDM/%s/EEG_RDM.mat',subj));
    load(eegfile);
    rdm(i_sub,:,:)=decodingAccuracy_matrixup;

end

% Load data and labels
data = rdm;
labels = [ones(1,21) 2*ones(1,22)];
time_points = length(data(1,1,:));

% Split subjects by age group
age_group_1_indices = find(labels == 1);
age_group_2_indices = find(labels == 2);

for test_subject_index = 1:43
    % Ensure the test subject is excluded from training selection
    remaining_indices = setdiff(1:43, test_subject_index);
    
    % Split remaining subjects by age group
    remaining_age_group_1_indices = intersect(remaining_indices, age_group_1_indices);
    remaining_age_group_2_indices = intersect(remaining_indices, age_group_2_indices);
    
    for perm = 1:n_perm

        rng('shuffle');
    
        % Randomly select 20 subjects from each age group for training
        train_indices_1 = randsample(remaining_age_group_1_indices, 20);
        train_indices_2 = randsample(remaining_age_group_2_indices, 20);
    
        % Combine training indices
        train_indices = [train_indices_1 train_indices_2];

        for time = 1:time_points

            % Prepare training and testing data
            train_data = data(train_indices, :, time);
            train_labels = labels(train_indices);
    
            test_data = squeeze(data(test_subject_index, :,time)); 
            test_labels = labels(test_subject_index);
    
            %disp('Train the SVM');
            train_param_str=  '-s 0 -t 0 -b 0 -c 1 -q';
            model= svmtrain(train_labels',train_data,train_param_str); 
        
            %disp('Test the SVM');
            [~, accuracy, ~] = svmpredict(test_labels',test_data,model,'-q');  
            decodingAccuracy(perm, time, test_subject_index)=accuracy(1);
        end
    end
end

%average over folds
decodingAccuracy_avgperm = squeeze(nanmean(decodingAccuracy,1)); 

end

