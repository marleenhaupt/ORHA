function [decodingAccuracy_avg]=image_sl_peak(i_sub)

if i_sub < 22
    t = 73; %young image peak
elseif i_sub > 21
    t = 77; %old image peak
end

subj = sprintf('sub%02d',i_sub);
n_permutations = 5; %originally 100
n_pseudotrials = 8;
n_conditions = 64; %exemplar/image to decode
n_channels = 64; 
chanIdx = 1:n_channels; 

%read in eeg data
filepath_preprocessed_data = sprintf('../data/EEG/EEG_preprocessed/%s/stepxx_timelock.mat',subj);
load(filepath_preprocessed_data)
data = data_rej_channel_interpolated_timelocked;

time_points = size(data.time,2);
cfg = [];
data = ft_selectdata(cfg, data);
    
% get minimum number of trials 
% minimum number of trials

[min_number_of_trials_all, individual_objects] = get_min_trial_per_object(data);
min_number_of_trials = min(min_number_of_trials_all);

%% Decoding
decodingAccuracy=NaN(n_permutations, n_conditions, n_conditions, n_channels);

for perm = 1:n_permutations
    %create data matrix for MVNN
    data_matrix_MVNN = create_data_matrix_MVNN_MH(n_conditions, min_number_of_trials, data, 'image',individual_objects); 
    % get inverted covariance matrix
    [data_catA_catB_MVNN, ~] = multivariate_noise_normalization(data_matrix_MVNN);
    num_trials_per_bin = round(min_number_of_trials/n_pseudotrials);
    pseudo_trials = create_pseudotrials(n_conditions, num_trials_per_bin, n_pseudotrials, data_catA_catB_MVNN);
    
    for objA = 1:n_conditions-1
        
        for objB = objA+1:n_conditions      
            
            for iChan = 1:64
                training_data =[pseudo_trials(objA,1:end-1,iChan,t),pseudo_trials(objB,1:end-1,iChan,t)]';
                testing_data  =[pseudo_trials(objA,end,iChan,t),pseudo_trials(objB,end,iChan,t)]';
                    
                labels_train  = [ones(1,n_pseudotrials-1) 2*ones(1,n_pseudotrials-1)];
                labels_test   = [1 2];

                %disp('Train the SVM');
                train_param_str= '-s 0 -t 0 -b 0 -c 1 -q';
                model=svmtrain(labels_train',training_data,train_param_str); 

                %disp('Test the SVM');
                [~, accuracy, ~] = svmpredict(labels_test',testing_data,model,'-q');  
                decodingAccuracy(perm,objA, objB, iChan)=accuracy(1);
                
            end

        end
        
    end
    
end

%% Save the decoding accuracies
decodingAccuracy = squeeze(mean(decodingAccuracy,1)); %average over permutations
decodingAccuracy_avg = squeeze(nanmean(decodingAccuracy,1));
decodingAccuracy_avg = squeeze(nanmean(decodingAccuracy_avg,1));
end


