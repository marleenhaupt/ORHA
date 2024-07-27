function [min_number_of_trials, individual_objects] = get_min_trial_per_object(data)

    [individual_objects, idx] = unique(data.trialinfo(:,1));
    [~, idx_category_sorted] = sort(data.trialinfo(idx,3));
    individual_objects = individual_objects(idx_category_sorted);
    
    min_number_of_trials = NaN(size(individual_objects,1),1);
    
    
    for idx = 1:size(individual_objects,1)
        min_number_of_trials(idx) = sum(data.trialinfo(:,1)==individual_objects(idx),'all');
    end

end

