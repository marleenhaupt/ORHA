function [] = statistics_time_time_MH(decoding, stats)

if ismac
    addpath('/Users/marleen/Documents/scratch/aging/eeg/code/stats');
    BASE = '/Users/marleen/Documents/scratch/aging/eeg/derivatives/grouplevel_timegen/';
elseif isunix
    addpath('/scratch/haum92/aging/eeg/code/stats');
    BASE = '/scratch/haum92/aging/eeg/derivatives/grouplevel_timegen_new/';
end

if decoding == 1
    decoding = 'animinanim';
elseif decoding == 2
    decoding = 'category';
elseif decoding == 3
    decoding = 'image'; 
end

load(sprintf('%s%s/decodingAcc_all.mat', BASE,decoding));
data= eval(sprintf('decodingAcc_all')); 

n_perm = 10000; %10000;
q_value = 0.05;

%average over dimensions 1 and 2 
%squeeze(mean(testmatrix,[1 2],'omitnan')) is the same as doing
%squeeze(mean(testmatrix,1,'omitnan')) two times after each other

if strcmp(stats, 'perm')
    out_path_results = sprintf('%s%s/statistic/fdr_perm/',BASE, decoding);

    if ~isfolder(out_path_results)
    mkdir(out_path_results);
    end
    
    [SignificantVariables_all,~,adjusted_pvalues_all] = fdr_corrected_perm_test_MH(data, n_perm, q_value,'right');
    save(sprintf('%ssignificant_variables_time_time_all.mat',out_path_results),'SignificantVariables_all');
    save(sprintf('%sadjusted_pvalues_time_time_all.mat',out_path_results),'adjusted_pvalues_all');
    [SignificantVariables_y,~,adjusted_pvalues_y] = fdr_corrected_perm_test_MH(data(1:21,:,:), n_perm, q_value,'right');
    save(sprintf('%ssignificant_variables_time_time_y.mat',out_path_results),'SignificantVariables_y');
    save(sprintf('%sadjusted_pvalues_time_time_y.mat',out_path_results),'adjusted_pvalues_y');
    [SignificantVariables_o,~,adjusted_pvalues_o] = fdr_corrected_perm_test_MH(data(22:43,:,:), n_perm, q_value,'right');
    save(sprintf('%ssignificant_variables_time_time_o.mat',out_path_results),'SignificantVariables_o');
    save(sprintf('%sadjusted_pvalues_time_time_o.mat',out_path_results),'adjusted_pvalues_o');
    [SignificantVariables_agediff,~,adjusted_pvalues_agediff] = fdr_corrected_perm_test_2samples_MH(data(1:21,:),data(22:43,:), n_perm, q_value,'both');
    save(sprintf('%ssignificant_variables_time_time_agediff.mat',out_path_results),'SignificantVariables_agediff');
    save(sprintf('%sadjusted_pvalues_time_time_agediff.mat',out_path_results),'adjusted_pvalues_agediff');
    
elseif strcmp(stats,'cluster')
    out_path_results = sprintf('%s%s/statistic/cluster_based_perm/',BASE, decoding);

    if ~isfolder(out_path_results)
    mkdir(out_path_results);
    end
    
    cluster_thr = 0.05;
    significance_thr = 0.05;
    
    [SignificantVariables_all,significantVarMax_all,pValWei_all,pValMax_all,clusters_all] = permutation_cluster_1sample_weight_alld(data, n_perm, cluster_thr, significance_thr,'right');
    save(sprintf('%ssignificant_variables_time_time_all.mat',out_path_results),'SignificantVariables_all','significantVarMax_all','pValWei_all','pValMax_all','clusters_all');
    [SignificantVariables_y,significantVarMax_y,pValWei_y,pValMax_y,clusters_y] = permutation_cluster_1sample_weight_alld(data(1:21,:,:), n_perm, cluster_thr, significance_thr,'right');
    save(sprintf('%ssignificant_variables_time_time_y.mat',out_path_results),'SignificantVariables_y','significantVarMax_y','pValWei_y','pValMax_y','clusters_y');
    [SignificantVariables_o,significantVarMax_o,pValWei_o,pValMax_o,clusters_o] = permutation_cluster_1sample_weight_alld(data(22:43,:,:), n_perm, cluster_thr, significance_thr,'right');
    save(sprintf('%ssignificant_variables_time_time_o.mat',out_path_results),'SignificantVariables_o','significantVarMax_o','pValWei_o','pValMax_o','clusters_o');  
    [SignificantVariables_agediff,significantVarMax_agediff,pValWei_agediff,pValMax_agediff,clusters_agediff] = permutation_cluster_2samples_weight_alld_MH(data(1:21,:,:),data(22:43,:,:), n_perm, cluster_thr, significance_thr,'both');
    save(sprintf('%ssignificant_variables_time_time_agediff.mat',out_path_results),'SignificantVariables_agediff','significantVarMax_agediff','pValWei_agediff','pValMax_agediff','clusters_agediff');

end
end

