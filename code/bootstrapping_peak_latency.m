function [outputArg1,outputArg2] = bootstrapping_peak_latency_MH(decoding,time)

%% set up prereqs
if ismac
   addpath('/Users/marleen/Documents/MATLAB/fieldtrip-20220519/')
   ft_defaults;
   BASE = '/Users/marleen/Documents/AgingProject/EEG/analysis';
elseif isunix
    addpath('/scratch/haum92/software/libsvm/matlab/');
    BASE = '/scratch/haum92/aging/eeg/derivatives/grouplevel_decoding/';
    code_dir = '/scratch/haum92/aging/eeg/code/';
    addpath(genpath('/scratch/haum92/aging/eeg/code/'));
    addpath('/scratch/haum92/aging/eeg/derivatives/grouplevel_decoding/');
end

if decoding == 1
    decoding = 'animinanim';
elseif decoding == 2
    decoding = 'category';
elseif decoding == 3
    decoding = 'image';
end

if time == 1

        load(sprintf('%s%s/decodingAcc_all.mat', BASE,decoding));
        data = decodingAcc_all;

        load(sprintf('%s%s/decodingAcc_young.mat', BASE,decoding));
        data_y = decodingAcc_young;

        load(sprintf('%s%s/decodingAcc_old.mat', BASE,decoding));
        data_o = decodingAcc_old;


    out_path_results = sprintf('%s%s/statistic/cluster_based_perm/',BASE, decoding);

% work on this part once I get to time-time
elseif time == 2
    train = 'time_time';
    if strcmp(fixcross,'difference') == 1
        load(sprintf('%s%s/statistic/cluster_based_perm/%s_time_time/time_time_diff_curve_%s.mat',BASE, decoding, method));
        data = diff_curve;
    else
        load(sprintf('%sdata/FixEyeEEG/main/results/%s_%s/%s_decodingAcc_%s_%s.mat', BASE,decoding, train, decoding,fixcross,method));
        data= eval(sprintf('decodingAcc_%s_all',fixcross));
    end
    out_path_results = sprintf('%sdata/FixEyeEEG/main/results/statistic/cluster_based_perm/%s_%s/',BASE, decoding,train);
end

bootstrap_samples = 10000; %10000;
peak_latency_samples = NaN(bootstrap_samples,2);

rng('shuffle')
for bs = 1:bootstrap_samples
    bootstrapped_data_y = datasample(data_y,length(data_y(:,1)),1);
    bootstrapped_data_o = datasample(data_o,length(data_o(:,1)),1);
    avg_data_y = squeeze(mean(bootstrapped_data_y,1));
    avg_data_o = squeeze(mean(bootstrapped_data_o,1));
    [peaks_x_y, peaks_y_y] = find(avg_data_y==max(avg_data_y,[],'all'));
    if size(peaks_x_y,2) > 1
        peaks_x_y = round(mean(peaks_x_y));
        peaks_y_y = round(mean(peaks_y_y));
    end
     [peaks_x_o, peaks_y_o] = find(avg_data_o==max(avg_data_o,[],'all'));
    if size(peaks_x_o,2) > 1
        peaks_x_o = round(mean(peaks_x_o));
        peaks_y_o = round(mean(peaks_y_o));
    end
    peak_latency_samples_y(bs,1)=peaks_x_y;
    peak_latency_samples_y(bs,2)=peaks_y_y;
    peak_latency_samples_o(bs,1)=peaks_x_o;
    peak_latency_samples_o(bs,2)=peaks_y_o;
    peak_latency_samples_diff(bs,1)=peaks_x_o -peaks_x_y;
    peak_latency_samples_diff(bs,2)=peaks_y_o -peaks_y_y;
end
peak_latency_y = round(mean(peak_latency_samples_y,1));
peak_latency_o = round(mean(peak_latency_samples_o,1));
peak_latency_diff = round(mean(peak_latency_samples_diff,1));

pvals = sum((peak_latency_samples_diff(:,2)<=0)) / (size(peak_latency_samples_diff,1));
%[p_W, h_W]  = signrank(peak_latency_samples_diff(:,2));
%[p_T, h_T]  = ttest(peak_latency_samples_diff(:,2),0, 'Alpha', 0.05);

CI_95_y = NaN(2,2);
CI_95_y(1,:) = prctile(peak_latency_samples_y,2.5);
CI_95_y(2,:) = prctile(peak_latency_samples_y,97.5);

CI_95_o = NaN(2,2);
CI_95_o(1,:) = prctile(peak_latency_samples_o,2.5);
CI_95_o(2,:) = prctile(peak_latency_samples_o,97.5);

CI_95_diff = NaN(2,2);
CI_95_diff(1,:) = prctile(peak_latency_samples_diff,2.5);
CI_95_diff(2,:) = prctile(peak_latency_samples_diff,97.5);

save(sprintf('%speak_latency_diff.mat',out_path_results),'peak_latency_diff','CI_95_diff','peak_latency_samples_diff');

end

