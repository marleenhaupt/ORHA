function [] = eeg_agegroup_decoding(timeshift)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EEG-based age group decoding script for:
%% Healthy aging delays and dedifferentiates high-level visual representations

% Marleen Haupt
% 31.07.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Required input
% please indicate whether analysis should account for timeshift 
%(peak delay of 21ms in older participants or not)

%timeshift: 'yes','no'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Required data
% EEG RDMs have to be downloaded and located in: 
% ../data/EEG/EEG_RDM

% Make sure your current directory is /code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path
addpath(genpath('../code/'));
%addpath('../toolboxes/fieldtrip');
ft_defaults;
%addpath('../toolboxes/libsvm');
addpath(genpath('/Users/marleen/Documents/MATLAB/libsvm-3.24/'));

if strcmp(timeshift, 'no')
    [result] = agegroup_decoding_rdm;
    save(sprintf('../output/agegroup_decoding_notimeshift.mat'),'result');
    load('../output/agegroup_decoding_notimeshift.mat');
    plot_statistics_agegroup_decoding(result);

elseif strcmp(timeshift, 'yes')
    [result] = agegroup_decoding_rdm_timeshift;
    save(sprintf('../output/agegroup_decoding_timeshift.mat'),'result');
    %load('../output/agegroup_decoding_timeshift.mat');
    plot_statistics_agegroup_decoding_timeshift(result);
end

end