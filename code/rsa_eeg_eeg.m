function [] = rsa_eeg_eeg()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EEG-EEG RSA script for:
%% Healthy aging delays and dedifferentiates high-level visual representations

% Marleen Haupt
% 31.07.2024
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
    
%% RSA
for i_sub = 22:43 %older participants
	[result]=rsa_eeg_eeg_youngavg(i_sub);
    i=i_sub-21;
	matrix_allsub(i,:,:) = result;
end

% save accuracy per time point for all subjects
save(sprintf('../output/RSA_EEG_EEG_alloldsub.mat'),'matrix_allsub');

%% Statistics and Plotting
% resulting plot will be saved in output folder

% load accuracy per time point for all subjects if you do not want to run
% time decoding above
load('../output/RSA_EEG_EEG_alloldsub.mat');

plot_statistics_eeg_eeg(matrix_allsub);
end