function [] = rsa_eeg_beh()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EEG-behavior RSA script for:
%% Healthy aging delays and dedifferentiates high-level visual representations

% Marleen Haupt
% 31.07.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Required data
% EEG and behavior RDMs have to be downloaded and located in: 
% ../data/EEG/EEG_RDM and ../data/beh/beh_RDM

% Make sure your current directory is /code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path
addpath(genpath('../code/'));
%addpath('../toolboxes/fieldtrip');
ft_defaults;
%addpath('../toolboxes/libsvm');
addpath(genpath('/Users/marleen/Documents/MATLAB/libsvm-3.24/'));

nsubs=43;

%% Time decoding
for i_sub = 1:nsubs
	[result]=rsa_time_beh_behavg(i_sub);
	matrix_allsub(i_sub,:) = result;
end

% save accuracy per time point for all subjects
save(sprintf('../output/RSA_EEG_beh_allsub.mat'),'matrix_allsub');

%% Statistics and Plotting
% resulting plot will be saved in output folder

% load accuracy per time point for all subjects if you do not want to run
% time decoding above
%load('../output/RSA_EEG_beh_allsub.mat');

plot_statistics_rsa_eeg_beh(matrix_allsub);

end