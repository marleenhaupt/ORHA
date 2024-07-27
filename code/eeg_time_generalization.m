function [] = eeg_time_generalization(decodinglevel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EEG time generalization script for:
%% Healthy aging delays and dedifferentiates high-level visual representations

% Marleen Haupt
% 31.07.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Required input
% decodinglevel: 'image', 'category', 'animacy'

%% Required data
% Preprocessed, downsampled (240 time points) EEG data has to be downloaded
% and located in: ../data/EEG/EEG_preprocessed

% Make sure your current directory is /code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path
addpath(genpath('../code/'));
%addpath('../toolboxes/fieldtrip');
%ft_defaults;
%addpath('../toolboxes/libsvm');
addpath(genpath('/Users/marleen/Documents/MATLAB/libsvm-3.24/'));

nsubs=43;

%% Time decoding
%for i_sub = 1:2
for i_sub = 1:nsubs
    
    if strcmp(decodinglevel, 'image')
        % further downsampled for time and memory efficiency
        [result]=image_timegen(i_sub);
    elseif strcmp(decodinglevel, 'category')
        [result]=category_timegen(i_sub);
    elseif strcmp(decodinglevel, 'animacy')
        [result]=animacy_timegen(i_sub);
    end
    
    matrix_avg_allsub(i_sub,:,:) = squeeze(result);
end

% save accuracy per time point for all subjects
save(sprintf('../output/EEG_%s_timegen_allsub.mat',decodinglevel),'matrix_avg_allsub');

%% Statistics and Plotting
% resulting plot will be saved in output folder

% load accuracy per time point for all subjects if you do not want to run
% time generalization above
load(sprintf('../output/EEG_%s_timegen_allsub.mat',decodinglevel));

plot_statistics_timegen(matrix_avg_allsub,decodinglevel);

end