function [] = eeg_time_decoding(decodinglevel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EEG time decoding script for:
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
ft_defaults;
%addpath('../toolboxes/libsvm');
addpath(genpath('/Users/marleen/Documents/MATLAB/libsvm-3.24/'));

nsubs=43;

%% Time decoding
%for i_sub = 1:2
for i_sub = 1:nsubs
    
    if strcmp(decodinglevel, 'image')
        [result]=image_decoding(i_sub);
    	matrix_allsubs(i_sub,:,:) = result;
    	save('../output/EEG_RDM_allsubs.mat','matrix_allsubs');
    elseif strcmp(decodinglevel, 'category')
        [result]=category_decoding(i_sub);
        matrix_allsubs(i_sub,:,:) = result;
        %matrix_avg_allsub(:,:) = squeeze(nanmean(matrix_allsubs,2));
    elseif strcmp(decodinglevel, 'animacy')
        [result]=animacy_decoding(i_sub);
        matrix_allsubs(i_sub,:,:) = result;
    end
    
end

matrix_avg_allsub(:,:) = squeeze(nanmean(matrix_allsubs,2));
% save accuracy per time point for all subjects
save(sprintf('../output/EEG_%s_acc_allsub.mat',decodinglevel),'matrix_avg_allsub');

%% Statistics and Plotting
% resulting plot will be saved in output folder

% load accuracy per time point for all subjects if you do not want to run
% time decoding above
%load(sprintf('../output/EEG_%s_acc_allsub.mat',decodinglevel));

plot_statistics_decoding(matrix_avg_allsub,decodinglevel);

end