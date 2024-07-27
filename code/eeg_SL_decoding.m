function [] = eeg_SL_decoding(decodinglevel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EEG searchlight decoding script for:
%% Healthy aging delays and dedifferentiates high-level visual representations

% Marleen Haupt
% 31.07.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Required input
% decodinglevel: 'image', 'category', 'animacy'
% agegroup: 'y','o'

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

for i_sub = 1:nsubs
    
    if strcmp(decodinglevel, 'image')
        [result]=image_sl_peak(i_sub);
    	matrix_avg_allsubs(i_sub,:) = squeeze(result);
    elseif strcmp(decodinglevel, 'category')
        [result]=category_sl_peak(i_sub);
        matrix_avg_allsubs(i_sub,:) = squeeze(result);
    elseif strcmp(decodinglevel, 'animacy')
        [result]=animacy_sl_peak(i_sub);
        matrix_avg_allsubs(i_sub,:) = squeeze(result);
    end
    
end

%save
save(sprintf('../output/EEG_SL_%s_acc_allsub.mat',decodinglevel),'matrix_avg_allsubs');

%% Statistics and Plotting
% resulting plot will be saved in output folder

% load result for all subjects if you do not want to run analysis above
load(sprintf('../output/EEG_SL_%s_acc_allsub.mat',decodinglevel));

plot_statistics_eeg_sl(matrix_avg_allsubs,decodinglevel,'y');
plot_statistics_eeg_sl(matrix_avg_allsubs,decodinglevel,'o');

end