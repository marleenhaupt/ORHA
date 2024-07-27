function [] = fmri_psa(roi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fMRI pattern similarity analysis (PSA) script for:
%% Healthy aging delays and dedifferentiates high-level visual representations

% Based on the input ROI, it calculates the pattern similarity of trials
% from the preferred category of this ROI (within, excluding correlation
% between identical trials) and of trials belonging to another category
% (between). In the end we get one file for all subjects per category
% comparison and per ROI, e.g. PPA with places vs. faces.

% Marleen Haupt
% 31.07.2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Required input
% roi: 'PPA','FFA','LOC'

%% Required data
% fMRI betas averaged over all 20 runs have to be downloaded
% and located in: ../data/fMRI/fMRI_betas

% Make sure your current directory is /code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../code/'));
%addpath('../toolboxes/spm');

output_path = '../output/';

n_subs=43;
allsub_avg_within=[];
allsub_avg_between_one=[];
allsub_avg_between_two=[];
allsub_avg_between_three=[];

if strcmpi(roi, 'PPA') 
    prefid =32; %places
    cat1 = 'faces';
    cat2 = 'objects';
elseif strcmpi(roi, 'FFA')
    prefid =0; %faces
    cat1 = 'places';
    cat2 = 'objects';
elseif strcmpi(roi, 'LOC')
    prefid=48; %objects
    cat1 = 'faces';
    cat2 = 'places';
end
cat3 = 'animals';

for i_sub = 1:n_subs
    subj = sprintf('sub%02d',i_sub);
    roimask = fullfile(sprintf('../data/fMRI/fMRI_preprocessed/%s/roi/%s_Kanwisher_nmost300.nii',subj,roi));
    
    if i_sub==40; fprintf('Excluding %i because of beh performance\n',i_sub);  continue;end;
    
    for i = 1:16
        
        idx=0+prefid+i; %place trials
        
        spmfile = fullfile(sprintf('../data/fMRI/fMRI_preprocessed/%s/PSA/Image_%d_beta_avg.nii,1',subj,idx));

        beta1=extract_roi_data_vector(roimask, spmfile)';
        
        for j = 1:16
        
            %within => we are only calculating within for pref category of
            %the roi
            jdw=0+prefid+j; %97 is last beta of run1
        
            spmfile_within = fullfile(sprintf('../data/fMRI/fMRI_preprocessed/%s/PSA/Image_%d_beta_avg.nii,1',subj,jdw));

            betaw=extract_roi_data_vector(roimask, spmfile_within)';
            
            if i==j
                within(i,j)=nan;
            else
               within(i,j)=corr(beta1(:,1),betaw(:,1),'rows','complete');
            end
            
            %between cat 1
            if strcmpi(cat1, 'faces')
                jdb=0+0+j; %97 is last beta of run1    
            elseif strcmpi(cat1, 'places')
                jdb=0+32+j; 
            elseif strcmpi(cat1, 'objects')
                jdb=0+48+j; 
            end
            spmfile_between = fullfile(sprintf('../data/fMRI/fMRI_preprocessed/%s/PSA/Image_%d_beta_avg.nii,1',subj,jdb));

            betab_one=extract_roi_data_vector(roimask, spmfile_between)';
            
            between_one(i,j)=corr(beta1(:,1),betab_one(:,1),'rows','complete');
            
            %between cat 2
            if strcmpi(cat2, 'faces')
                jdbb=0+0+j; %97 is last beta of run1    
            elseif strcmpi(cat2, 'places')
                jdbb=0+32+j; 
            elseif strcmpi(cat2, 'objects')
                jdbb=0+48+j; 
            end
            spmfile_between = fullfile(sprintf('../data/fMRI/fMRI_preprocessed/%s/PSA/Image_%d_beta_avg.nii,1',subj,jdbb));

            betab_two=extract_roi_data_vector(roimask, spmfile_between)';
            
            between_two(i,j)=corr(beta1(:,1),betab_two(:,1),'rows','complete');
            
            %between cat 3
            jdbbb=0+16+j;
            spmfile_between = fullfile(sprintf('../data/fMRI/fMRI_preprocessed/%s/PSA/Image_%d_beta_avg.nii,1',subj,jdbbb));

            betab_three=extract_roi_data_vector(roimask, spmfile_between)';
            
            between_three(i,j)=corr(beta1(:,1),betab_three(:,1),'rows','complete');

            
        end
    end
    %Z transformation
    within_Z=atanh(within);
    between_one_Z=atanh(between_one);
    between_two_Z=atanh(between_two); 
    between_three_Z=atanh(between_three); 
    %average
    avg_within=mean(within_Z,'all','omitnan');
    avg_between_one=mean(between_one_Z,'all','omitnan');
    avg_between_two=mean(between_two_Z,'all','omitnan');
    avg_between_three=mean(between_three_Z,'all','omitnan');
    allsub_avg_within = [allsub_avg_within; avg_within];
    allsub_avg_between_one = [allsub_avg_between_one; avg_between_one];
    allsub_avg_between_two = [allsub_avg_between_two; avg_between_two];
    allsub_avg_between_three = [allsub_avg_between_three; avg_between_three];
    allsub_avg_diff_one = allsub_avg_within - allsub_avg_between_one;
    allsub_avg_diff_two = allsub_avg_within - allsub_avg_between_two;
    allsub_avg_diff_three = allsub_avg_within - allsub_avg_between_three;
end

mean_y_one = mean(allsub_avg_diff_one(1:21));
mean_o_one = mean(allsub_avg_diff_one(22:42));
mean_y_two = mean(allsub_avg_diff_two(1:21));
mean_o_two = mean(allsub_avg_diff_two(22:42));
mean_y_three = mean(allsub_avg_diff_three(1:21));
mean_o_three = mean(allsub_avg_diff_three(22:42));

std_y_one = std(allsub_avg_diff_one(1:21));
std_o_one = std(allsub_avg_diff_one(22:42));
std_y_two = std(allsub_avg_diff_two(1:21));
std_o_two = std(allsub_avg_diff_two(22:42));
std_y_three = std(allsub_avg_diff_three(1:21));
std_o_three = std(allsub_avg_diff_three(22:42));

%[SignificantVariables_one,~,adjusted_pvalues_one] = fdr_corrected_perm_test_2sample(allsub_avg_diff_one(1:21,:),allsub_avg_diff_one(22:42,:),10000, 0.05,'right');
%[SignificantVariables_two,~,adjusted_pvalues_two] = fdr_corrected_perm_test_2sample(allsub_avg_diff_two(1:21,:),allsub_avg_diff_two(22:42,:),10000, 0.05,'right');
%[SignificantVariables_three,~,adjusted_pvalues_three] = fdr_corrected_perm_test_2sample(allsub_avg_diff_three(1:21,:),allsub_avg_diff_three(22:42,:),10000, 0.05,'right');

%% Simple plot

means = [
    mean(allsub_avg_diff_one(1:21)), mean(allsub_avg_diff_one(22:42));
    mean(allsub_avg_diff_two(1:21)), mean(allsub_avg_diff_two(22:42));
    mean(allsub_avg_diff_three(1:21)), mean(allsub_avg_diff_three(22:42))
];

std_devs = [
    std(allsub_avg_diff_one(1:21)), std(allsub_avg_diff_one(22:42));
    std(allsub_avg_diff_two(1:21)), std(allsub_avg_diff_two(22:42));
    std(allsub_avg_diff_three(1:21)), std(allsub_avg_diff_three(22:42))
];


x_labels = {'Category1', 'Category2', 'Category3'};
y_color = [62/255 73/255 137/255]; 
o_color = [181/255 222/255 43/255]; 


figure;
b = bar(means, 'grouped');
b(1).FaceColor = y_color;
b(2).FaceColor = o_color;
hold on;
numGroups = size(means, 1);
numBars = size(means, 2);
groupwidth = min(0.8, numBars/(numBars + 1.5));
for i = 1:numBars
    % X position for each bar
    x = (1:numGroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numBars);
    
    % Plot error bars
    errorbar(x, means(:, i), std_devs(:, i), 'k', 'linestyle', 'none', 'LineWidth', 1.5);
end
set(gca, 'XTickLabel', x_labels);
xlabel('Condition');
ylabel('Within-Between');
legend({'younger', 'older'}, 'Location', 'NorthEast');
hold off;
saveas(gca,sprintf('%sfMRI_PSA_%s.png',output_path, roi));

end