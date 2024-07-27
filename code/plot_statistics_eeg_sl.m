function plot_statistics_eeg_sl(matrix_avg_allsubs,decodinglevel,agegroup)

output_path = '../output/';

% load channel locations
load('../code/helpers/acticap-64ch-standard2.mat'); 

% delete reference and ground electrode
lay.label = lay.label([1:14,16:21,23:66],:);
lay.pos = lay.pos([1:14,16:21,23:66],:);
lay.width = lay.width([1:14,16:21,23:66],:);
lay.height = lay.height([1:14,16:21,23:66],:);
layout = lay;


if strcmp(agegroup,'y')
    n1=1;
    n2=21;
elseif strcmp(agegroup,'o')
    n1=22;
    n2=43;
end
    
%% Stats
for_stats_data = matrix_avg_allsubs(n1:n2,:);
    
num_perms = 100; %originally 10000
tail = 'right';
qvalue = 0.05;
pvalue = 0.05;
[significant_channels,pvalues,crit_p, adjusted_pvalues]= permutation_cluster_1sample(for_stats_data,num_perms,qvalue,pvalue,tail);

%% Average over subjects 
avg_over_conditions_all_subjects = squeeze(nanmean(for_stats_data,1)); 
        
%% Plot searchlight results per time point
searchlight_patterns = avg_over_conditions_all_subjects-50;
        
%check whether electrode display is correct
%ft_plot_layout(layout);

load('../code/helpers/electrodelabelorder.mat');

%ft_topoplotER requires 'comp', 'timelock' or 'freq' data as input
%create fake timelock data
%based off this: https://mailman.science.ru.nl/pipermail/fieldtrip/2021-March/040728.html
data = [];
data.weights = searchlight_patterns';
data.label = electrodelabelorder.label; %electrode label order from recording
data.time = 0; % give it the actual time so it gets displayed in plot correctly or 0 as random value
data.dimord = 'chan_time'; %just give it this name

%actual topoplot
cfg = [];
cfg.xlim ='maxmin'; %time dimension
cfg.zlim = [0 40]; %acc dimension
cfg.layout = layout; %
cfg.parameter = 'weights'; 
cfg.style = 'both_imsat'; 
cfg.style = 'straight'; %plot without contour lines
cfg.colormap = 'viridis';
cfg.colorbar =  'SouthOutside';
cfg.comment = 'no';
cfg.marker = 'off';
if any(significant_channels==1)
    cfg.highlight = 'on';
    cfg.highlightsymbol = 'o';
    cfg.highlightsize = 8;
    significant_channels = logical(significant_channels);
    cfg.highlightchannel =data.label(significant_channels);
else 
    cfg.highlight = 'off';
end
figure
ft_topoplotER(cfg,data);
        
saveas(gca,sprintf('%s%s_eeg_sl_%s.png',output_path, decodinglevel, agegroup));

end