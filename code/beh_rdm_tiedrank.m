function [] = beh_rdm_tiedrank(agegroup)

output_path = '../output/';

%% load average RDM

if strcmp(agegroup, 'y')
    behfile = fullfile('../data/beh/beh_RDM/beh_rdm_young_avg.mat');
    load(behfile);  
    beh_rdm=beh_rdm_young_avg;
elseif strcmp(agegroup, 'o')
    behfile = fullfile('../data/beh/beh_RDM/beh_rdm_old_avg.mat');
    load(behfile);  
    beh_rdm=beh_rdm_old_avg;
end

%% convert to tied ranks
ranks= tiedrank(beh_rdm);

%% fill RDM
ranked_matrix= zeros(64, 64);
ranked_matrix(tril(true(64),-1)) = ranks;
ranked_matrix = ranked_matrix + tril(ranked_matrix, -1)';


%% plot RDM
figure;
imagesc(ranked_matrix);
colormap(viridis);
colorbar;
set(gca, 'XTick', [16, 32, 48, 64]);
set(gca, 'YTick', [16, 32, 48, 64]);
axis square;

%% save
saveas(gca,sprintf('%sBehavior_RDM_tiedrank_%s.png',output_path,agegroup));

