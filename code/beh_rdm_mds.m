function [] = beh_rdm_mds()

output_path = '../output/';

%% load average RDMs younger and older

behfile = fullfile('../data/beh/beh_RDM/beh_rdm_young_avg.mat');
load(behfile);  
rdm_young=beh_rdm_young_avg;

behfile = fullfile('../data/beh/beh_RDM/beh_rdm_old_avg.mat');
load(behfile);  
rdm_old=beh_rdm_old_avg;

%Data parameters
cond1 = 1:16;
cond2 = 17:32;
cond3 = 33:48;
cond4 = 49:64;

%Procrustes
mds_young = cmdscale(rdm_young,2);
mds_old = cmdscale(rdm_old,2);
[d, Z, transform] = procrustes(mds_young, mds_old);

%Plotting parameters
dotsize = 150;
color1 = [225/255,190/255,106/255];
color2 = [230/255,97/255,0/255];
color3 = [64/255,176/255,166/255];
color4 = [93/255,58/255,155/255];


% plot with dots 
figure
subplot(1, 2, 1);
points1 = squeeze(mds_young(cond1,:));
points2 = squeeze(mds_young(cond2,:));
points3 = squeeze(mds_young(cond3,:));
points4 = squeeze(mds_young(cond4,:));
scatter(points1(:,1),points1(:,2),dotsize,color1,'filled');
hold on;
scatter(points2(:,1),points2(:,2),dotsize,color2,'filled');
hold on;
scatter(points3(:,1),points3(:,2),dotsize,color3,'filled');
hold on;
scatter(points4(:,1),points4(:,2),dotsize,color4,'filled');
hold on;
xlim([-0.3,0.3]);
ylim([-0.3,0.3]);
hold on;
subplot(1, 2, 2);
points1 = squeeze(Z(cond1,:));
points2 = squeeze(Z(cond2,:));
points3 = squeeze(Z(cond3,:));
points4 = squeeze(Z(cond4,:));
scatter(points1(:,1),points1(:,2),dotsize,color1,'filled');
hold on;
scatter(points2(:,1),points2(:,2),dotsize,color2,'filled');
hold on;
scatter(points3(:,1),points3(:,2),dotsize,color3,'filled');
hold on;
scatter(points4(:,1),points4(:,2),dotsize,color4,'filled');
hold on;
xlim([-0.3,0.3]);
ylim([-0.3,0.3]);
hold on;
saveas(gca,sprintf('%sBehavior_RDM_MDS_dots.png',output_path));

% plot with images
figure
subplot(1, 2, 1);
for c = 1:64
	point = squeeze(mds_young(c,:));
	image_name = fullfile('../data/stimuli',sprintf('%d.jpg',c));
	M = imread(image_name);
    image_width = 0.03; % adjust this value as needed to make images smaller/larger
    image_height = 0.03; % adjust this value as needed to make images smaller/larger
    % Display the image in the plot at the specified coordinates
    imagesc('XData', [point(1) - image_width/2, point(1) + image_width/2], ...
            'YData', [point(2) + image_height/2, point(2) - image_height/2], ...
            'CData', M);
	%imagesc('XData',[point(1)-5 point(1)+5] ,'YData', [point(2)+5 point(2)-5 ],'CData', M);
	daspect([1 1 1]); 
    xlim([-0.3,0.3]);
    ylim([-0.3,0.3]);
    hold on
	drawnow
end
subplot(1, 2, 2);
for c = 1:64
	point = squeeze(Z(c,:));
	image_name = fullfile('../data/stimuli',sprintf('%d.jpg',c));
	M = imread(image_name);
    image_width = 0.03; % adjust this value as needed to make images smaller/larger
    image_height = 0.03; % adjust this value as needed to make images smaller/larger
    % Display the image in the plot at the specified coordinates
    imagesc('XData', [point(1) - image_width/2, point(1) + image_width/2], ...
            'YData', [point(2) + image_height/2, point(2) - image_height/2], ...
            'CData', M);
	daspect([1 1 1]); 
    xlim([-0.3,0.3]);
    ylim([-0.3,0.3]);
    hold on
	drawnow
end
saveas(gca,sprintf('%sBehavior_RDM_MDS_images.png',output_path));