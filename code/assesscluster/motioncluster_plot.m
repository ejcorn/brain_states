% makes Fig S7

addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
savedir = fullfile(masterdir,'analyses','centroids','motion_scrub');
mkdir(savedir);
motion_thresh = 0.1;
load(fullfile(savedir,['kMeansMotionScrub',num2str(motion_thresh),'mm_k',num2str(numClusters),'.mat']),'partition','centroids');
scrubbed_centroids = centroids;

% reorder motion-scrubbed centroids by similarity to non-motion-scrubbed centroids
load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
PNCcentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
PNCnames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
PNCvsScrub = corr(scrubbed_centroids,PNCcentroids);
[~,shuffleIdx] = max(PNCvsScrub,[],1);

if length(unique(shuffleIdx)) == length(shuffleIdx)
	scrubbed_centroidsPNCorder = scrubbed_centroids(:,shuffleIdx);
else
	scrubbed_centroidsPNCorder = scrubbed_centroids;
end
clusterNames = NAME_CLUSTERS_ANGLE(scrubbed_centroidsPNCorder);
[clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(scrubbed_centroidsPNCorder);

f = figure;
PNCvsScrubSpatialCorr = PNCvsScrub(shuffleIdx,:);
imagesc(PNCvsScrubSpatialCorr); colormap('plasma');
ylabel('Motion Scrubbed'); xlabel('Full Sample'); axis square
yticks(1:numClusters); xticks(1:numClusters);
yticklabels(clusterNames); xticklabels(PNCnames);
xtickangle(90);
COLOR_TICK_LABELS(true,true,numClusters);
h=colorbar; caxis([-1 1]); h.Ticks = [-1 0 1]; h.TickLabels = [-1 0 1];
title({'Spatial Correlation',['All F.D. < ',num2str(motion_thresh),' mm']});
set(gca,'FontSize',8);
f.PaperUnits = 'inches';
f.PaperSize = [2.7 2.7];
f.PaperPosition = [0 0 2.7 2.7];
saveas(f,fullfile(savedir,['MotionScrub',num2str(motion_thresh),'mm_SpatialCorr_k',num2str(numClusters),name_root,'.pdf']),'pdf');

kClusterCentroids = scrubbed_centroidsPNCorder;	% rename centroid variable to use with one plotting script
save(fullfile(savedir,['CentroidsMotionScrubbed',num2str(motion_thresh),'mm_k',num2str(numClusters),'.mat']),'kClusterCentroids','scrubbed_centroids','clusterNames','clusterNamesUp','clusterNamesDown','PNCvsScrubSpatialCorr');

% make scatter plots of task-regressed centroids vs. PNC centroids

purple = [0.2941 0 0.5098];	% rgb values
msz = 3; % marker size
f = figure;
for K = 1:numClusters
	subplot(1,numClusters,K);
	scatter(zscore(PNCcentroids(:,K)),zscore(scrubbed_centroidsPNCorder(:,K)),msz,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0);
	prettifyEJC;
	ylabel([clusterNames{K}]); xlabel([PNCnames{K}]); axis square
	title(['r = ',num2str(round(corr(PNCcentroids(:,K),scrubbed_centroidsPNCorder(:,K)),2,'significant'))]);
end
f.PaperUnits = 'centimeters';
f.PaperSize = [21 6];
f.PaperPosition = [0 0 21 6];
saveas(f,fullfile(savedir,['MotionScrubbedVsOriginal_Scatter_k',num2str(numClusters),name_root,'.pdf']),'pdf');