addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
savedir = fullfile(masterdir,'analyses','centroids','task_regress');
mkdir(savedir);
load(fullfile(savedir,['kMeansnBackTaskRegressed_k',num2str(numClusters),'.mat']),'partition','centroids');
task_regressed_centroids = centroids;

% reorder task-regressed data centroids by similarity to untouched data centroids
load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
PNCcentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
PNCnames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
PNCvsScrub = corr(task_regressed_centroids,PNCcentroids);
[~,shuffleIdx] = max(PNCvsScrub,[],1);
if length(unique(shuffleIdx)) == length(shuffleIdx)
	task_regressed_centroidsPNCorder = task_regressed_centroids(:,shuffleIdx);
else
	task_regressed_centroidsPNCorder = task_regressed_centroids;
end
clusterNames = NAME_CLUSTERS_ANGLE(task_regressed_centroidsPNCorder);
[clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(task_regressed_centroidsPNCorder);

f = figure;
imagesc(PNCvsScrub(shuffleIdx,:)); colormap('plasma');
ylabel('Task-Regressed'); xlabel('Full Sample'); axis square
yticks(1:numClusters); xticks(1:numClusters);
yticklabels(clusterNames); xticklabels(PNCnames);
xtickangle(90);
COLOR_TICK_LABELS(true,true,numClusters);
h=colorbar; caxis([-1 1]); h.Ticks = [-1 0 1]; h.TickLabels = [-1 0 1];
title({'Spatial Correlation'});
set(gca,'FontSize',8);
f.PaperUnits = 'inches';
f.PaperSize = [2.7 2.7];
f.PaperPosition = [0 0 2.7 2.7];
saveas(f,fullfile(savedir,['TaskRegress_SpatialCorr_k',num2str(numClusters),name_root,'.pdf']),'pdf');

kClusterCentroids = task_regressed_centroidsPNCorder;	% rename centroid variable to use with one plotting script
save(fullfile(savedir,['CentroidsnBackTaskRegressed_k',num2str(numClusters),'.mat']),'kClusterCentroids','task_regressed_centroids','clusterNames','clusterNamesUp','clusterNamesDown');

% make scatter plots of task-regressed centroids vs. PNC centroids

purple = [0.2941 0 0.5098];	% rgb values
msz = 3; % marker size
f = figure;
for K = 1:numClusters
	subplot(1,numClusters,K);
	scatter(PNCcentroids(:,K),task_regressed_centroidsPNCorder(:,K),msz,'MarkerFaceColor',purple,'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0);
	prettifyEJC;
	ylabel([clusterNames{K}]); xlabel([PNCnames{K}]); axis square
	title(['r = ',num2str(round(corr(PNCcentroids(:,K),task_regressed_centroidsPNCorder(:,K)),2,'significant'))]);
end
f.PaperUnits = 'centimeters';
f.PaperSize = [22 6];
f.PaperPosition = [0 0 22 6];
saveas(f,fullfile(savedir,['TaskRegressVsOriginal_Scatter_k',num2str(numClusters),name_root,'.pdf']),'pdf');