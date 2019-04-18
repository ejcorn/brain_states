addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
load(['data/iprtimeseries',name_root,'.mat'],'iprTS');
savedir = fullfile(masterdir,'analyses','centroids','ipr_null');
mkdir(savedir);

load(fullfile(masterdir,'clusterAssignments',['IPRpartition_k',num2str(numClusters)]),'partition_ipr')

kClusterCentroids = zeros(nparc,numClusters);
for K = 1:numClusters
	kClusterCentroids(:,K) = mean(iprTS(partition_ipr == K,:),1);
end

clusterNames = cellstr(num2str([1:numClusters]'));
clusterNamesUp = repmat({' '},numClusters,1);
clusterNamesDown = repmat({' '},numClusters,1);

save(fullfile(savedir,['CentroidsIPRNull_k',num2str(numClusters),'.mat']),'kClusterCentroids','clusterNames','clusterNamesUp','clusterNamesDown');

load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
PNCcentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
PNCnames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
PNCvsIPR = corr(kClusterCentroids,PNCcentroids);
[~,shuffleIdx] = max(PNCvsIPR,[],1);

if length(unique(shuffleIdx)) == length(shuffleIdx)
	IPR_centroidsPNCorder = kClusterCentroids(:,shuffleIdx);
else
	IPR_centroidsPNCorder = kClusterCentroids;
end

IPRNullSpatialCorr = PNCvsIPR(shuffleIdx,:);
save(fullfile(savedir,['FigS4f__IPRNullSpatialCorr_k',num2str(numClusters),name_root,'.mat']),'IPRNullSpatialCorr');

f = figure;
imagesc(IPRNullSpatialCorr); colormap('plasma');
ylabel('IPR Null'); xlabel('Full Sample'); axis square
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
saveas(f,fullfile(savedir,['IPRNullSpatialCorr_k',num2str(numClusters),name_root,'.pdf']),'pdf');