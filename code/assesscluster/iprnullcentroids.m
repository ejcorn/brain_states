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
clusterNamesUp = cellstr(strings(numClusters,1));
clusterNamesDown = cellstr(strings(numClusters,1));


save(fullfile(savedir,['CentroidsIPRNull_k',num2str(numClusters),'.mat']),'kClusterCentroids','clusterNames','clusterNamesUp','clusterNamesDown');
