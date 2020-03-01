addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
savedir = fullfile(masterdir,'clusterAssignments');

nreps = 10;
load(['data/iprtimeseries',name_root,'.mat'],'iprTS');
[partition_ipr] = kmeans(iprTS,numClusters,'Distance',distanceMethod,'Replicates',nreps);
save(fullfile(savedir,['IPRpartition_k',num2str(numClusters)]),'partition_ipr')