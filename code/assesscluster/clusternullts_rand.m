addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
savedir = fullfile(masterdir,'clusterAssignments');

nreps = 10;
load(['data/randtimeseries',name_root,'.mat'],'randTS');
[partition_randn] = kmeans(randTS,numClusters,'Distance',distanceMethod,'Replicates',nreps);;
save(fullfile(savedir,['randnpartition_k',num2str(numClusters)]),'partition_randn')
