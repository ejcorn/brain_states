addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));

rTR = 120; nTR = 225;
downSampleTask = false(nTR,1); downSampleTask(1:rTR) = true; % construct mask to downsample task so that there are equal numbers of task and rest TRs
downSampleTask = [true(rTR*nobs,1); repmat(downSampleTask,[nobs 1])];
nreps = 20;

% cluster only on first 120 TRs of task
[partition,~,sumd] = kmeans(concTS(downSampleTask,:),numClusters,'Distance',distanceMethod,'Replicates',nreps);

savedir = fullfile(masterdir,'analyses','choosing_k','rvn');
mkdir(savedir);
save(fullfile(savedir,['DownSampleTaskCluster_k',num2str(numClusters),name_root,'.mat']),'partition','sumd');

