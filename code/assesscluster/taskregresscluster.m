% cluster on task-regressed time series

a=clock;
rng(a(6));

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
savedir = fullfile(masterdir,'analyses','centroids','task_regress');
mkdir(savedir);
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));

rTR = 120; nTR = 225;
task_regress_ts = zeros(nTR*nobs,nparc);
for N = 1:nobs
	disp(['Subject ',num2str(N)])
	% declare filenames for motion files
	restfname = ['data/taskregressedts/',num2str(demoLTN.bblid(N)),'_',num2str(demoLTN.scanid(N)),'_Lausanne_ROIv_scale250_timeseries.csv'];
	% store rest motion in vector corresponding to each TR in concTS
	ts = dlmread(restfname);
	task_regress_ts((1+nTR*(N-1)):(nTR*N),:) = ts(7:end,1:nparc);	% truncate initial acquisition
end

% replace task data in concTS with task-regressed data
concTS((1+rTR*nobs):(rTR*nobs + nTR*nobs),:) = task_regress_ts;

nreps = 15;		% 20 reps is plenty to hit the low error solutions

disp('begin clustering')
[partition,centroids] = kmeans(concTS,numClusters,'Distance',distanceMethod,'Replicates',nreps);
centroids = centroids';	% transpose centroids for consistency with plotting script
save(fullfile(savedir,['kMeansnBackTaskRegressed_k',num2str(numClusters),'.mat']),'partition','centroids');