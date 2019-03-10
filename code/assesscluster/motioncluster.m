% scrub high motion frames

a=clock;
rng(a(6));

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
savedir = fullfile(masterdir,'analyses','centroids','motion_scrub');
mkdir(savedir);
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));

% load motion files

rTR = 120; nTR = 225;
motion_concTS = zeros(nobs*(rTR+nTR),1);
for N = 1:nobs
	% declare filenames for motion files
	restfname = ['data/motion_files/',num2str(demoLTN.bblid(N)),'_',num2str(demoLTN.scanid(N)),'_relRMS.1D'];
	% store rest motion in vector corresponding to each TR in concTS
	motion_concTS((1+rTR*(N-1)):(rTR*N)) = dlmread(restfname);
end

motion_thresh = 0.1;		% remove frames with motion > 2mm
disp(['Removing ',num2str(sum(motion_concTS > motion_thresh)),' TRs with > ',num2str(motion_thresh),' mm framewise displacement']);
concTS = concTS(motion_concTS < motion_thresh,:);
nreps = 20;		% 20 reps is plenty to hit the low error solutions

[partition,centroids] = kmeans(concTS,numClusters,'Distance',distanceMethod,'Replicates',nreps);
centroids = centroids';	% transpose centroids for consistency with plotting script
save(fullfile(savedir,['kMeansMotionScrub',num2str(motion_thresh),'mm_k',num2str(numClusters),'.mat']),'partition','centroids');