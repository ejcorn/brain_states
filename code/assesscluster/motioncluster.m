% scrub high motion frames

a=clock;
rng(a(6));

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
savedir = fullfile(masterdir,'analyses','centroids','motion_scrub');
mkdir(savedir);
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));

nback_fd = dlmread('data/framewisediplacement_nback.csv');
nback_fd = sortrows(nback_fd,2);	% sort by scan id
nback_fd = nback_fd(ismember(nback_fd(:,2),demoLTN.scanid),:);	% exclude subjects
nback_fd = nback_fd(:,10:234);	% remove scan ID and bbl id, column of 0's, and first 6 volumes

% load rest motion files and concatenate into group motion time series

rTR = 120; nTR = 225;
motion_concTS = zeros(nobs*(rTR+nTR),1);
for N = 1:nobs
	% declare filenames for motion files
	restfname = ['data/motion_files/',num2str(demoLTN.bblid(N)),'_',num2str(demoLTN.scanid(N)),'_relRMS.1D'];
	% store rest motion in vector corresponding to each TR in concTS
	motion_concTS((1+rTR*(N-1)):(rTR*N)) = dlmread(restfname);
	% store n-back motion in vector corresponding to each TR in concTS
	motion_concTS((1+rTR*nobs+nTR*(N-1)):(rTR*nobs+nTR*N)) = nback_fd(N,:);
end

motion_thresh = 0.1;		% remove frames with motion > 2mm
disp(['Removing ',num2str(sum(motion_concTS > motion_thresh)),' TRs with > ',num2str(motion_thresh),' mm framewise displacement']);
concTS = concTS(motion_concTS < motion_thresh,:);
nreps = 20;		% 20 reps is plenty to hit the low error solutions

[partition,centroids] = kmeans(concTS,numClusters,'Distance',distanceMethod,'Replicates',nreps);
centroids = centroids';	% transpose centroids for consistency with plotting script
save(fullfile(savedir,['kMeansMotionScrub',num2str(motion_thresh),'mm_k',num2str(numClusters),'.mat']),'partition','centroids');