addpaths;

load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));

if scan == 'N'
    scanlab = {'nBack'}; numTRs = 225;
    scanInd = scanInd * 0;  %because scanInd is all 1's 
elseif scan == 'C'
    scanlab = {'nBackComb'}; numTRs = 225; % index rest num TRs (120) and nback to loop through
end

savedir = [masterdir,'/analyses/nbackblocks/'];
mkdir(savedir);
load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
kClusterAssignments = clusterAssignments.(['k',num2str(numClusters)]).partition;
kClusterAssignments = kClusterAssignments(scanInd == 1);	% only look at n-back state labels
subjInd = subjInd(scanInd==1);
kClusterCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

% load n-back blocks
load(fullfile(datadir,'nbackBlocks.mat'));

numBlocks = length(unique(nbackBlocks));
BlockLabels = [0 1 2];
BlockNames = {'0back','1back','2back'};

cd(savedir);
for B = 1:numBlocks	% iterate through blocks and calculate dwell time for each subject in each cluster
	BlockDwellTime = zeros(nobs,numClusters);
	BlockDuration = sum(nbackBlocks == BlockLabels(B));
	for N = 1:nobs
		subjMask = subjInd == N;
		tmpAssignments = kClusterAssignments(subjMask);	% select state labels for one subject at a time
		for K = 1:numClusters			
				% calculate time spent in each state as percentage of time spent in each block
			BlockDwellTime(N,K) = sum(tmpAssignments(nbackBlocks == BlockLabels(B)) == K) / BlockDuration;
		end
	end
	save(['DwellTime',BlockNames{B},'_k',num2str(numClusters),name_root,'.mat'],'BlockDwellTime');
end

%% calculate dwell times for n-back excluding rest blocks

BlockDwellTime = zeros(nobs,numClusters);
BlockDuration = sum(nbackBlocks ~= BlockLabels(end));
for N = 1:nobs	
	subjMask = subjInd == N;
	tmpAssignments = kClusterAssignments(subjMask);	% select state labels for one subject at a time
	for K = 1:numClusters			
			% calculate time spent in each state as percentage of time spent in non-rest blocks
		BlockDwellTime(N,K) = sum(tmpAssignments(nbackBlocks ~= BlockLabels(end)) == K) / BlockDuration;
	end
	save(['nBackRestExcludeDwellTime_k',num2str(numClusters),name_root,'.mat'],'BlockDwellTime');
end


% *** doesn't address case where blocks are not contiguous in time ***
for B = 1:numBlocks
	% will go: 0-back, 1-back, 2-back, rest
	BlockMask = repmat(nbackBlocks == BlockLabels(B),[nobs 1]);	%select only time points in each block
	BlockSubjInd = subjInd(BlockMask);
	tmpAssignments = kClusterAssignments(BlockMask);
	[transitionProbability, transitionProbabilityMats] = GET_TRANS_PROBS(tmpAssignments,BlockSubjInd);
	save(['TransProbs',BlockNames{B},'_k',num2str(numClusters),name_root,'.mat'],'transitionProbability','transitionProbabilityMats');

end

