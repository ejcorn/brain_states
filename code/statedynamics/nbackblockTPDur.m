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

BlockLabels = [0 1 2];	% only iterate through 0-back, 1-back, 2-back. Ignore rest and instructions
numBlocks = length(BlockLabels);
BlockNames = {'0back','1back','2back'};

%% Compute fractional occupancy for each state in each n-back block

for B = 1:numBlocks	% iterate through blocks and calculate dwell time for each subject in each cluster
	BlockFractionalOccupancy = zeros(nobs,numClusters);
	BlockDuration = sum(nbackBlocks == BlockLabels(B));
	for N = 1:nobs
		subjMask = subjInd == N;
		tmpAssignments = kClusterAssignments(subjMask);	% select state labels for one subject at a time
		for K = 1:numClusters			
				% calculate time spent in each state as percentage of time spent in each block
			BlockFractionalOccupancy(N,K) = sum(tmpAssignments(nbackBlocks == BlockLabels(B)) == K) / BlockDuration;
		end
	end
	save(fullfile(savedir,['FractionalOccupancy',BlockNames{B},'_k',num2str(numClusters),name_root,'.mat']),'BlockFractionalOccupancy');
end

%% Compute fractional occupancy for whole n-back task, excluding rest or instruction components

BlockFractionalOccupancy = zeros(nobs,numClusters);
BlockDuration = sum(ismember(nbackBlocks,BlockLabels));		% number of TRs in 0,1,2 back blocks, i.e. not rest or instructions
for N = 1:nobs	
	subjMask = subjInd == N;
	tmpAssignments = kClusterAssignments(subjMask);	% select state labels for one subject at a time
	for K = 1:numClusters			
			% calculate time spent in each state as percentage of time spent in non-rest blocks
		BlockFractionalOccupancy(N,K) = sum(tmpAssignments(ismember(nbackBlocks,BlockLabels)) == K) / BlockDuration;
	end
	save(fullfile(savedir,['nBackRestExcludeFractionalOccupancy_k',num2str(numClusters),name_root,'.mat']),'BlockFractionalOccupancy');
end

%% Compute dwell times for each state in each n-back block

TR = 3;		%3 seconds/TR in PNC data

for B = 1:numBlocks	% iterate through blocks and calculate dwell time for each subject in each cluster
	BlockDwellTime = zeros(nobs,numClusters);
	BlockDuration = sum(nbackBlocks == BlockLabels(B));
	for N = 1:nobs
		subjMask = subjInd == N;
		tmpAssignments = kClusterAssignments(subjMask);	% select state labels for one subject at a time			
		% calculate Dwell Time as mean length of runs of each state
		mean_dt = GET_BLOCK_DWELL_TIME(tmpAssignments,nbackBlocks == BlockLabels(B),numClusters);
		BlockDwellTime(N,:) = mean_dt*TR;	% store dwell time *in seconds*		
	end
	save(fullfile(savedir,['DwellTime',BlockNames{B},'_k',num2str(numClusters),name_root,'.mat']),'BlockDwellTime');
end

%% Transition probabilities within each block

for B = 1:numBlocks
	BlockTransitionProbability = zeros(nobs,numClusters^2);
	BlockTransitionProbabilityMats = zeros(nobs,numClusters,numClusters);	% preallocate transition probability matrices
	for N = 1:nobs
		subjPartition = kClusterAssignments(subjInd == N);
		[BlockTransitionProbabilityMats(N,:,:),BlockTransitionProbability(N,:)] = GET_BLOCK_TRANS_PROBS(subjPartition,nbackBlocks == BlockLabels(B),numClusters);
	end	
	save(fullfile(savedir,['TransProbs',BlockNames{B},'_k',num2str(numClusters),name_root,'.mat']),'BlockTransitionProbability','BlockTransitionProbabilityMats');
end

%% Transition probabilities while excluding persistence

for B = 1:numBlocks
	BlockTransitionProbability = zeros(nobs,numClusters^2);
	BlockTransitionProbabilityMats = zeros(nobs,numClusters,numClusters);	% preallocate transition probability matrices
	for N = 1:nobs
		subjPartition = kClusterAssignments(subjInd == N);
		[BlockTransitionProbabilityMats(N,:,:),BlockTransitionProbability(N,:)] = GET_BLOCK_TRANS_PROBS_NO_PERSIST(subjPartition,nbackBlocks == BlockLabels(B),numClusters);
	end	
	save(fullfile(savedir,['TransProbsNoPersist',BlockNames{B},'_k',num2str(numClusters),name_root,'.mat']),'BlockTransitionProbability','BlockTransitionProbabilityMats');
end

%% Transition probabilities only at block boundaries

BetweenBlockNumTransitions = zeros(nobs,numClusters^2);
BetweenBlockTransitionProbability = zeros(nobs,numClusters^2);
for N = 1:nobs
	subjPartition = kClusterAssignments(subjInd == N);
	[~,BetweenBlockTransitionProbability(N,:),~,BetweenBlockNumTransitions(N,:)] = GET_BETWEEN_BLOCK_TRANS_PROBS(subjPartition,nbackBlocks,numClusters);
end
save(fullfile(savedir,['BetweenBlockTransProbs_k',num2str(numClusters),name_root,'.mat']),'BetweenBlockNumTransitions','BetweenBlockTransitionProbability');