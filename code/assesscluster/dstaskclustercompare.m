% makes Fig S4

addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
concTS = csvread(fullfile(datadir,['ConcTSCSV_',name_root,'.csv']));

savedir = fullfile(masterdir,'analyses','choosing_k','rvn');
load(fullfile(savedir,['DownSampleTaskCluster_k',num2str(numClusters),name_root,'.mat']));

rTR = 120; nTR = 225;
downSampleTask = false(nTR,1); downSampleTask(1:rTR) = true;
downSampleTask = [true(rTR*nobs,1); repmat(downSampleTask,[nobs 1])];
dsTS = concTS(downSampleTask,:);
subjInd = [repelem(1:nobs,rTR),repelem(1:nobs,rTR)];
scanInd = false(rTR*nobs*2,1); scanInd(rTR*nobs+1:end) = true;

dsCentroids = zeros(nparc,numClusters);
for K = 1:numClusters
    dsCentroids(:,K) = mean(dsTS(partition == K,:),1);
end


load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
overallPartition = clusterAssignments.(['k',num2str(numClusters)]).partition;
overallCentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
overallNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

dsVsOverall = corr(dsCentroids,overallCentroids);
[~,shuffleIdx] = max(dsVsOverall,[],1);
dsVsOverall = dsVsOverall(shuffleIdx,:);
dsCentroids = dsCentroids(:,shuffleIdx);
[dsNames] = NAME_CLUSTERS_ANGLE(dsCentroids);
[dsNamesUp,dsNamesDown] = NAME_CLUSTERS_UP_DOWN(dsCentroids);

% rename for plotting script
kClusterCentroids = dsCentroids;
clusterNames = dsNames;
[clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(dsCentroids);

save(fullfile(savedir,['DownSampleTaskClusterCentroids_k',num2str(numClusters),'.mat']),'dsCentroids',...
	'dsNames','dsNamesUp','dsNamesDown','kClusterCentroids','clusterNames','clusterNamesUp','clusterNamesDown');

% Fig S4d
f = figure;
imagesc(dsVsOverall); colormap('plasma'); axis square
xticks(1:numClusters); xticklabels(overallNames); xtickangle(90); xlabel('Full Sample')
yticks(1:numClusters); yticklabels(dsNames); ylabel('6 min. Task');
COLOR_TICK_LABELS(true,true,numClusters);
h=colorbar; caxis([-1 1]); h.Ticks = [-1 0 1]; h.TickLabels = [-1 0 1];
title('Spatial Correlation');
set(gca,'FontSize',8);
f.PaperUnits = 'inches';
f.PaperSize = [2.7 2.7];
f.PaperPosition = [0 0 2.7 2.7];

saveas(f,fullfile(savedir,['DownSampleTaskvsOverallClusterSpatialCorr_k',num2str(numClusters),name_root,'.pdf']),'pdf');

%% reorder centroids based on similarity to original centroids

[~,shuffleIdx] = sort(shuffleIdx);	% reorder to index indiv. cluster not overall

idx = false(length(partition),numClusters);
% rest cluster K needs to become comb cluster shuffleIdx(K)

for K = 1:numClusters 
	idx(:,K) = partition == K;  % find each rest cluster
end

for K = 1:numClusters
	partition(idx(:,K)) = shuffleIdx(K);  % replace with cluster # most similar to comb clusters
end

%% dwell times
numTRs = [rTR; nTR];
scanlab = {'RestDS','nBackDS'};

for i = 1:numel(scanlab)
	stateDuration = zeros(nobs,numClusters);

	for N = 1:nobs
	    for K = 1:numClusters
	        stateDuration(N,K) = sum(and(and(partition == K,subjInd' == N),scanInd == (i-1))) / numTRs(i);
	    end
	end

	save(fullfile(savedir,[scanlab{i},'StateDurations_k',num2str(numClusters),name_root,'.mat']),'stateDuration')
end

%% transition probabilities

[~,restTransProbs] = GET_TRANS_PROBS_NO_PERSIST(partition(scanInd == 0),subjInd(scanInd == 0));

load(fullfile(datadir,'nbackBlocks.mat'));
BlockLabels = [0 1 2];	% only iterate through 0-back, 1-back, 2-back. Ignore rest and instructions
numBlocks = length(BlockLabels);
BlockNames = {'0back','1back','2back'};
for B = 1:numBlocks
	nBackTransProbs = zeros(nobs,numClusters,numClusters);	% preallocate transition probability matrices
	for N = 1:nobs
		subjPartition = partition(scanInd' == 1 & subjInd == N);
		% using only first rTR TRs of task
		[nBackTransProbs(N,:,:),~] = GET_BLOCK_TRANS_PROBS_NO_PERSIST(subjPartition,nbackBlocks(1:rTR) == BlockLabels(B),numClusters);
	end		
end

grpAvgRest = squeeze(mean(restTransProbs,1)) .* ~eye(numClusters);
grpAvgnBack = squeeze(nanmean(nBackTransProbs,1)) .* ~eye(numClusters);
maxVal = max(max([grpAvgRest,grpAvgnBack]));

% Fig S4d
f = figure;

subplot(1,2,1);
imagesc(grpAvgRest);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(dsNames); xtickangle(90); yticklabels(dsNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title('Rest');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,2,2);
imagesc(grpAvgnBack);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(dsNames); xtickangle(90); yticklabels(dsNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title('n-back');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

f.PaperUnits = 'inches';
f.PaperSize = [16/3 4];
f.PaperPosition = [0 0 16/3 4];
saveas(f,fullfile(savedir,['DownSampleTaskRest2BackTransProbs_k',num2str(numClusters),'.pdf']));

%% compare to overall trans probs

load([masterdir,'/analyses/transitionprobabilities/RestCombTransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']);
restTransitionProbabilityMats = transitionProbabilityMats; clear transitionProbabilityMats
overallRestTransProbs = squeeze(mean(restTransitionProbabilityMats,1)).* ~eye(numClusters);
load([masterdir,'/analyses/transitionprobabilities/nBackCombTransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']);
nBackTransitionProbabilityMats = transitionProbabilityMats; clear transitionProbabilityMats
overallnBackTransProbs = squeeze(mean(nBackTransitionProbabilityMats,1)).* ~eye(numClusters);

onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:25; offDiag(onDiag) = [];

grpAvgRest = reshape(grpAvgRest',numClusters^2,1);
overallRestTransProbs = reshape(overallRestTransProbs',numClusters^2,1);
corr(grpAvgRest(offDiag),overallRestTransProbs(offDiag))

grpAvgnBack = reshape(grpAvgnBack',numClusters^2,1);
overallnBackTransProbs = reshape(overallnBackTransProbs',numClusters^2,1);
corr(grpAvgnBack(offDiag),overallnBackTransProbs(offDiag))