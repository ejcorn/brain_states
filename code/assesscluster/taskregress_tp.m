addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
savedir = fullfile(masterdir,'analyses','centroids','task_regress');
mkdir(savedir);
load(fullfile(savedir,['kMeansnBackTaskRegressed_k',num2str(numClusters),'.mat']),'partition','centroids');
task_regressed_centroids = centroids;

% reorder task-regressed data centroids by similarity to untouched data centroids
load([masterdir,'/clusterAssignments/k',num2str(numClusters),name_root,'.mat'])
PNCcentroids = clusterAssignments.(['k',num2str(numClusters)]).bestCentroid;
PNCnames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;
PNCvsScrub = corr(task_regressed_centroids,PNCcentroids);

[~,shuffleIdx] = max(PNCvsScrub,[],1);
if length(unique(shuffleIdx)) == length(shuffleIdx)
	task_regressed_centroidsPNCorder = task_regressed_centroids(:,shuffleIdx);
else
	task_regressed_centroidsPNCorder = task_regressed_centroids;
end
clusterNames = NAME_CLUSTERS_ANGLE(task_regressed_centroidsPNCorder);
[clusterNamesUp,clusterNamesDown] = NAME_CLUSTERS_UP_DOWN(task_regressed_centroidsPNCorder);

idx = false(length(partition),numClusters);
% task regress cluster K needs to become original cluster shuffleIdx(K)

partitionPNCorder = zeros(length(partition),1);
for K = 1:numClusters 
	idx(:,K) = partition == K;  % find each task regress cluster cluster
end

for K = 1:numClusters
	partitionPNCorder(idx(:,K)) = shuffleIdx(K);  % replace with cluster # most similar to original clusters
end

[taskregressrTP,taskregressresttransprobs] = GET_TRANS_PROBS(partitionPNCorder(scanInd == 0),subjInd(scanInd == 0));
[taskregressnTP,taskregressnBacktransprobs] = GET_TRANS_PROBS(partitionPNCorder(scanInd == 1),subjInd(scanInd == 1));
taskregressrTP = mean(taskregressrTP,1);
taskregressnTP = mean(taskregressnTP,1);

% compare HCP and PNC 
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = [];

load([masterdir,'/analyses/transitionprobabilities/RestCombTransitionProbabilities_k',num2str(numClusters),name_root,'.mat'])
pncrTP = mean(transitionProbability,1);
load([masterdir,'/analyses/transitionprobabilities/nBackCombTransitionProbabilities_k',num2str(numClusters),name_root,'.mat'])
pncnTP = mean(transitionProbability,1);

% correlate TPs and PPs between HCP and PNC
restTPsim = round(corr(pncrTP(offDiag)',taskregressrTP(offDiag)'),2,'significant');
restPPsim = round(corr(pncrTP(onDiag)',taskregressrTP(onDiag)'),2,'significant');
nbackTPsim = round(corr(pncnTP(offDiag)',taskregressnTP(offDiag)'),2,'significant');
nbackPPsim = round(corr(pncnTP(onDiag)',taskregressnTP(onDiag)'),2,'significant');


% plot
taskregressgrpAvgRest = squeeze(mean(taskregressresttransprobs,1)) .* ~eye(numClusters);
taskregressgrpAvgnBack = squeeze(mean(taskregressnBacktransprobs,1)) .* ~eye(numClusters);
maxVal = max(max([taskregressgrpAvgRest,taskregressgrpAvgnBack]));

f = figure;

subplot(1,3,1);
imagesc(taskregressgrpAvgRest);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title({'Rest',['r = ',num2str(restTPsim),' with original']});
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,3,2);
imagesc(taskregressgrpAvgnBack);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);

ylabel('State at t'); xlabel('State at t + 1');
title({'n-back',['r = ',num2str(nbackTPsim),' with original']});
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

f.PaperUnits = 'inches';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];

saveas(f,fullfile(savedir,['TaskRegressedTransProbs_k',num2str(numClusters),name_root,'.pdf']),'pdf');
