addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

scanlab = {'RestComb','nBackComb'}; 

savedir = [masterdir,'/analyses/nbackblocks'];


ZeroBack = load(fullfile(savedir,['TransProbs0back_k',num2str(numClusters),name_root,'.mat']),'BlockTransitionProbability','BlockTransitionProbabilityMats');
OneBack = load(fullfile(savedir,['TransProbs1back_k',num2str(numClusters),name_root,'.mat']),'BlockTransitionProbability','BlockTransitionProbabilityMats');
TwoBack = load(fullfile(savedir,['TransProbs2back_k',num2str(numClusters),name_root,'.mat']),'BlockTransitionProbability','BlockTransitionProbabilityMats');

% calculate group average trans probs for each block, remove persistence
% use nanmean because if a state was not found in a block for a subject you don't observe transitions
% but you don't know if true transition probability is 0 or just unable to estimate how low it was
%
ZeroBackTP = squeeze(nanmean(ZeroBack.BlockTransitionProbabilityMats,1));
OneBackTP = squeeze(nanmean(OneBack.BlockTransitionProbabilityMats,1));
TwoBackTP = squeeze(nanmean(TwoBack.BlockTransitionProbabilityMats,1));
%}
%{
ZeroBack.BlockTransitionProbabilityMats(isnan(ZeroBack.BlockTransitionProbabilityMats)) = 0;
OneBack.BlockTransitionProbabilityMats(isnan(OneBack.BlockTransitionProbabilityMats)) = 0;
TwoBack.BlockTransitionProbabilityMats(isnan(TwoBack.BlockTransitionProbabilityMats)) = 0;

ZeroBackTP = squeeze(mean(ZeroBack.BlockTransitionProbabilityMats,1));
OneBackTP = squeeze(mean(OneBack.BlockTransitionProbabilityMats,1));
TwoBackTP = squeeze(mean(TwoBack.BlockTransitionProbabilityMats,1));
%}

maxval = max(max([ZeroBackTP .* ~eye(numClusters),OneBackTP .* ~eye(numClusters),TwoBackTP .* ~eye(numClusters)]));

%% plot transition probabilities
f=figure;
subplot(1,3,1);
imagesc(ZeroBackTP .* ~eye(numClusters)); 
caxis([0 maxval]); colorbar
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title('0-back');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

subplot(1,3,2); 
imagesc(OneBackTP .* ~eye(numClusters)); 
caxis([0 maxval]); colorbar
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title('1-back');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

subplot(1,3,3);
imagesc(TwoBackTP .* ~eye(numClusters)); 
caxis([0 maxval]); colorbar
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title('2-back');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];

saveas(f,fullfile(savedir,['nBackBlockTransProbs_k',num2str(numClusters),'.pdf']));

%% plot persistence probabilities

maxval = max(max([ZeroBackTP(~~eye(numClusters)),OneBackTP(~~eye(numClusters)),TwoBackTP(~~eye(numClusters))]));
f=figure;
imagesc(ZeroBackTP(~~eye(numClusters))');
xticks([]); yticks([]);
caxis([0 maxval]); colormap('plasma'); %colorbar
%h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];

saveas(f,fullfile(savedir,['ZeroBackPersistenceProbs_k',num2str(numClusters),'.pdf']));

f=figure;
imagesc(OneBackTP(~~eye(numClusters))');
xticks([]); yticks([]);
caxis([0 maxval]); colormap('plasma'); %colorbar
%h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];

saveas(f,fullfile(savedir,['OneBackPersistenceProbs_k',num2str(numClusters),'.pdf']));

f=figure;
imagesc(TwoBackTP(~~eye(numClusters))');
xticks([]); yticks([]);
caxis([0 maxval]); colormap('plasma'); %colorbar
%h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2]; 

saveas(f,fullfile(savedir,['TwoBackPersistenceProbs_k',num2str(numClusters),'.pdf']));

f=figure; h=colorbar; caxis([0 maxval]); h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
saveas(f,fullfile(savedir,['nBackBlcokPersistenceProbsColorbar_k',num2str(numClusters),'.pdf']))