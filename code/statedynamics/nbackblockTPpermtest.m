addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));
load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

scanlab = {'RestComb','nBackComb'}; 

savedir = [masterdir,'/analyses/nbackblocks'];

%% load block-specific transition probabilities

ZeroBack = load(fullfile(savedir,['TransProbs0back_k',num2str(numClusters),name_root,'.mat']),'BlockTransitionProbability','BlockTransitionProbabilityMats');
OneBack = load(fullfile(savedir,['TransProbs1back_k',num2str(numClusters),name_root,'.mat']),'BlockTransitionProbability','BlockTransitionProbabilityMats');
TwoBack = load(fullfile(savedir,['TransProbs2back_k',num2str(numClusters),name_root,'.mat']),'BlockTransitionProbability','BlockTransitionProbabilityMats');

% calculate group average trans probs for each block, remove persistence
% use nanmean because if a state was not found in a block for a subject you don't observe transitions
% but you don't know if true transition probability is 0 or just unable to estimate how low it was

ZeroBackTP = squeeze(nanmean(ZeroBack.BlockTransitionProbabilityMats,1));
OneBackTP = squeeze(nanmean(OneBack.BlockTransitionProbabilityMats,1));
TwoBackTP = squeeze(nanmean(TwoBack.BlockTransitionProbabilityMats,1));

%% non-parametric permutation testing to compare 0-back, 1-back, 2-back transition probabilities

nperms = 100000;
pvals_2b0b = PERM_TEST(TwoBack.BlockTransitionProbabilityMats,ZeroBack.BlockTransitionProbabilityMats,nperms);

%% plot
f=figure;
subplot(1,3,2);
TwoBackMinusZeroBackTP = (TwoBackTP-ZeroBackTP);
imagesc(TwoBackMinusZeroBackTP.*~eye(numClusters)); colormap('plasma');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('State at t'); xlabel('State at t + 1');
sig_thresh = 0.05 / (numClusters^2);      % bonferroni correction, for two-tailed p-values so only
[y,x] = find(pvals_2b0b.*~eye(numClusters) < sig_thresh);
text(x-.12,y+.12,'*','Color','w');
caxis_bound = max(max(abs(TwoBackMinusZeroBackTP.*~eye(numClusters))));
h = colorbar; ylabel(h,'2back - 0back'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('2-back > 0-back');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];

saveas(f,fullfile(savedir,['TwoBackvsZeroBack_nonpar_k',num2str(numClusters),'.pdf']));

%% Persistence probabilities: rest vs n-back

[y,x] = find(diag(pvals_2b0b)' < sig_thresh);
f=figure;
imagesc(diag(TwoBackMinusZeroBackTP)');
caxis_bound = max(max(abs(diag(TwoBackMinusZeroBackTP))));
xticks([]); yticks([]);
text(x-.12,y+.12,'*','Color','w');
caxis([-caxis_bound caxis_bound]); colormap('plasma'); %colorbar
%h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
disp(['Persistence probability n-back - rest colorbar bounds: +/-',num2str(caxis_bound)])
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];
saveas(f,fullfile(savedir,['PersistTwoBackvsZeroBack_nonpar_k',num2str(numClusters),'.pdf']));

% for source data file
save(fullfile(savedir,['FigS8d__nBackMinusRest_k',num2str(numClusters),'.mat']),'TwoBackMinusZeroBackTP','pvals_2b0b');
