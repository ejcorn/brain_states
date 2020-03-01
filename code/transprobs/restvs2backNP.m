addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));

scanlab = {'RestComb','TwoBackComb'}; 

savedir = [masterdir,'/analyses/transitionprobabilities'];

rng('shuffle');

%% load transition probabilities
disp('loading transition probabilities')

Rest = load([masterdir,'/analyses/transitionprobabilities/',scanlab{1},'TransitionProbabilitiesNoPersist_k',num2str(numClusters),name_root,'.mat']);
restTransitionProbabilityMats = Rest.transitionProbabilityMats;
TwoBack = load(fullfile(masterdir,'analyses/nbackblocks',['TransProbsNoPersist2back_k',num2str(numClusters),name_root,'.mat']),'BlockTransitionProbability','BlockTransitionProbabilityMats');
TwoBackTransitionProbabilityMats = TwoBack.BlockTransitionProbabilityMats;

load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

%% permute rest and 2-back within subjects

nperms = 100000;
pvals_twotail = PERM_TEST(TwoBackTransitionProbabilityMats,restTransitionProbabilityMats,nperms);

%% Plot 
disp('plot')

cd(savedir);

grpAvgRest = squeeze(nanmean(restTransitionProbabilityMats,1));
grpAvgRest = grpAvgRest .* ~eye(numClusters);
grpAvgTwoBack = squeeze(nanmean(TwoBackTransitionProbabilityMats,1)); %restGrpAvg(logical(eye(numClusters))) = 0;
grpAvgTwoBack = grpAvgTwoBack .* ~eye(numClusters);
maxval = max(max([grpAvgRest,grpAvgTwoBack]))

f=figure;
subplot(1,3,1);
imagesc(grpAvgRest); 
caxis([0 maxval]); colorbar
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next State');
title('Rest');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

subplot(1,3,2); 
imagesc(grpAvgTwoBack); 
caxis([0 maxval]); colorbar
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next State');
title('2-back');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

subplot(1,3,3);
TwoBackMinusRestTPMat = (grpAvgTwoBack-grpAvgRest);
imagesc(TwoBackMinusRestTPMat.*~eye(numClusters)); colormap('plasma');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Current State'); xlabel('Next State');
sig_thresh = 0.05 / numClusters^2;      % bonferroni correction, for two-tailed p-values so only
[y,x] = find(pvals_twotail.*~eye(numClusters) < sig_thresh);
text(x-.12,y+.12,'*','Color','w');
%h=colorbar; ylabel(h,'-log_{10}(p)'); caxis([0 log10(nperms)]); h.Ticks = [0 ut lt log10(nperms)]; h.TickLabels = [0 round(ut,2,'significant') round(lt,2,'significant') log10(nperms)];
caxis_bound = max(max(abs(TwoBackMinusRestTPMat.*~eye(numClusters))));
h = colorbar; ylabel(h,'2back - rest'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('2-back > Rest');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];
save('Fig4a-c,e-f__SourceData.mat','grpAvgRest','grpAvgTwoBack','TwoBackMinusRestTPMat','pvals_twotail')
saveas(f,['RestvsTwoBackNoPersist_nonpar_k',num2str(numClusters),'.pdf']);

%% Persistence probabilities: rest vs n-back
%{
% retrieve diagonal again
grpAvgRest = squeeze(nanmean(restTransitionProbabilityMats,1));    
grpAvgTwoBack = squeeze(nanmean(TwoBackTransitionProbabilityMats,1));
TwoBackMinusRestTPMat = (grpAvgTwoBack-grpAvgRest);

[y,x] = find(diag(pvals_twotail)' < sig_thresh);
f=figure;
imagesc(diag(TwoBackMinusRestTPMat)');
caxis_bound = max(max(abs(diag(TwoBackMinusRestTPMat))));
xticks([]); yticks([]);
text(x-.12,y+.12,'*','Color','w');
caxis([-caxis_bound caxis_bound]); colormap('plasma'); %colorbar
%h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
disp(['Persistence probability 2-back - rest colorbar bounds: +/-',num2str(caxis_bound)])
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];
saveas(f,['PersistRestvsTwoBack_nonpar_k',num2str(numClusters),'.pdf']);

%}