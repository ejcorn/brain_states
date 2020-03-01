% makes Figure S9
addpaths;
load(fullfile(basedir,['data/Demographics',name_root,'.mat']));
load(fullfile(datadir,['TimeSeriesIndicators',name_root,'.mat']));

scanlab = {'RestComb','nBackComb'}; 

savedir = [masterdir,'/analyses/transitionprobabilities'];

rng('shuffle');

%% load transition probabilities
disp('loading transition probabilities')

load([masterdir,'/analyses/transitionprobabilities/',scanlab{1},'TransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']);
restTransitionProbabilityMats = transitionProbabilityMats; clear transitionProbabilityMats
load([masterdir,'/analyses/transitionprobabilities/',scanlab{2},'TransitionProbabilityMatrices_k',num2str(numClusters),name_root,'.mat']);
nBackTransitionProbabilityMats = transitionProbabilityMats; clear transitionProbabilityMats

load(fullfile(masterdir,['clusterAssignments/k',num2str(numClusters),name_root,'.mat']));
clusterNames = clusterAssignments.(['k',num2str(numClusters)]).clusterNames;

%% permute rest and n-back within subjects

nperms = 100000;
pvals_twotail = PERM_TEST(nBackTransitionProbabilityMats,restTransitionProbabilityMats,nperms);

%% Plot 
disp('plot')

cd(savedir);

grpAvgRest = squeeze(mean(restTransitionProbabilityMats,1));
grpAvgRest = grpAvgRest .* ~eye(numClusters);
grpAvgnBack = squeeze(mean(nBackTransitionProbabilityMats,1)); %restGrpAvg(logical(eye(numClusters))) = 0;
grpAvgnBack = grpAvgnBack .* ~eye(numClusters);
maxval = max(max([grpAvgRest,grpAvgnBack]))

% add asterisks for significance relative to null
lt = 0.025/(2*numClusters^2); ut = 1-lt;   % two-tailed Bonferroni corrected thresholds over rest and n-back
load(['randnull/',scanlab{1},'Null_ShuffledStatesD_k',num2str(numClusters),'.mat'],'probExceedNull')
[ygt,xgt] = find((probExceedNull > ut).*~eye(numClusters));
[ylt,xlt] = find((probExceedNull < lt).*~eye(numClusters));

f=figure;
subplot(1,3,1);
imagesc(grpAvgRest); 
text(xgt,ygt,'+','Color','w');
text(xlt,ylt,'-','Color','w');
caxis([0 maxval]); colorbar
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title('Rest');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

clear probExceedNull
load(['randnull/',scanlab{2},'Null_ShuffledStatesD_k',num2str(numClusters),'.mat'],'probExceedNull')
[ygt,xgt] = find((probExceedNull > ut).*~eye(numClusters));
[ylt,xlt] = find((probExceedNull < lt).*~eye(numClusters));
cd(savedir);
subplot(1,3,2); 
imagesc(grpAvgnBack); 
text(xgt,ygt,'+','Color','w');
text(xlt,ylt,'-','Color','w');
caxis([0 maxval]); colorbar
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma');
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title('n-back');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

subplot(1,3,3);
nBackMinusRestTPMat = (grpAvgnBack-grpAvgRest);
imagesc(nBackMinusRestTPMat.*~eye(numClusters)); colormap('plasma');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('State at t'); xlabel('State at t + 1');
sig_thresh = 0.05 / numClusters^2;      % bonferroni correction, for two-tailed p-values so only
[y,x] = find(pvals_twotail.*~eye(numClusters) < sig_thresh);
text(x-.12,y+.12,'*','Color','w');
%h=colorbar; ylabel(h,'-log_{10}(p)'); caxis([0 log10(nperms)]); h.Ticks = [0 ut lt log10(nperms)]; h.TickLabels = [0 round(ut,2,'significant') round(lt,2,'significant') log10(nperms)];
caxis_bound = max(max(abs(nBackMinusRestTPMat.*~eye(numClusters))));
h = colorbar; ylabel(h,'nback - rest'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('n-back > Rest');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];

saveas(f,['RestvsnBack_nonpar_k',num2str(numClusters),'.pdf']);

%% Persistence probabilities: rest vs n-back
% retrieve diagonal again
grpAvgRest = squeeze(mean(restTransitionProbabilityMats,1));    
grpAvgnBack = squeeze(mean(nBackTransitionProbabilityMats,1));
nBackMinusRestTPMat = (grpAvgnBack-grpAvgRest);

[y,x] = find(diag(pvals_twotail)' < sig_thresh);
f=figure;
imagesc(diag(nBackMinusRestTPMat)');
caxis_bound = max(max(abs(diag(nBackMinusRestTPMat))));
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
saveas(f,['PersistRestvsnBack_nonpar_k',num2str(numClusters),'.pdf']);
