addpaths;
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','hcpLR');
nparc = 462;
rTR = round(405); nTR = 405;

load(fullfile(savedir,['HCP_XHcentroids_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.mat']));
load(fullfile(savedir,['HCP_XHpartition_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.mat']));

[hcprTP,HCPresttransprobs] = GET_TRANS_PROBS(HCPpartitionPNCorder(HCPscanInd == 0),HCPsubjInd(HCPscanInd == 0));
[hcpnTP,HCPnBacktransprobs] = GET_TRANS_PROBS(HCPpartitionPNCorder(HCPscanInd == 1),HCPsubjInd(HCPscanInd == 1));
hcprTP = mean(hcprTP,1);
hcpnTP = mean(hcpnTP,1);

% compare HCP and PNC 
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = [];

load([masterdir,'/analyses/transitionprobabilities/RestCombTransitionProbabilities_k',num2str(numClusters),name_root,'.mat'])
pncrTP = mean(transitionProbability,1);
load([masterdir,'/analyses/transitionprobabilities/nBackCombTransitionProbabilities_k',num2str(numClusters),name_root,'.mat'])
pncnTP = mean(transitionProbability,1);

% correlate TPs and PPs between HCP and PNC
restTPsim = round(corr(pncrTP(offDiag)',hcprTP(offDiag)'),2,'significant');
restPPsim = round(corr(pncrTP(onDiag)',hcprTP(onDiag)'),2,'significant'); disp(['Rest PP corr: ',num2str(restPPsim)]);
nbackTPsim = round(corr(pncnTP(offDiag)',hcpnTP(offDiag)'),2,'significant');
nbackPPsim = round(corr(pncnTP(onDiag)',hcpnTP(onDiag)'),2,'significant'); disp(['n-back PP corr: ',num2str(nbackPPsim)]);

% fold difference in TP and PP
mean(pncrTP(offDiag) ./ hcprTP(offDiag))
mean(hcprTP(onDiag) ./ pncrTP(onDiag))
mean(pncnTP(offDiag) ./ hcpnTP(offDiag))
mean(hcpnTP(onDiag) ./ pncnTP(onDiag))

%% test rest vs. n-back transitions

restTransitionProbabilityMats = HCPresttransprobs;
nBackTransitionProbabilityMats = HCPnBacktransprobs;

%% permute rest and n-back within subjects
disp('start permutation testing')
nperms = 100000;
pvals_twotail = PERM_TEST(nBackTransitionProbabilityMats,restTransitionProbabilityMats,nperms);

% plot
HCPgrpAvgRest = squeeze(mean(HCPresttransprobs,1)) .* ~eye(numClusters);
HCPgrpAvgnBack = squeeze(mean(HCPnBacktransprobs,1)) .* ~eye(numClusters);
maxVal = max(max([HCPgrpAvgRest,HCPgrpAvgnBack]));

f = figure;

subplot(1,3,1);
imagesc(HCPgrpAvgRest);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('State at t'); xlabel('State at t + 1');
title({'Rest',['r = ',num2str(restTPsim),' with PNC']});
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,3,2);
imagesc(HCPgrpAvgnBack);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);

ylabel('State at t'); xlabel('State at t + 1');
title({'n-back',['r = ',num2str(nbackTPsim),' with PNC']});
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,3,3);
grpAvgRest = squeeze(mean(restTransitionProbabilityMats,1));    
grpAvgnBack = squeeze(mean(nBackTransitionProbabilityMats,1));

subplot(1,3,3);
HCPnBackMinusRestTP = (grpAvgnBack-grpAvgRest);
imagesc(HCPnBackMinusRestTP.*~eye(numClusters)); colormap('plasma');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('State at t'); xlabel('State at t + 1');
sig_thresh = 0.05 / numClusters^2;      % bonferroni correction, for two-tailed p-values so only
[y,x] = find(pvals_twotail.*~eye(numClusters) < sig_thresh);
text(x-.12,y+.12,'*','Color','w');
caxis_bound = max(max(abs(HCPnBackMinusRestTP.*~eye(numClusters))));
h = colorbar; ylabel(h,'nback - rest'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('n-back > Rest');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];

saveas(f,fullfile(savedir,['HCP_XH_RestnBackTransProbs_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.pdf']),'pdf');

grpAvgRest = squeeze(mean(restTransitionProbabilityMats,1));    
grpAvgnBack = squeeze(mean(nBackTransitionProbabilityMats,1));

[y,x] = find(diag(pvals_twotail)' < sig_thresh);
f=figure;
imagesc(diag(HCPnBackMinusRestTP)');
caxis_bound = max(max(abs(diag(HCPnBackMinusRestTP))));
xticks([]); yticks([]);
text(x-.12,y+.12,'*','Color','w');
caxis([-caxis_bound caxis_bound]); colormap('plasma'); %colorbar

set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];
saveas(f,fullfile(savedir,['HCP_XH_RestvsnBackPersistProbs_k',num2str(numClusters),'.pdf']));

%% Plot persistence probabilities
rPP = diag(squeeze(mean(restTransitionProbabilityMats,1)))';
nPP = diag(squeeze(mean(nBackTransitionProbabilityMats,1)))';
maxval = max([rPP,nPP])

f=figure;
imagesc(rPP);
xticks([]); yticks([]);
caxis([0 maxval]); colormap('plasma'); %colorbar
disp(['Persistence probability n-back - rest colorbar bounds: +/-',num2str(caxis_bound)])
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];

saveas(f,fullfile(savedir,['HCP_XH_RestPersistenceProbs_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),'.pdf']),'pdf');

f=figure;
imagesc(nPP);
xticks([]); yticks([]);
caxis([0 maxval]); colormap('plasma'); %colorbar
%h=colorbar; h.Ticks = [0 maxval]; h.TickLabels = [0, round(maxval,2,'significant')];
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
f.PaperUnits = 'inches';
f.PaperSize = [1 .2];
f.PaperPosition = [0 0 1 .2];

saveas(f,fullfile(savedir,['HCP_XH_nBackPersistenceProbs_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),'.pdf']),'pdf');

% for source data file
save(fullfile(savedir,['FigS6d-f__HCPTransProbs_k',num2str(numClusters),'.mat']),'HCPgrpAvgRest','HCPgrpAvgnBack','HCPnBackMinusRestTP','pvals_twotail','rPP','nPP');