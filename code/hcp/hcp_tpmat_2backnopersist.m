addpaths;
masterdir = fullfile(basedir,'results',name_root);
savedir = fullfile(masterdir,'analyses','hcpLR');
nparc = 462;
rTR = round(405); nTR = 405;

load(fullfile(savedir,['HCP_XHcentroids_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.mat']));
load(fullfile(savedir,['HCP_XHpartition_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.mat']));

[hcprTP,HCPresttransprobs] = GET_TRANS_PROBS_NO_PERSIST(HCPpartitionPNCorder(HCPscanInd == 0),HCPsubjInd(HCPscanInd == 0));

% figure out the timing of 2-back block onsets
TR = 0.72; % in sec
BlockDurationInTRs = round(27.5/TR); % in TRs
BlockOnsetsInTRs = round([7.985 79.195 150.5 178.582]/TR); % onset of tools, faces, places, etc. in TRs
TwoBackBlock = zeros(nTR,1);
for B = 1:length(BlockOnsetsInTRs)
	TwoBackBlock(BlockOnsetsInTRs(B):BlockOnsetsInTRs(B)+BlockDurationInTRs) = 1; % set all TRs during 2-back equal to 1
end

hcpnTP = zeros(nsubjs,numClusters^2);
HCPnBacktransprobs = zeros(nsubjs,numClusters,numClusters);	% preallocate transition probability matrices
for N = 1:nsubjs
	subjPartition = HCPpartitionPNCorder(HCPsubjInd' == N & HCPscanInd == 1);
	[HCPnBacktransprobs(N,:,:),hcpnTP(N,:)] = GET_BLOCK_TRANS_PROBS_NO_PERSIST(subjPartition,TwoBackBlock,numClusters);
end	
hcprTP = mean(hcprTP,1);
hcpnTP = nanmean(hcpnTP,1);

% compare HCP and PNC 
onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
offDiag = 1:(numClusters^2); offDiag(onDiag) = [];

load([masterdir,'/analyses/transitionprobabilities/RestCombTransitionProbabilitiesNoPersist_k',num2str(numClusters),name_root,'.mat'])
pncrTP = mean(transitionProbability,1);
TwoBack = load(fullfile(masterdir,'analyses/nbackblocks',['TransProbsNoPersist2back_k',num2str(numClusters),name_root,'.mat']),'BlockTransitionProbability','BlockTransitionProbabilityMats');
pncnTP = nanmean(TwoBack.BlockTransitionProbabilityMats,1);

% correlate TPs and PPs between HCP and PNC
restTPsim = round(corr(pncrTP(offDiag)',hcprTP(offDiag)'),2,'significant');
nbackTPsim = round(corr(pncnTP(offDiag)',hcpnTP(offDiag)'),2,'significant');


%% test rest vs. n-back transitions

restTransitionProbabilityMats = HCPresttransprobs;
nBackTransitionProbabilityMats = HCPnBacktransprobs;

%% permute rest and n-back within subjects
disp('start permutation testing')
nperms = 100000;
pvals_twotail = PERM_TEST(nBackTransitionProbabilityMats,restTransitionProbabilityMats,nperms);

% plot
HCPgrpAvgRest = squeeze(mean(HCPresttransprobs,1)) .* ~eye(numClusters);
HCPgrpAvgnBack = squeeze(nanmean(HCPnBacktransprobs,1)) .* ~eye(numClusters);
maxVal = max(max([HCPgrpAvgRest,HCPgrpAvgnBack]));

f = figure;

subplot(1,3,1);
imagesc(HCPgrpAvgRest);
xticks(1:numClusters); yticks(1:numClusters); colormap('plasma'); 
xticklabels(clusterNames); xtickangle(90); yticklabels(clusterNames); axis square;
COLOR_TICK_LABELS(true,true,numClusters);
ylabel('Current State'); xlabel('Next New State');
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

ylabel('Current State'); xlabel('Next New State');
title({'2-back',['r = ',num2str(nbackTPsim),' with PNC']});
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');
caxis([0 maxVal]); colorbar

subplot(1,3,3);
grpAvgRest = squeeze(mean(restTransitionProbabilityMats,1));    
grpAvgnBack = squeeze(nanmean(nBackTransitionProbabilityMats,1));
HCPnBackMinusRestTP = (grpAvgnBack-grpAvgRest);
imagesc(HCPnBackMinusRestTP.*~eye(numClusters)); colormap('plasma');
xticks(1:numClusters); xticklabels(clusterNames); xtickangle(90);
yticks(1:numClusters); yticklabels(clusterNames); axis square
ylabel('Current State'); xlabel('Next New State');
sig_thresh = 0.05 / numClusters^2;      % bonferroni correction, for two-tailed p-values so only
[y,x] = find(pvals_twotail.*~eye(numClusters) < sig_thresh);
text(x-.12,y+.12,'*','Color','w');
caxis_bound = max(max(abs(HCPnBackMinusRestTP.*~eye(numClusters))));
h = colorbar; ylabel(h,'2-back - rest'); caxis([-caxis_bound caxis_bound]); h.Ticks = [-caxis_bound 0 caxis_bound]; h.TickLabels = [round(-caxis_bound,2,'significant') 0 round(caxis_bound,2,'significant')];
COLOR_TICK_LABELS(true,true,numClusters);
title('n-back > Rest');
set(gca,'FontSize',8);
set(gca,'TickLength',[0 0]);
set(gca,'Fontname','arial');

f.PaperUnits = 'inches';
f.PaperSize = [8 4];
f.PaperPosition = [0 0 8 4];

saveas(f,fullfile(savedir,['HCP_XH_Rest2BackTransProbs_k',num2str(numClusters),'_R',num2str(rTR),'N',num2str(nTR),name_root,'.pdf']),'pdf');
